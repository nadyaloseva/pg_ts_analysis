#include "postgres.h"
#include "fmgr.h"
#include "funcapi.h"
#include "utils/array.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

PG_MODULE_MAGIC;

// ДЕКЛАРАЦИИ ФУНКЦИЙ
PG_FUNCTION_INFO_V1(acf_array);
PG_FUNCTION_INFO_V1(pacf_array);
PG_FUNCTION_INFO_V1(stl_decompose);

// ACF
static double acf_lag(const double *data, int n, int lag) {
    if (lag >= n || n < 2) return 0.0;
    
    double mean = 0;
    for (int i = 0; i < n; i++) mean += data[i];
    mean /= n;
    
    double var0 = 0, cov_lag = 0;
    for (int i = 0; i < n - lag; i++) {
        double d1 = data[i] - mean;
        double d2 = data[i + lag] - mean;
        var0 += d1 * d1;
        cov_lag += d1 * d2;
    }
    return var0 > 1e-10 ? cov_lag / var0 : 0.0;
}

Datum acf_array(PG_FUNCTION_ARGS) {
    ArrayType *input_array = PG_GETARG_ARRAYTYPE_P(0);
    int32 max_lag = PG_GETARG_INT32(1);
    
    int n = ArrayGetNItems(ARR_NDIM(input_array), ARR_DIMS(input_array));
    if (n < 2 || max_lag < 1 || max_lag >= n) PG_RETURN_NULL();
    
    double *data = (double *) ARR_DATA_PTR(input_array);
    Datum *result = palloc(max_lag * sizeof(Datum));
    
    for (int lag = 1; lag <= max_lag; lag++) {
        result[lag-1] = Float8GetDatum(acf_lag(data, n, lag));
    }
    ArrayType *acf_result = construct_array(result, max_lag, FLOAT8OID, sizeof(double), true, 'd');
    PG_RETURN_ARRAYTYPE_P(acf_result);
}

//  PACF 
static void pacf_yw(const double *x, int n, int max_lag, double *out) {
    if (n <= 1) return;
    double r[max_lag + 1];
    double phi[max_lag + 1][max_lag + 1];
    double sigma[max_lag + 1];

    for (int i = 0; i <= max_lag; i++) { 
        sigma[i] = 0.0;
        for (int j = 0; j <= max_lag; j++)
            phi[i][j] = 0.0;}

    for (int k = 0; k <= max_lag; k++) {
        double mean = 0.0;
        for (int i = 0; i < n; i++)
            mean += x[i];
        mean /= n;
        double sum = 0.0;
        for (int i = 0; i < n - k; i++) {sum += (x[i] - mean) * (x[i + k] - mean);}
        r[k] = (n > 0) ? sum / n : 0.0;
    }
    if (fabs(r[0]) < 1e-12) {for (int k = 1; k <= max_lag; k++) out[k-1] = NAN; return;}
    sigma[0] = r[0];
    for (int k = 1; k <= max_lag; k++) {
        double sum = 0.0;
        for (int j = 1; j < k; j++) {sum += phi[k-1][j] * r[k - j];}
        phi[k][k] = (r[k] - sum) / sigma[k-1];
        for (int j = 1; j < k; j++) {phi[k][j] = phi[k-1][j] - phi[k][k] * phi[k-1][k-j];}
        sigma[k] = sigma[k-1] * (1.0 - phi[k][k] * phi[k][k]);
        out[k-1] = phi[k][k];}
}

Datum pacf_array(PG_FUNCTION_ARGS) {
    ArrayType *input_array = PG_GETARG_ARRAYTYPE_P(0);
    int32 max_lag = PG_GETARG_INT32(1);
    int n = ArrayGetNItems(ARR_NDIM(input_array), ARR_DIMS(input_array));

    if (n < 2 || max_lag < 1 || max_lag >= n)
        PG_RETURN_NULL();

    double *data = (double *) ARR_DATA_PTR(input_array);
    Datum *result = palloc(max_lag * sizeof(Datum));
    double *tmp = palloc(max_lag * sizeof(double));

    pacf_yw(data, n, max_lag, tmp);
    for (int i = 0; i < max_lag; i++) {
        result[i] = Float8GetDatum(tmp[i]);
    }
    ArrayType *res = construct_array(result, max_lag, FLOAT8OID, sizeof(double), true,'d');
    PG_RETURN_ARRAYTYPE_P(res);
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}
 
static double tricube(double x) {
    x = fabs(x);
    if (x >= 1.0) return 0.0;
    double t = 1.0 - x * x * x;
    return t * t * t;
}
 
static double bisquare(double u) {
    double au = fabs(u);
    if (au >= 1.0) return 0.0;
    double t = 1.0 - u * u;
    return t * t;
}
 
static double median_of(const double *a, int n) {
    double *tmp = malloc(n * sizeof(double));
    memcpy(tmp, a, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_double);
    double m = (n & 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    free(tmp);
    return m;
}
 
static double mad_vec(const double *r, int n) {
    double med = median_of(r, n);
    double *tmp = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) tmp[i] = fabs(r[i] - med);
    double m = median_of(tmp, n);
    free(tmp);
    return m + 1e-12;
}
 
static int next_odd_gt(double x) {
    int v = (int)ceil(x);
    /* "strictly greater" */
    if ((double)v <= x) v++;
    if (v % 2 == 0) v++;
    return v;
}

static void loess_fit(
    const double *xs, int nq,
    const double *xd, const double *yd, int nd,
    int bandwidth, int deg,
    const double *wrob, int jump,
    double *out)
{
    if (bandwidth < 1) bandwidth = 1;
    if (bandwidth > nd) bandwidth = nd;
    if (jump < 1) jump = 1;
 
    double *dist = malloc(nd * sizeof(double));
    double *dtmp = malloc(nd * sizeof(double));
 
    int max_samples = nq / jump + 2;
    int   *sidx = malloc(max_samples * sizeof(int));
    double *sval = malloc(max_samples * sizeof(double));
    int ns = 0; /* ns >= 1 guaranteed since nq >= 1 */
    for (int i = 0; i < nq; i += jump)
        sidx[ns++] = i;
    if (sidx[ns - 1] != nq - 1)
        sidx[ns++] = nq - 1; 
 
    for (int k = 0; k < ns; k++) {
        double x0 = xs[sidx[k]];
 
        for (int j = 0; j < nd; j++)
            dist[j] = fabs(xd[j] - x0);
 
        memcpy(dtmp, dist, nd * sizeof(double));
        qsort(dtmp, nd, sizeof(double), cmp_double);
        double dmax = dtmp[bandwidth - 1];
        if (dmax < 1e-12) dmax = 1e-12;
 
        double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
        for (int j = 0; j < nd; j++) {
            double w = tricube(dist[j] / dmax);
            if (wrob) w *= wrob[j];
            sw   += w;
            swx  += w * xd[j];
            swy  += w * yd[j];
            swxx += w * xd[j] * xd[j];
            swxy += w * xd[j] * yd[j];
        }
 
        double val;
        if (sw < 1e-12) {
            val = 0.0;
        } else if (deg == 0) {
            val = swy / sw;
        } else {
            double denom = sw * swxx - swx * swx;
            if (fabs(denom) > 1e-12) {
                double b = (sw * swxy - swx * swy) / denom;
                double a = (swy - b * swx) / sw;
                val = a + b * x0;
            } else {
                val = swy / sw; 
            }
        }
        sval[k] = val;
    }
 
    int seg = 0;
    for (int i = 0; i < nq; i++) {
        while (seg < ns - 1 && i >= sidx[seg + 1])
            seg++;
 
        if (i == sidx[seg]) {
            out[i] = sval[seg];
        } else if (seg + 1 < ns) {
            double lo_x = xs[sidx[seg]];
            double hi_x = xs[sidx[seg + 1]];
            double t = (xs[i] - lo_x) / (hi_x - lo_x);
            out[i] = sval[seg] + t * (sval[seg + 1] - sval[seg]);
        } else {
            out[i] = sval[seg];
        }
    }
 
    free(dist); free(dtmp); free(sidx); free(sval);
}
 
static void seasonal_step(
    const double *detr, int n, int period,
    int n_s, int n_l,
    int deg_s, int deg_l,
    int jump_s, int jump_l,
    const double *wrob,
    double *season)
{
    double *cv = calloc(n, sizeof(double));
 
    for (int p = 0; p < period; p++) {
        
        int m = 0;
        for (int i = p; i < n; i += period) m++;
 
        double *xsub = malloc(m * sizeof(double));
        double *ysub = malloc(m * sizeof(double));
        double *wsub = wrob ? malloc(m * sizeof(double)) : NULL;
        double *yhat = malloc(m * sizeof(double));
 
        int idx = 0;
        for (int i = p; i < n; i += period) {
            xsub[idx] = (double)idx;
            ysub[idx] = detr[i];
            if (wsub) wsub[idx] = wrob[i];
            idx++;
        }
 
        int bw = (n_s < m) ? n_s : m;
        loess_fit(xsub, m, xsub, ysub, m, bw, deg_s, wsub, jump_s, yhat);
 
        idx = 0;
        for (int i = p; i < n; i += period)
            cv[i] = yhat[idx++];
 
        free(xsub); free(ysub); free(yhat);
        if (wsub) free(wsub);
    }
 
    
    double *lp = malloc(n * sizeof(double));
    double *xs = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) xs[i] = (double)i;
 
    loess_fit(xs, n, xs, cv, n, n_l, deg_l, NULL, jump_l, lp);
 
    for (int i = 0; i < n; i++)
        season[i] = cv[i] - lp[i];
 
    free(cv); free(lp); free(xs);
}
 

static void trend_step(
    const double *desea, int n,
    int n_t, int deg_t, int jump_t,
    const double *wrob, double *trend)
{
    double *xs = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) xs[i] = (double)i;
    loess_fit(xs, n, xs, desea, n, n_t, deg_t, wrob, jump_t, trend);
    free(xs);
}
 

static void compute_robust_weights(
    const double *y, const double *trend, const double *season,
    int n, double *rw)
{
    double *r = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) r[i] = y[i] - trend[i] - season[i];
    double mad = mad_vec(r, n);
    double c   = 6.0 * mad;
    for (int i = 0; i < n; i++) rw[i] = bisquare(r[i] / c);
    free(r);
}

typedef struct {
    int period;
    int seasonal;      
    int trend;         
    int low_pass;     
    int seasonal_deg;
    int trend_deg;
    int low_pass_deg;
    int robust;
    int seasonal_jump;
    int trend_jump;
    int low_pass_jump;
    int inner_iter;   
    int outer_iter;   
} STLConfig;
 
static void stl_resolve_defaults(STLConfig *c) {
    if (c->trend == 0) {
        /* smallest odd > 1.5*period / (1 - 1.5/seasonal) */
        double t = 1.5 * c->period / (1.0 - 1.5 / (double)c->seasonal);
        c->trend = next_odd_gt(t);
    }
    if (c->low_pass == 0) {
        /* smallest odd > period */
        c->low_pass = next_odd_gt((double)c->period);
    }
    if (c->outer_iter == 0)
        c->outer_iter = c->robust ? 15 : 0;
}
 

void stl_fit(
    const double *y, int n,
    STLConfig *cfg,
    double *trend, double *season, double *residual)
{
    stl_resolve_defaults(cfg);
 
    double *detr  = malloc(n * sizeof(double));
    double *desea = malloc(n * sizeof(double));
    double *rw    = malloc(n * sizeof(double));
 
    for (int i = 0; i < n; i++) { trend[i] = 0.0; season[i] = 0.0; rw[i] = 1.0; }

    int n_outer = cfg->robust ? cfg->outer_iter + 1 : 1;
 
    for (int o = 0; o < n_outer; o++) {
        const double *wrob = (o == 0) ? NULL : rw;
 
        for (int inner = 0; inner < cfg->inner_iter; inner++) {
            for (int i = 0; i < n; i++) detr[i] = y[i] - trend[i];
 
            seasonal_step(detr, n, cfg->period,
                          cfg->seasonal, cfg->low_pass,
                          cfg->seasonal_deg, cfg->low_pass_deg,
                          cfg->seasonal_jump, cfg->low_pass_jump,
                          wrob, season);
 
            for (int i = 0; i < n; i++) desea[i] = y[i] - season[i];
 
            trend_step(desea, n,
                       cfg->trend, cfg->trend_deg, cfg->trend_jump,
                       wrob, trend);
        }
 
        if (cfg->robust && o < n_outer - 1)
            compute_robust_weights(y, trend, season, n, rw);
    }
 
    for (int i = 0; i < n; i++)
        residual[i] = y[i] - trend[i] - season[i];
 
    free(detr); free(desea); free(rw);
}

PG_MODULE_MAGIC;

typedef struct {
    int period;
    int seasonal;
    int trend;
    int low_pass;
    int seasonal_deg;
    int trend_deg;
    int low_pass_deg;
    int robust;
    int seasonal_jump;
    int trend_jump;
    int low_pass_jump;
    int inner_iter;
    int outer_iter;
} STLConfig;

extern void stl_fit(
    const double *y, int n,
    STLConfig *cfg,
    double *trend,
    double *season,
    double *residual
);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);

    int period        = PG_GETARG_INT32(1);
    int seasonal      = PG_GETARG_INT32(2);
    int trend         = PG_GETARG_INT32(3);
    int low_pass      = PG_GETARG_INT32(4);
    int seasonal_deg  = PG_GETARG_INT32(5);
    int trend_deg     = PG_GETARG_INT32(6);
    int low_pass_deg  = PG_GETARG_INT32(7);
    int robust        = PG_GETARG_INT32(8);
    int seasonal_jump = PG_GETARG_INT32(9);
    int trend_jump    = PG_GETARG_INT32(10);
    int low_pass_jump = PG_GETARG_INT32(11);
    int inner_iter    = PG_GETARG_INT32(12);
    int outer_iter    = PG_GETARG_INT32(13);

    int n;
    Datum *values;
    bool *nulls;
    int nelems;

    deconstruct_array(arr,
                      FLOAT8OID,
                      sizeof(double),
                      true,
                      'd',
                      &values,
                      &nulls,
                      &nelems);

    n = nelems;

    if (n <= 0 || n < period)
        PG_RETURN_NULL();

    double *y = (double *) palloc(sizeof(double) * n);

    for (int i = 0; i < n; i++)
    {
        if (nulls[i])
            y[i] = 0.0;
        else
            y[i] = DatumGetFloat8(values[i]);
    }

    STLConfig cfg;
    cfg.period        = period;
    cfg.seasonal      = seasonal;
    cfg.trend         = trend;
    cfg.low_pass      = low_pass;
    cfg.seasonal_deg  = seasonal_deg;
    cfg.trend_deg     = trend_deg;
    cfg.low_pass_deg  = low_pass_deg;
    cfg.robust        = robust;
    cfg.seasonal_jump = seasonal_jump;
    cfg.trend_jump    = trend_jump;
    cfg.low_pass_jump = low_pass_jump;
    cfg.inner_iter    = inner_iter;
    cfg.outer_iter    = outer_iter;

    double *trend_out    = (double *) palloc(sizeof(double) * n);
    double *season_out   = (double *) palloc(sizeof(double) * n);
    double *resid_out    = (double *) palloc(sizeof(double) * n);

    if (!trend_out || !season_out || !resid_out)
        elog(ERROR, "Memory allocation failed");

    stl_fit(y, n, &cfg, trend_out, season_out, resid_out);

    Datum *trend_datums = (Datum *) palloc(sizeof(Datum) * n);
    for (int i = 0; i < n; i++)
        trend_datums[i] = Float8GetDatum(trend_out[i]);

    ArrayType *trend_arr =
        construct_array(trend_datums, n, FLOAT8OID, sizeof(double), true, 'd');

    Datum *season_datums = (Datum *) palloc(sizeof(Datum) * n);
    for (int i = 0; i < n; i++)
        season_datums[i] = Float8GetDatum(season_out[i]);

    ArrayType *season_arr =
        construct_array(season_datums, n, FLOAT8OID, sizeof(double), true, 'd');

    Datum *res_datums = (Datum *) palloc(sizeof(Datum) * n);
    for (int i = 0; i < n; i++)
        res_datums[i] = Float8GetDatum(resid_out[i]);

    ArrayType *res_arr =
        construct_array(res_datums, n, FLOAT8OID, sizeof(double), true, 'd');

    TupleDesc tupdesc;
    HeapTuple tuple;
    Datum result_values[3];
    bool result_nulls[3] = {false, false, false};

    tupdesc = RelationNameGetTupleDesc("stl_result");

    result_values[0] = PointerGetDatum(trend_arr);
    result_values[1] = PointerGetDatum(season_arr);
    result_values[2] = PointerGetDatum(res_arr);

    tuple = heap_form_tuple(tupdesc, result_values, result_nulls);

    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

void _PG_init(void) {}
