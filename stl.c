/*
 * pg_stl — PostgreSQL extension: STL time-series decomposition
 * Совместим с statsmodels.tsa.seasonal.STL (Cleveland 1990).
 *
 * Исправления vs предыдущей версии:
 *   1. Все объявления вынесены в начало блоков (C89/C90 compliance).
 *   2. DirectFunctionCall1(float8, ...) заменён на float8in +
 *      OidFunctionCall1 — корректный способ конвертации в PG 16.
 */

#include "postgres.h"
#include "fmgr.h"
#include "funcapi.h"
#include "utils/array.h"
#include "utils/builtins.h"
#include "catalog/pg_type.h"

#include <math.h>
#include <string.h>

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(acf_array);
PG_FUNCTION_INFO_V1(pacf_array);

/* ================================================================== */
/* ACF                                                                  */
/* ================================================================== */
static double acf_lag(const double *data, int n, int lag)
{
    double mean;
    double var0;
    double cov_lag;
    double d1, d2;
    int    i;

    if (lag >= n || n < 2) return 0.0;

    mean = 0.0;
    for (i = 0; i < n; i++) mean += data[i];
    mean /= n;

    var0 = 0.0;
    cov_lag = 0.0;
    for (i = 0; i < n - lag; i++)
    {
        d1 = data[i]       - mean;
        d2 = data[i + lag] - mean;
        var0    += d1 * d1;
        cov_lag += d1 * d2;
    }
    return var0 > 1e-10 ? cov_lag / var0 : 0.0;
}

Datum acf_array(PG_FUNCTION_ARGS)
{
    ArrayType *input_array;
    int32      max_lag;
    int        n;
    double    *data;
    Datum     *result;
    ArrayType *acf_result;
    int        lag;

    input_array = PG_GETARG_ARRAYTYPE_P(0);
    max_lag     = PG_GETARG_INT32(1);
    n = ArrayGetNItems(ARR_NDIM(input_array), ARR_DIMS(input_array));

    if (n < 2 || max_lag < 1 || max_lag >= n)
        PG_RETURN_NULL();

    data   = (double *) ARR_DATA_PTR(input_array);
    result = (Datum *) palloc(max_lag * sizeof(Datum));

    for (lag = 1; lag <= max_lag; lag++)
        result[lag - 1] = Float8GetDatum(acf_lag(data, n, lag));

    acf_result = construct_array(result, max_lag, FLOAT8OID,
                                 sizeof(double), true, 'd');
    PG_RETURN_ARRAYTYPE_P(acf_result);
}

/* ================================================================== */
/* PACF (Yule-Walker)                                                   */
/* ================================================================== */
static void pacf_yw(const double *x, int n, int max_lag, double *out)
{
    double *r;
    double *phi;   /* max_lag+1 x max_lag+1, row-major */
    double *sigma;
    double  mean;
    double  sum;
    int     i, j, k;

    if (n <= 1) return;

    r     = (double *) palloc((max_lag + 1) * sizeof(double));
    phi   = (double *) palloc((max_lag + 1) * (max_lag + 1) * sizeof(double));
    sigma = (double *) palloc((max_lag + 1) * sizeof(double));

    /* инициализация */
    for (i = 0; i <= max_lag; i++)
    {
        sigma[i] = 0.0;
        for (j = 0; j <= max_lag; j++)
            phi[i * (max_lag + 1) + j] = 0.0;
    }

    /* вычисляем автокорреляции r[0..max_lag] */
    mean = 0.0;
    for (i = 0; i < n; i++) mean += x[i];
    mean /= n;

    for (k = 0; k <= max_lag; k++)
    {
        sum = 0.0;
        for (i = 0; i < n - k; i++)
            sum += (x[i] - mean) * (x[i + k] - mean);
        r[k] = (n > 0) ? sum / n : 0.0;
    }

    if (fabs(r[0]) < 1e-12)
    {
        for (k = 1; k <= max_lag; k++) out[k - 1] = 0.0;
        pfree(r); pfree(phi); pfree(sigma);
        return;
    }

    /* алгоритм Дурбина-Левинсона */
    sigma[0] = r[0];
    for (k = 1; k <= max_lag; k++)
    {
        sum = 0.0;
        for (j = 1; j < k; j++)
            sum += phi[(k-1) * (max_lag + 1) + j] * r[k - j];
        phi[k * (max_lag + 1) + k] = (r[k] - sum) / sigma[k - 1];
        for (j = 1; j < k; j++)
            phi[k * (max_lag + 1) + j] =
                phi[(k-1) * (max_lag + 1) + j]
                - phi[k * (max_lag + 1) + k] * phi[(k-1) * (max_lag + 1) + k - j];
        sigma[k] = sigma[k-1] * (1.0 - phi[k*(max_lag+1)+k] * phi[k*(max_lag+1)+k]);
        out[k - 1] = phi[k * (max_lag + 1) + k];
    }

    pfree(r); pfree(phi); pfree(sigma);
}

Datum pacf_array(PG_FUNCTION_ARGS)
{
    ArrayType *input_array;
    int32      max_lag;
    int        n;
    double    *data;
    Datum     *result;
    double    *tmp;
    ArrayType *res;
    int        i;

    input_array = PG_GETARG_ARRAYTYPE_P(0);
    max_lag     = PG_GETARG_INT32(1);
    n = ArrayGetNItems(ARR_NDIM(input_array), ARR_DIMS(input_array));

    if (n < 2 || max_lag < 1 || max_lag >= n)
        PG_RETURN_NULL();

    data   = (double *) ARR_DATA_PTR(input_array);
    result = (Datum *) palloc(max_lag * sizeof(Datum));
    tmp    = (double *) palloc(max_lag * sizeof(double));

    pacf_yw(data, n, max_lag, tmp);
    for (i = 0; i < max_lag; i++)
        result[i] = Float8GetDatum(tmp[i]);

    res = construct_array(result, max_lag, FLOAT8OID,
                          sizeof(double), true, 'd');
    PG_RETURN_ARRAYTYPE_P(res);
}



/* ================================================================== */
/* Вспомогательные функции                                             */
/* ================================================================== */

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static double tricube(double x)
{
    double t;
    x = fabs(x);
    if (x >= 1.0) return 0.0;
    t = 1.0 - x * x * x;
    return t * t * t;
}

static double bisquare(double u)
{
    double t;
    if (fabs(u) >= 1.0) return 0.0;
    t = 1.0 - u * u;
    return t * t;
}

static double median_of(const double *a, int n)
{
    double *tmp;
    double  m;
    tmp = (double *) palloc(n * sizeof(double));
    memcpy(tmp, a, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_double);
    m = (n & 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
    pfree(tmp);
    return m;
}

static double mad_vec(const double *r, int n)
{
    double  med;
    double *tmp;
    double  m;
    int     i;
    tmp = (double *) palloc(n * sizeof(double));
    med = median_of(r, n);
    for (i = 0; i < n; i++) tmp[i] = fabs(r[i] - med);
    m = median_of(tmp, n);
    pfree(tmp);
    return m + 1e-12;
}

static int next_odd_gt(double x)
{
    int v = (int) floor(x) + 1;
    if (v % 2 == 0) v++;
    return v;
}

/* ================================================================== */
/* LOESS в одной точке x0                                             */
/*   Предиктор центрирован: dx = xi - x0 (числовая стабильность)     */
/* ================================================================== */
static double loess_at(double x0,
                       const double *xd, const double *yd, int nd,
                       int bandwidth, int deg,
                       const double *wrob)
{
    double *dist;
    double *dtmp;
    double  dmax;
    double  sw, sw_dx, sw_dx2, sw_dxy, sw_y;
    double  denom, b, a;
    int     j;

    dist = (double *) palloc(nd * sizeof(double));
    dtmp = (double *) palloc(nd * sizeof(double));

    for (j = 0; j < nd; j++)
        dist[j] = fabs(xd[j] - x0);

    memcpy(dtmp, dist, nd * sizeof(double));
    qsort(dtmp, nd, sizeof(double), cmp_double);

    dmax = dtmp[bandwidth - 1];
    if (dmax < 1e-12) dmax = 1e-12;

    sw = sw_dx = sw_dx2 = sw_dxy = sw_y = 0.0;
    for (j = 0; j < nd; j++)
    {
        double w  = tricube(dist[j] / dmax);
        double dx = xd[j] - x0;
        if (wrob) w *= wrob[j];
        sw     += w;
        sw_dx  += w * dx;
        sw_dx2 += w * dx * dx;
        sw_dxy += w * dx * yd[j];
        sw_y   += w * yd[j];
    }

    pfree(dist);
    pfree(dtmp);

    if (sw < 1e-12) return 0.0;
    if (deg == 0)   return sw_y / sw;

    denom = sw * sw_dx2 - sw_dx * sw_dx;
    if (fabs(denom) > 1e-12)
    {
        b = (sw * sw_dxy - sw_dx * sw_y) / denom;
        a = (sw_y - b * sw_dx) / sw;
        return a;   /* при x = x0: a + b*0 = a */
    }
    return sw_y / sw;
}

/* ================================================================== */
/* LOESS по массиву запросных точек с jump-интерполяцией              */
/*   Последняя точка nq-1 всегда вычисляется точно.                  */
/* ================================================================== */
static void loess_fit(const double *xs, int nq,
                      const double *xd, const double *yd, int nd,
                      int bandwidth, int deg,
                      const double *wrob, int jump,
                      double *out)
{
    int    *sidx;
    double *sval;
    int     max_s;
    int     ns, seg, i, k;
    int     bw;

    bw = (bandwidth < nd) ? bandwidth : nd;
    if (bw < 1)   bw   = 1;
    if (jump < 1) jump  = 1;

    max_s = nq / jump + 2;
    sidx  = (int *)    palloc(max_s * sizeof(int));
    sval  = (double *) palloc(max_s * sizeof(double));

    ns = 0;
    for (i = 0; i < nq; i += jump)
        sidx[ns++] = i;
    if (ns == 0 || sidx[ns - 1] != nq - 1)
        sidx[ns++] = nq - 1;

    for (k = 0; k < ns; k++)
        sval[k] = loess_at(xs[sidx[k]], xd, yd, nd, bw, deg, wrob);

    seg = 0;
    for (i = 0; i < nq; i++)
    {
        while (seg < ns - 1 && i >= sidx[seg + 1]) seg++;
        if (i == sidx[seg])
        {
            out[i] = sval[seg];
        }
        else if (seg + 1 < ns)
        {
            double t = (double)(i - sidx[seg]) / (sidx[seg + 1] - sidx[seg]);
            out[i] = sval[seg] + t * (sval[seg + 1] - sval[seg]);
        }
        else
        {
            out[i] = sval[seg];
        }
    }

    pfree(sidx);
    pfree(sval);
}

/* ================================================================== */
/* Low-pass seasonal фильтр: три прохода LOESS                        */
/*   pass 1: LOESS(n_l)                                               */
/*   pass 2: LOESS(n_l)                                               */
/*   pass 3: LOESS(3)                                                 */
/* ================================================================== */
static void low_pass_seasonal(const double *cv, int n,
                               int n_l, int deg_l, int jump_l,
                               double *lp)
{
    double *tmp;
    double *xs;
    int     i;

    tmp = (double *) palloc(n * sizeof(double));
    xs  = (double *) palloc(n * sizeof(double));
    for (i = 0; i < n; i++) xs[i] = (double) i;

    loess_fit(xs, n, xs, cv,  n, n_l, deg_l, NULL, jump_l, tmp);
    loess_fit(xs, n, xs, tmp, n, n_l, deg_l, NULL, jump_l, lp);
    loess_fit(xs, n, xs, lp,  n, 3,   deg_l, NULL, jump_l, tmp);

    memcpy(lp, tmp, n * sizeof(double));
    pfree(tmp);
    pfree(xs);
}

/* ================================================================== */
/* Сезонный шаг                                                        */
/* ================================================================== */
static void seasonal_step(const double *detr, int n, int period,
                           int n_s, int n_l,
                           int deg_s, int deg_l,
                           int jump_s, int jump_l,
                           const double *wrob,
                           double *season)
{
    double *cv;
    double *lp;
    double *xs;
    int     p, i;

    cv = (double *) palloc0(n * sizeof(double));

    for (p = 0; p < period; p++)
    {
        double *xsub;
        double *ysub;
        double *wsub;
        double *yhat;
        int     m, bw, idx;

        m = 0;
        for (i = p; i < n; i += period) m++;

        xsub = (double *) palloc(m * sizeof(double));
        ysub = (double *) palloc(m * sizeof(double));
        wsub = wrob ? (double *) palloc(m * sizeof(double)) : NULL;
        yhat = (double *) palloc(m * sizeof(double));

        idx = 0;
        for (i = p; i < n; i += period)
        {
            xsub[idx] = (double) idx;
            ysub[idx] = detr[i];
            if (wsub) wsub[idx] = wrob[i];
            idx++;
        }

        bw = (n_s < m) ? n_s : m;
        loess_fit(xsub, m, xsub, ysub, m, bw, deg_s, wsub, jump_s, yhat);

        idx = 0;
        for (i = p; i < n; i += period)
            cv[i] = yhat[idx++];

        pfree(xsub); pfree(ysub); pfree(yhat);
        if (wsub) pfree(wsub);
    }

    lp = (double *) palloc(n * sizeof(double));
    xs = (double *) palloc(n * sizeof(double));
    for (i = 0; i < n; i++) xs[i] = (double) i;

    low_pass_seasonal(cv, n, n_l, deg_l, jump_l, lp);

    for (i = 0; i < n; i++)
        season[i] = cv[i] - lp[i];

    pfree(cv);
    pfree(lp);
    pfree(xs);
}

/* ================================================================== */
/* Шаг тренда                                                          */
/* ================================================================== */
static void trend_step(const double *desea, int n,
                        int n_t, int deg_t, int jump_t,
                        const double *wrob, double *trend)
{
    double *xs;
    int     i;
    xs = (double *) palloc(n * sizeof(double));
    for (i = 0; i < n; i++) xs[i] = (double) i;
    loess_fit(xs, n, xs, desea, n, n_t, deg_t, wrob, jump_t, trend);
    pfree(xs);
}

/* ================================================================== */
/* Робастные веса (bisquare)                                           */
/* ================================================================== */
static void compute_robust_weights(const double *y,
                                    const double *trend,
                                    const double *season,
                                    int n, double *rw)
{
    double *r;
    double  mad, c;
    int     i;
    r = (double *) palloc(n * sizeof(double));
    for (i = 0; i < n; i++) r[i] = y[i] - trend[i] - season[i];
    mad = mad_vec(r, n);
    c   = 6.0 * mad;
    for (i = 0; i < n; i++) rw[i] = bisquare(r[i] / c);
    pfree(r);
}

/* ================================================================== */
/* Конфигурация STL                                                    */
/* ================================================================== */
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

static void stl_resolve(STLConfig *c)
{
    double t;
    if (c->trend == 0)
    {
        t = 1.5 * c->period / (1.0 - 1.5 / (double) c->seasonal);
        c->trend = next_odd_gt(t);
    }
    if (c->low_pass == 0)
        c->low_pass = next_odd_gt((double) c->period);
    if (c->robust && c->outer_iter == 0)
        c->outer_iter = 15;
}

/* ================================================================== */
/* Ядро декомпозиции                                                   */
/* ================================================================== */
static void stl_core(const double *y, int n,
                      STLConfig *cfg,
                      double *trend, double *season, double *residual)
{
    double      *detr;
    double      *desea;
    double      *rw;
    int          n_outer;
    int          o, inner, i;
    const double *wrob;

    stl_resolve(cfg);

    detr  = (double *) palloc(n * sizeof(double));
    desea = (double *) palloc(n * sizeof(double));
    rw    = (double *) palloc(n * sizeof(double));

    for (i = 0; i < n; i++) { trend[i] = 0.0; season[i] = 0.0; rw[i] = 1.0; }

    n_outer = cfg->robust ? cfg->outer_iter + 1 : 1;

    for (o = 0; o < n_outer; o++)
    {
        wrob = (o == 0) ? NULL : rw;

        for (inner = 0; inner < cfg->inner_iter; inner++)
        {
            for (i = 0; i < n; i++) detr[i] = y[i] - trend[i];

            seasonal_step(detr, n, cfg->period,
                          cfg->seasonal, cfg->low_pass,
                          cfg->seasonal_deg, cfg->low_pass_deg,
                          cfg->seasonal_jump, cfg->low_pass_jump,
                          wrob, season);

            for (i = 0; i < n; i++) desea[i] = y[i] - season[i];

            trend_step(desea, n,
                       cfg->trend, cfg->trend_deg, cfg->trend_jump,
                       wrob, trend);
        }

        if (cfg->robust && o < n_outer - 1)
            compute_robust_weights(y, trend, season, n, rw);
    }

    for (i = 0; i < n; i++)
        residual[i] = y[i] - trend[i] - season[i];

    pfree(detr);
    pfree(desea);
    pfree(rw);
}

/* ================================================================== */
/* Упаковка массива double[] в PostgreSQL ArrayType                    */
/* ================================================================== */
static ArrayType *double_array_to_pg(const double *data, int n)
{
    Datum *elems;
    int    i;
    elems = (Datum *) palloc(n * sizeof(Datum));
    for (i = 0; i < n; i++)
        elems[i] = Float8GetDatum(data[i]);
    return construct_array(elems, n, FLOAT8OID,
                           sizeof(double), FLOAT8PASSBYVAL, 'd');
}

/* ================================================================== */
/* Datum: stl_decompose                                                */
/* ================================================================== */
PG_FUNCTION_INFO_V1(stl_decompose);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType  *arr;
    int32       period;
    int32       seasonal;
    bool        robust;
    int32       trend_w;
    int32       low_pass_w;
    int32       inner_iter;
    int32       outer_iter;
    int         n;
    double     *y;
    STLConfig   cfg;
    double     *trend_out;
    double     *seasonal_out;
    double     *residual_out;
    TupleDesc   tupdesc;
    Datum       values[3];
    bool        isnull[3];
    HeapTuple   tuple;

    arr        = PG_GETARG_ARRAYTYPE_P(0);
    period     = PG_GETARG_INT32(1);
    seasonal   = PG_ARGISNULL(2) ? 7    : PG_GETARG_INT32(2);
    robust     = PG_ARGISNULL(3) ? true : PG_GETARG_BOOL(3);
    trend_w    = PG_ARGISNULL(4) ? 0    : PG_GETARG_INT32(4);
    low_pass_w = PG_ARGISNULL(5) ? 0    : PG_GETARG_INT32(5);
    inner_iter = PG_ARGISNULL(6) ? 2    : PG_GETARG_INT32(6);
    outer_iter = PG_ARGISNULL(7) ? 0    : PG_GETARG_INT32(7);

    if (ARR_NDIM(arr) != 1)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("stl_decompose: input must be a 1-D array")));

    n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));

    if (period < 2)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("stl_decompose: period must be >= 2")));

    if (n < 2 * period)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("stl_decompose: series length (%d) must be >= 2 * period (%d)",
                        n, 2 * period)));

    if (seasonal < 3 || seasonal % 2 == 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("stl_decompose: seasonal must be an odd integer >= 3")));

    if (ARR_HASNULL(arr))
        ereport(ERROR,
                (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
                 errmsg("stl_decompose: input array must not contain NULLs")));

    if (inner_iter < 1) inner_iter = 1;
    if (outer_iter < 0) outer_iter = 0;

    /* Принимаем только float8[], остальные типы отклоняем с подсказкой */
    if (ARR_ELEMTYPE(arr) != FLOAT8OID)
        ereport(ERROR,
                (errcode(ERRCODE_DATATYPE_MISMATCH),
                 errmsg("stl_decompose: input array must be double precision[]"),
                 errhint("Cast your array: ARRAY[...]::double precision[]")));

    y = (double *) ARR_DATA_PTR(arr);

    cfg.period        = period;
    cfg.seasonal      = seasonal;
    cfg.trend         = trend_w;
    cfg.low_pass      = low_pass_w;
    cfg.seasonal_deg  = 1;
    cfg.trend_deg     = 1;
    cfg.low_pass_deg  = 1;
    cfg.robust        = robust ? 1 : 0;
    cfg.seasonal_jump = 1;
    cfg.trend_jump    = 1;
    cfg.low_pass_jump = 1;
    cfg.inner_iter    = inner_iter;
    cfg.outer_iter    = outer_iter;

    trend_out    = (double *) palloc(n * sizeof(double));
    seasonal_out = (double *) palloc(n * sizeof(double));
    residual_out = (double *) palloc(n * sizeof(double));

    stl_core(y, n, &cfg, trend_out, seasonal_out, residual_out);

    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("stl_decompose: must be called as a record-returning function")));

    tupdesc = BlessTupleDesc(tupdesc);

    isnull[0] = isnull[1] = isnull[2] = false;
    values[0] = PointerGetDatum(double_array_to_pg(trend_out,    n));
    values[1] = PointerGetDatum(double_array_to_pg(seasonal_out, n));
    values[2] = PointerGetDatum(double_array_to_pg(residual_out, n));

    tuple = heap_form_tuple(tupdesc, values, isnull);
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/* ================================================================== */
/* Общий хелпер для stl_trend / stl_seasonal / stl_residual           */
/* ================================================================== */
static ArrayType *stl_component(PG_FUNCTION_ARGS, int component)
{
    ArrayType  *arr;
    int32       period;
    int32       seasonal;
    bool        robust;
    int32       trend_w;
    int32       low_pass_w;
    int32       inner_iter;
    int32       outer_iter;
    int         n;
    double     *y;
    STLConfig   cfg;
    double     *trend_out;
    double     *seasonal_out;
    double     *residual_out;

    arr        = PG_GETARG_ARRAYTYPE_P(0);
    period     = PG_GETARG_INT32(1);
    seasonal   = PG_ARGISNULL(2) ? 7    : PG_GETARG_INT32(2);
    robust     = PG_ARGISNULL(3) ? true : PG_GETARG_BOOL(3);
    trend_w    = PG_ARGISNULL(4) ? 0    : PG_GETARG_INT32(4);
    low_pass_w = PG_ARGISNULL(5) ? 0    : PG_GETARG_INT32(5);
    inner_iter = PG_ARGISNULL(6) ? 2    : PG_GETARG_INT32(6);
    outer_iter = PG_ARGISNULL(7) ? 0    : PG_GETARG_INT32(7);

    if (ARR_NDIM(arr) != 1)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("input must be a 1-D array")));
    if (ARR_HASNULL(arr))
        ereport(ERROR,
                (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
                 errmsg("input array must not contain NULLs")));
    if (ARR_ELEMTYPE(arr) != FLOAT8OID)
        ereport(ERROR,
                (errcode(ERRCODE_DATATYPE_MISMATCH),
                 errmsg("input array must be double precision[]"),
                 errhint("Cast: ARRAY[...]::double precision[]")));

    n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));

    if (period < 2 || n < 2 * period)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("invalid period or series too short")));

    y = (double *) ARR_DATA_PTR(arr);

    cfg.period        = period;
    cfg.seasonal      = (seasonal < 3) ? 7 : seasonal;
    cfg.trend         = trend_w;
    cfg.low_pass      = low_pass_w;
    cfg.seasonal_deg  = 1;
    cfg.trend_deg     = 1;
    cfg.low_pass_deg  = 1;
    cfg.robust        = robust ? 1 : 0;
    cfg.seasonal_jump = 1;
    cfg.trend_jump    = 1;
    cfg.low_pass_jump = 1;
    cfg.inner_iter    = (inner_iter < 1) ? 1 : inner_iter;
    cfg.outer_iter    = (outer_iter < 0) ? 0 : outer_iter;

    trend_out    = (double *) palloc(n * sizeof(double));
    seasonal_out = (double *) palloc(n * sizeof(double));
    residual_out = (double *) palloc(n * sizeof(double));

    stl_core(y, n, &cfg, trend_out, seasonal_out, residual_out);

    switch (component)
    {
        case 0:  return double_array_to_pg(trend_out,    n);
        case 1:  return double_array_to_pg(seasonal_out, n);
        default: return double_array_to_pg(residual_out, n);
    }
}

PG_FUNCTION_INFO_V1(stl_trend);
Datum stl_trend(PG_FUNCTION_ARGS)
{
    PG_RETURN_ARRAYTYPE_P(stl_component(fcinfo, 0));
}

PG_FUNCTION_INFO_V1(stl_seasonal);
Datum stl_seasonal(PG_FUNCTION_ARGS)
{
    PG_RETURN_ARRAYTYPE_P(stl_component(fcinfo, 1));
}

PG_FUNCTION_INFO_V1(stl_residual);
Datum stl_residual(PG_FUNCTION_ARGS)
{
    PG_RETURN_ARRAYTYPE_P(stl_component(fcinfo, 2));
}