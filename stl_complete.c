#include "postgres.h"
#include "fmgr.h"
#include "funcapi.h"
#include "utils/array.h"
#include <math.h>

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


// LOESS (tricube kernel, local regression approximation) 

static double loess_local(const double *x, const double *w,
                          int n, int i, int span)
{
    double num = 0.0, den = 0.0;
    int half = span / 2;
    for (int j = i - half; j <= i + half; j++)
    {
        if (j >= 0 && j < n)
        {   double weight = w ? w[j] : 1.0;
            double dist = fabs(j - i);
            double d = dist / (double)half;
            double tricube = (d >= 1.0) ? 0.0 : pow(1.0 - pow(d, 3), 3);
            double wt = weight * tricube;
            num += wt * x[j];
            den += wt;
        }
    }
    return (den > 0.0) ? num / den : x[i];
}

// INITIAL TREND (moving average like STL init)

static void init_trend(const double *y, double *trend, int n)
{
    for (int i = 0; i < n; i++)
        trend[i] = y[i];
}

//SEASONAL COMPONENT (LOESS over seasonal indices)

static void stl_seasonal(const double *y, const double *trend,
                         double *seasonal,
                         int n, int period)
{
    for (int i = 0; i < period; i++)
    {
        double sum = 0.0;
        int count = 0;
        for (int j = i; j < n; j += period)
        {   sum += (y[j] - trend[j]);
            count++;}
        double s = (count > 0) ? sum / count : 0.0;
        for (int j = i; j < n; j += period)
            seasonal[j] = s;
    }
}

// ROBUST WEIGHTS (Cleveland STL idea using MAD approx)
static void stl_update_weights(const double *residual,
                               double *weights,
                               int n)
{
    double mean_abs = 0.0;
    for (int i = 0; i < n; i++) mean_abs += fabs(residual[i]);
    mean_abs /= n;

    for (int i = 0; i < n; i++)
    {  double r = fabs(residual[i]) / (6.0 * mean_abs + 1e-12);
        if (r < 1.0)
            weights[i] = pow(1.0 - r * r, 2);
        else
            weights[i] = 0.0;
    }
}

//MAIN STL FUNCTION (PostgreSQL)
Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);

    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));

    if (n < period || period < 3)
        PG_RETURN_NULL();

    double *y = (double *) ARR_DATA_PTR(arr);

    double *trend = palloc(sizeof(double) * n);
    double *seasonal = palloc(sizeof(double) * n);
    double *residual = palloc(sizeof(double) * n);
    double *weights = palloc(sizeof(double) * n);
    double *detrended = palloc(sizeof(double) * n);

    if (!trend || !seasonal || !residual || !weights || !detrended)
        elog(ERROR, "Memory allocation failed");

    for (int i = 0; i < n; i++)
        weights[i] = 1.0;

    init_trend(y, trend, n);

    int outer_iter = 3;  

    for (int iter = 0; iter < outer_iter; iter++)
    {
        stl_seasonal(y, trend, seasonal, n, period);
        for (int i = 0; i < n; i++)
            detrended[i] = y[i] - seasonal[i];
        int span = period * 2 + 1;
        for (int i = 0; i < n; i++)
            trend[i] = loess_local(detrended, weights, n, i, span);
        for (int i = 0; i < n; i++)
            residual[i] = y[i] - trend[i] - seasonal[i];
        stl_update_weights(residual, weights, n);
    }
    for (int i = 0; i < n; i++)
        residual[i] = y[i] - trend[i] - seasonal[i];
    Datum *result = palloc(sizeof(Datum) * n);
    for (int i = 0; i < n; i++)
        result[i] = Float8GetDatum(trend[i]);
    ArrayType *res_arr = construct_array(
        result, n, FLOAT8OID, sizeof(double), true, 'd'
    );

    PG_RETURN_ARRAYTYPE_P(res_arr);
}

//  STL_old
Datum stl_decompose_old(PG_FUNCTION_ARGS) {
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    
    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));
    if (n < 2) PG_RETURN_NULL();
    
    double *data = (double *) ARR_DATA_PTR(arr);
    Datum *result = palloc(n * sizeof(Datum));
    
    // ЭКСПОНЕНЦИАЛЬНОЕ СГЛАЖИВАНИE
    double alpha = 0.3;
    double trend = data[0];
    
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            trend = alpha * data[i] + (1.0 - alpha) * trend;
        }
        result[i] = Float8GetDatum(trend);
    }
    
    ArrayType *res_arr = construct_array(result, n, FLOAT8OID, 
                                        sizeof(double), true, 'd');
    PG_RETURN_ARRAYTYPE_P(res_arr);
}

void _PG_init(void) {}
