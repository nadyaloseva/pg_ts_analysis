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

    if (fabs(r[0]) < 1e-12) {
        for (int k = 1; k <= max_lag; k++) out[k-1] = NAN; return;}

    sigma[0] = r[0];

    for (int k = 1; k <= max_lag; k++) {

        double sum = 0.0;

        for (int j = 1; j < k; j++) {sum += phi[k-1][j] * r[k - j];}

        phi[k][k] = (r[k] - sum) / sigma[k-1];

        for (int j = 1; j < k; j++) {phi[k][j] = phi[k-1][j] - phi[k][k] * phi[k-1][k-j];}
        sigma[k] = sigma[k-1] * (1.0 - phi[k][k] * phi[k][k]);
        out[k-1] = phi[k][k];
    }
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

//  STL
Datum stl_decompose(PG_FUNCTION_ARGS) {
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
