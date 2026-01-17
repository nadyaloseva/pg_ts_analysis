#include "postgres.h"
#include "fmgr.h"
#include "funcapi.h"
#include "utils/array.h"
#include <math.h>

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(acf_array);
PG_FUNCTION_INFO_V1(pacf_array);

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

// PACF (без изменений)
static double pacf_lag(const double *data, int n, int lag) {
    if (lag >= n || n < 2) return 0.0;
    
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
    int count = 0;
    for (int i = lag; i < n; i++) {
        double x = data[i - lag];
        double y = data[i];
        sum_x += x; sum_y += y; sum_xy += x * y; sum_xx += x * x;
        count++;
    }
    if (count == 0) return 0.0;
    double mean_x = sum_x / count, mean_y = sum_y / count;
    double cov_xy = sum_xy / count - mean_x * mean_y;
    double var_x = sum_xx / count - mean_x * mean_x;
    return (var_x > 1e-10) ? cov_xy / var_x : 0.0;
}

Datum pacf_array(PG_FUNCTION_ARGS) {
    ArrayType *input_array = PG_GETARG_ARRAYTYPE_P(0);
    int32 max_lag = PG_GETARG_INT32(1);
    
    int n = ArrayGetNItems(ARR_NDIM(input_array), ARR_DIMS(input_array));
    if (n < 2 || max_lag < 1 || max_lag >= n) PG_RETURN_NULL();
    
    double *data = (double *) ARR_DATA_PTR(input_array);
    Datum *result = palloc(max_lag * sizeof(Datum));
    
    for (int lag = 1; lag <= max_lag; lag++) {
        result[lag-1] = Float8GetDatum(pacf_lag(data, n, lag));
    }
    ArrayType *pacf_result = construct_array(result, max_lag, FLOAT8OID, sizeof(double), true, 'd');
    PG_RETURN_ARRAYTYPE_P(pacf_result);
}
