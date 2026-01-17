#include "postgres.h"
#include "fmgr.h"
#include "utils/array.h"
#include "catalog/pg_type.h"
#include "math.h"

PG_FUNCTION_INFO_V1(stl_decompose);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));
    double *data = (double *) ARR_DATA_PTR(arr);
    
    Datum *result = palloc(n * 3 * sizeof(Datum));
    
    double trend_sum = 0;
    for (int i = 0; i < n; i++) {
        trend_sum += data[i];
        double trend_val = trend_sum / (i + 1.0);
        double seasonal_val = sin(6.28318 * i / period) * 0.1;
        double remainder_val = data[i] - trend_val - seasonal_val;
        
        result[i*3]     = Float8GetDatum(trend_val);
        result[i*3 + 1] = Float8GetDatum(seasonal_val);
        result[i*3 + 2] = Float8GetDatum(remainder_val);
    }
    
    ArrayType *res_arr = construct_array(result, n*3, FLOAT8OID,
                                         sizeof(double), true, 'i');
    PG_RETURN_ARRAYTYPE_P(res_arr);
}
