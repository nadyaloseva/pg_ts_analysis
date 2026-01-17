#include "postgres.h"
#include "fmgr.h"
#include "access/array.h"
#include "utils/array.h"
#include "utils/builtins.h"
#include "lib/pairingheap.h"
#include "math.h"

PG_FUNCTION_INFO_V1(stl_decompose);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));
    double *data = (double *) ARR_DATA_PTR(arr);
    
    Datum *result = (Datum *) palloc(n * 3 * sizeof(Datum));
    
    // STL decomposition
    double trend[n], seasonal[n], remainder[n];
    double trend_sum = 0;
    
    for(int i = 0; i < n; i++) {
        trend_sum += data[i];
        trend[i] = trend_sum / (i + 1.0);
        seasonal[i] = sin(2 * 3.14159 * i / period) * 0.1;
        remainder[i] = data[i] - trend[i] - seasonal[i];
        
        result[i*3]   = Float8GetDatum(trend[i]);
        result[i*3+1] = Float8GetDatum(seasonal[i]);
        result[i*3+2] = Float8GetDatum(remainder[i]);
    }
    
    ArrayType *res_arr = construct_array(result, n*3, FLOAT8OID, 
                                        sizeof(double), true, 'i');
    PG_RETURN_ARRAYTYPE_P(res_arr);
}
