#include "postgres.h"
#include "fmgr.h"
#include "lib/statslin.h"

PG_FUNCTION_INFO_V1(stl_decompose);
Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType  *input_array = PG_GETARG_ARRAYTYPE_P(0);
    int32       period = PG_GETARG_INT32(1);
    int         n = ARR_DIMS(input_array)[0];
    
    double     *data = (double *) ARR_DATA_PTR(input_array);
    double      trend[1000], seasonal[1000], remainder[1000];
    
    // STL: trend = moving average, seasonal = detrend, remainder = residual
    double sum = 0;
    for(int i = 0; i < n; i++) {
        sum += data[i];
        trend[i] = sum / (i+1);
        seasonal[i] = data[i] - trend[i];
        remainder[i] = data[i] - trend[i] - seasonal[i];
    }
    
    Datum *result = palloc(sizeof(Datum) * n * 3);
    for(int i = 0; i < n; i++) {
        result[i*3]   = Float8GetDatum(trend[i]);
        result[i*3+1] = Float8GetDatum(seasonal[i]);
        result[i*3+2] = Float8GetDatum(remainder[i]);
    }
    
    PG_RETURN_ARRAYTYPE_P(construct_array(result, n*3, FLOAT8OID, 
                                          sizeof(double), true, 'i'));
}

    PG_RETURN_ARRAYTYPE_P(construct_array(result, n*3, FLOAT8OID, 
                                          sizeof(double), true, 'i'));
}
