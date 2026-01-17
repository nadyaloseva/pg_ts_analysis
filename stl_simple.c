#include "postgres.h"
#include "fmgr.h"

PG_FUNCTION_INFO_V1(stl_decompose);
Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType  *input = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    int n = ARR_DIMS(input)[0];
    double *data = (double *)ARR_DATA_PTR(input);
    
    Datum *result = palloc(sizeof(Datum) * n * 3);
    
    // STL: trend (cumulative avg), seasonal, remainder
    double trend_sum = 0;
    for(int i = 0; i < n; i++) {
        trend_sum += data[i];
        double trend_val = trend_sum / (i + 1);
        double seasonal_val = sin(2 * M_PI * i / period) * 0.1; // dummy seasonal
        double remainder_val = data[i] - trend_val - seasonal_val;
        
        result[i*3]   = Float8GetDatum(trend_val);
        result[i*3+1] = Float8GetDatum(seasonal_val);
        result[i*3+2] = Float8GetDatum(remainder_val);
    }
    
    PG_RETURN_ARRAYTYPE_P(construct_array(result, n*3, FLOAT8OID,
                                          sizeof(double), true, 'i'));
}
