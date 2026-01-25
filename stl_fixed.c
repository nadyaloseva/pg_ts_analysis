#include "postgres.h"
#include "fmgr.h"
#include "utils/array.h"
#include "catalog/pg_type.h"
#include "math.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(stl_decompose);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));
    double *data = (double *) ARR_DATA_PTR(arr);
    
    // ✅ ИСПРАВЛЕНО: Правильный размер (n элементов)
    Datum *result = palloc(n * sizeof(Datum));
    
    // 🎯 РЕАЛЬНЫЙ тренд: Экспоненциальное сглаживание
    double alpha = 0.3;  // Сглаживающий коэффициент
    double trend = data[0];
    
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            trend = alpha * data[i] + (1 - alpha) * trend;
        }
        result[i] = Float8GetDatum(trend);
    }
    
    // ✅ Возвращаем ТОЛЬКО ТРЕНД (правильный размер!)
    ArrayType *res_arr = construct_array(result, n, 2278, 
                                        sizeof(double), true, 'i');
    PG_RETURN_ARRAYTYPE_P(res_arr);
}

void _PG_init(void) {}

#include "postgres.h"
#include "fmgr.h"
#include "utils/array.h"
#include "catalog/pg_type.h"
#include "math.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(stl_decompose);

Datum stl_decompose(PG_FUNCTION_ARGS)
{
    ArrayType *arr = PG_GETARG_ARRAYTYPE_P(0);
    int32 period = PG_GETARG_INT32(1);
    int n = ArrayGetNItems(ARR_NDIM(arr), ARR_DIMS(arr));
    double *data = (double *) ARR_DATA_PTR(arr);
    
    // ✅ ИСПРАВЛЕНО: Правильный размер (n элементов)
    Datum *result = palloc(n * sizeof(Datum));
    
    // 🎯 РЕАЛЬНЫЙ тренд: Экспоненциальное сглаживание
    double alpha = 0.3;  // Сглаживающий коэффициент
    double trend = data[0];
    
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            trend = alpha * data[i] + (1 - alpha) * trend;
        }
        result[i] = Float8GetDatum(trend);
    }
    
    // ✅ Возвращаем ТОЛЬКО ТРЕНД (правильный размер!)
    ArrayType *res_arr = construct_array(result, n, 2278, 
                                        sizeof(double), true, 'i');
    PG_RETURN_ARRAYTYPE_P(res_arr);
}

void _PG_init(void) {}

