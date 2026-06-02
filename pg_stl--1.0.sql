CREATE FUNCTION acf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'acf_array' LANGUAGE C STRICT IMMUTABLE;

CREATE FUNCTION pacf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'pacf_array' LANGUAGE C STRICT IMMUTABLE;

-- Составной тип результата для stl_decompose
CREATE TYPE stl_result AS (
    trend    double precision[],
    seasonal double precision[],
    residual double precision[]
);
 
-- Функция алгоритма STL
CREATE OR REPLACE FUNCTION stl_decompose(
    y            double precision[],
    period       integer,
    seasonal     integer  DEFAULT 7,
    robust       boolean  DEFAULT true,
    trend        integer  DEFAULT 0,
    low_pass     integer  DEFAULT 0,
    inner_iter   integer  DEFAULT 2,
    outer_iter   integer  DEFAULT 0
)
RETURNS stl_result
AS '$libdir/pg_stl', 'stl_decompose'
LANGUAGE C STRICT IMMUTABLE PARALLEL SAFE;

-- Вспомогательные функции для удобства
CREATE OR REPLACE FUNCTION stl_trend(
    y            double precision[],
    period       integer,
    seasonal     integer  DEFAULT 7,
    robust       boolean  DEFAULT true,
    trend        integer  DEFAULT 0,
    low_pass     integer  DEFAULT 0,
    inner_iter   integer  DEFAULT 2,
    outer_iter   integer  DEFAULT 0
)
RETURNS double precision[]
AS '$libdir/pg_stl', 'stl_trend'
LANGUAGE C STRICT IMMUTABLE PARALLEL SAFE;
 

CREATE OR REPLACE FUNCTION stl_seasonal(
    y            double precision[],
    period       integer,
    seasonal     integer  DEFAULT 7,
    robust       boolean  DEFAULT true,
    trend        integer  DEFAULT 0,
    low_pass     integer  DEFAULT 0,
    inner_iter   integer  DEFAULT 2,
    outer_iter   integer  DEFAULT 0
)
RETURNS double precision[]
AS '$libdir/pg_stl', 'stl_seasonal'
LANGUAGE C STRICT IMMUTABLE PARALLEL SAFE;
 
CREATE OR REPLACE FUNCTION stl_residual(
    y            double precision[],
    period       integer,
    seasonal     integer  DEFAULT 7,
    robust       boolean  DEFAULT true,
    trend        integer  DEFAULT 0,
    low_pass     integer  DEFAULT 0,
    inner_iter   integer  DEFAULT 2,
    outer_iter   integer  DEFAULT 0
)
RETURNS double precision[]
AS '$libdir/pg_stl', 'stl_residual'
LANGUAGE C STRICT IMMUTABLE PARALLEL SAFE;
 

-- Агрегатная функция, которая собирает значения колонки в упорядоченный массив

CREATE OR REPLACE FUNCTION stl_collect_ordered(
    tbl  regclass,
    val  text,
    ord  text
)
RETURNS double precision[]
LANGUAGE plpgsql STABLE AS $$
DECLARE
    result double precision[];
BEGIN
    EXECUTE format(
        'SELECT array_agg(%I ORDER BY %I) FROM %s',
        val, ord, tbl
    ) INTO result;
    RETURN result;
END;
$$;

-- ============================================================
-- Holt-Winters forecasting (PL/pgSQL)
-- Author of algorithm: Pugolovok Alexander (tg: @alexander_pu)
-- Adapted for pg_stl extension: data is passed as array parameter
-- instead of reading from a fixed table.
-- Requires PostgreSQL 16+ (uses random() * x for coefficients).
-- ============================================================
 
CREATE OR REPLACE FUNCTION holt_winters_mse(
    seasonal_type         text,
    start_data_array      real[],
    row_data_length       int,
    period_length         int,
    mean_data_first_period real,
    levell                real[],
    trend                 real[],
    seasonal              real[],
    forecast              real[],
    alpha                 real,
    beta                  real,
    gamma                 real
)
/*
  Compute total squared error of the Holt-Winters forecast vs. the given series.
 
  Parameters:
    seasonal_type          : 'mult' or 'add'
    start_data_array       : input time series
    row_data_length        : number of elements in the series
    period_length          : length of one seasonal cycle
    mean_data_first_period : mean of the first period (used for initialisation)
    levell                 : array of level values (pre-initialised)
    trend                  : array of trend values (pre-initialised)
    seasonal               : array of seasonal values (pre-initialised)
    forecast               : array of forecast values (pre-initialised)
    alpha                  : smoothing coefficient for level  (0 <= alpha < 1)
    beta                   : smoothing coefficient for trend  (beta < alpha)
    gamma                  : smoothing coefficient for seasonality (gamma < 1 - alpha)
 
  Returns:
    Total squared error (real).
*/
RETURNS real
LANGUAGE plpgsql
AS $$
DECLARE
    mse_error real := 0.0;
BEGIN
    IF seasonal_type = 'mult' THEN
        FOR i IN period_length + 1 .. row_data_length LOOP
            levell[i]  := alpha * (start_data_array[i] / seasonal[i - period_length])
                          + (1 - alpha) * (levell[i-1] + trend[i-1]);
            trend[i]   := beta  * (levell[i] - levell[i-1]) + (1 - beta)  * trend[i-1];
            seasonal[i]:= gamma * (start_data_array[i] / levell[i])
                          + (1 - gamma) * seasonal[i - period_length];
            forecast[i]:= (levell[i-1] + trend[i-1]) * seasonal[i - period_length];
        END LOOP;
 
        FOR i IN period_length + 1 .. row_data_length LOOP
            mse_error := mse_error + (forecast[i] - start_data_array[i])^2;
        END LOOP;
    END IF;
 
    IF seasonal_type = 'add' THEN
        FOR i IN period_length + 1 .. row_data_length LOOP
            levell[i]  := alpha * (start_data_array[i] - seasonal[i - period_length])
                          + (1 - alpha) * (levell[i-1] + trend[i-1]);
            trend[i]   := beta  * (levell[i] - levell[i-1]) + (1 - beta)  * trend[i-1];
            seasonal[i]:= gamma * (start_data_array[i] - levell[i])
                          + (1 - gamma) * seasonal[i - period_length];
            forecast[i]:= levell[i-1] + trend[i-1] + seasonal[i - period_length];
        END LOOP;
 
        FOR i IN period_length + 1 .. row_data_length LOOP
            mse_error := mse_error + (forecast[i] - start_data_array[i])^2;
        END LOOP;
    END IF;
 
    RETURN mse_error;
END;
$$;
 
 
CREATE OR REPLACE FUNCTION holt_winters_predict(
    seasonal_type     text,
    period_length     int,
    start_data_array  real[]
)
/*
  Forecast one full seasonal cycle ahead using the Holt-Winters method.
 
  Parameters:
    seasonal_type    : 'mult' (multiplicative) or 'add' (additive)
    period_length    : length of one seasonal cycle (e.g. 4 for quarterly, 12 for monthly)
    start_data_array : the observed time series as a real[] array
 
  Returns:
    Table of forecasted values (one row per period step ahead).
 
  Example:
    SELECT * FROM holt_winters_predict(
        'mult',
        4,
        ARRAY[362,385,432,341,382,409,498,387,473,513,582,474,544,582,681,557,628,707,773,592,627,725,854,661]::real[]
    );
 
  Notes:
    - Exponential smoothing coefficients (alpha, beta, gamma) are selected automatically
      by minimising the total squared error via 500 random initialisations followed by
      gradient-style refinement in steps of 0.001.
    - Constraints: alpha in [0, 1), beta < alpha, gamma < 1 - alpha.
    - Compatible with PostgreSQL 16+.
*/
RETURNS TABLE(val real)
LANGUAGE plpgsql
AS $$
DECLARE
    row_data_length        integer;
    mean_data_first_period real    := 0.0;
    alpha                  numeric := 0.0;
    beta                   numeric := 0.0;
    gamma                  numeric := 0.0;
    levell                 real[];
    trend                  real[];
    seasonal               real[];
    forecast               real[];
    forecast_ahead         real[];
    err_alf_bet_gam_arr    real[][];
    array_error_args       real[][];
    tmp_mse_error          real    := 0.0;
    mse_error              real    := 0.0;
    tmp_for_trend_1        real    := 0.0;
    tmp_arg_dif            real    := 0.0;
    cnt_in_loop            integer := 0;
BEGIN
    row_data_length := array_length(start_data_array, 1);
 
    -- Mean of first period
    mean_data_first_period := (
        SELECT ROUND(AVG(x)::numeric, 3)
        FROM unnest(start_data_array[1:period_length]) AS x
    );
 
    -- Initial seasonality
    IF seasonal_type = 'mult' THEN
        FOR i IN 1 .. period_length LOOP
            seasonal := array_append(seasonal, start_data_array[i] / mean_data_first_period);
        END LOOP;
    END IF;
    IF seasonal_type = 'add' THEN
        FOR i IN 1 .. period_length LOOP
            seasonal := array_append(seasonal, start_data_array[i] - mean_data_first_period);
        END LOOP;
    END IF;
    FOR i IN period_length + 1 .. row_data_length LOOP
        seasonal := array_append(seasonal, NULL);
    END LOOP;
 
    -- Initial level
    FOR i IN 1 .. period_length - 1 LOOP
        levell := array_append(levell, NULL);
    END LOOP;
    levell := array_append(levell, mean_data_first_period);
    FOR i IN period_length + 1 .. row_data_length LOOP
        levell := array_append(levell, NULL);
    END LOOP;
 
    -- Initial trend
    FOR i IN 1 .. period_length - 1 LOOP
        trend := array_append(trend, NULL);
    END LOOP;
    FOR i IN 1 .. period_length LOOP
        tmp_for_trend_1 := tmp_for_trend_1
            + (start_data_array[period_length + i] - start_data_array[i]) / period_length^2;
    END LOOP;
    trend := array_append(trend, tmp_for_trend_1);
    FOR i IN period_length + 1 .. row_data_length LOOP
        trend := array_append(trend, NULL);
    END LOOP;
 
    -- Initial forecast array
    forecast := array_fill(0.0, ARRAY[row_data_length]);
    FOR i IN 1 .. period_length LOOP
        forecast[i] := NULL;
    END LOOP;
 
    -- Step 1: random search for best alpha/beta/gamma (500 iterations)
    FOR _ IN 1 .. 500 LOOP
        alpha     := random() * 0.999;
        beta      := random() * alpha;
        gamma     := random() * (1.0 - alpha);
        mse_error := holt_winters_mse(
            seasonal_type, start_data_array, row_data_length, period_length,
            mean_data_first_period, levell, trend, seasonal, forecast,
            alpha::real, beta::real, gamma::real
        );
        err_alf_bet_gam_arr := array_cat(
            err_alf_bet_gam_arr,
            ARRAY[ARRAY[mse_error, alpha::real, beta::real, gamma::real]]
        );
    END LOOP;
 
    -- Pick best from random search
    SELECT e[1], e[2], e[3], e[4]
    INTO mse_error, alpha, beta, gamma
    FROM (
        SELECT ARRAY[t[1], t[2], t[3], t[4]] AS e
        FROM unnest(err_alf_bet_gam_arr) WITH ORDINALITY AS u(t, ord)
        ORDER BY t[1]
        LIMIT 1
    ) sub;
 
    -- Recompute with best coefficients
    mse_error := holt_winters_mse(
        seasonal_type, start_data_array, row_data_length, period_length,
        mean_data_first_period, levell, trend, seasonal, forecast,
        alpha::real, beta::real, gamma::real
    );
 
    -- Step 2: gradient refinement in steps of 0.001
    LOOP
        cnt_in_loop := cnt_in_loop + 1;
        EXIT WHEN cnt_in_loop >= 500;
 
        array_error_args := NULL;
 
        -- alpha + 0.001
        IF alpha < 0.999 THEN
            tmp_arg_dif   := alpha + 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, tmp_arg_dif::real, beta::real, gamma::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := alpha;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, tmp_arg_dif::real, beta::real, gamma::real]]);
 
        -- alpha - 0.001
        IF alpha >= 0.001 THEN
            tmp_arg_dif   := alpha - 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, tmp_arg_dif::real, beta::real, gamma::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := alpha;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, tmp_arg_dif::real, beta::real, gamma::real]]);
 
        -- beta + 0.001
        IF beta < alpha - 0.001 THEN
            tmp_arg_dif   := beta + 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, alpha::real, tmp_arg_dif::real, gamma::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := beta;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, alpha::real, tmp_arg_dif::real, gamma::real]]);
 
        -- beta - 0.001
        IF beta >= 0.001 THEN
            tmp_arg_dif   := beta - 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, alpha::real, tmp_arg_dif::real, gamma::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := beta;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, alpha::real, tmp_arg_dif::real, gamma::real]]);
 
        -- gamma + 0.001
        IF gamma < 0.999 - alpha THEN
            tmp_arg_dif   := gamma + 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, alpha::real, beta::real, tmp_arg_dif::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := gamma;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, alpha::real, beta::real, tmp_arg_dif::real]]);
 
        -- gamma - 0.001
        IF gamma >= 0.001 THEN
            tmp_arg_dif   := gamma - 0.001;
            tmp_mse_error := holt_winters_mse(seasonal_type, start_data_array, row_data_length, period_length,
                mean_data_first_period, levell, trend, seasonal, forecast, alpha::real, beta::real, tmp_arg_dif::real);
        ELSE
            tmp_mse_error := mse_error; tmp_arg_dif := gamma;
        END IF;
        array_error_args := array_cat(array_error_args, ARRAY[ARRAY[tmp_mse_error, alpha::real, beta::real, tmp_arg_dif::real]]);
 
        -- Pick best candidate
        SELECT e[1], e[2], e[3], e[4]
        INTO tmp_mse_error, tmp_arg_dif, beta, gamma
        FROM (
            SELECT ARRAY[t[1], t[2], t[3], t[4]] AS e
            FROM unnest(array_error_args) WITH ORDINALITY AS u(t, ord)
            ORDER BY t[1]
            LIMIT 1
        ) sub;
 
        EXIT WHEN tmp_mse_error >= mse_error;
 
        mse_error := tmp_mse_error;
        alpha     := tmp_arg_dif;
    END LOOP;
 
    -- Final pass: compute level/trend/seasonal with optimised coefficients
    IF seasonal_type = 'mult' THEN
        FOR i IN period_length + 1 .. row_data_length LOOP
            levell[i]  := alpha::real * (start_data_array[i] / seasonal[i - period_length])
                          + (1 - alpha::real) * (levell[i-1] + trend[i-1]);
            trend[i]   := beta::real  * (levell[i] - levell[i-1]) + (1 - beta::real)  * trend[i-1];
            seasonal[i]:= gamma::real * (start_data_array[i] / levell[i])
                          + (1 - gamma::real) * seasonal[i - period_length];
        END LOOP;
 
        FOR i IN row_data_length + 1 .. row_data_length + period_length LOOP
            forecast_ahead := array_append(forecast_ahead,
                (levell[row_data_length] + trend[row_data_length] * (i - row_data_length))
                * seasonal[i - period_length]);
        END LOOP;
    END IF;
 
    IF seasonal_type = 'add' THEN
        FOR i IN period_length + 1 .. row_data_length LOOP
            levell[i]  := alpha::real * (start_data_array[i] - seasonal[i - period_length])
                          + (1 - alpha::real) * (levell[i-1] + trend[i-1]);
            trend[i]   := beta::real  * (levell[i] - levell[i-1]) + (1 - beta::real)  * trend[i-1];
            seasonal[i]:= gamma::real * (start_data_array[i] - levell[i])
                          + (1 - gamma::real) * seasonal[i - period_length];
        END LOOP;
 
        FOR i IN row_data_length + 1 .. row_data_length + period_length LOOP
            forecast_ahead := array_append(forecast_ahead,
                levell[row_data_length] + trend[row_data_length] * (i - row_data_length)
                + seasonal[i - period_length]);
        END LOOP;
    END IF;
 
    -- Return forecast
    FOR i IN 1 .. period_length LOOP
        val := forecast_ahead[i];
        RETURN NEXT;
    END LOOP;
END;
$$;
 
