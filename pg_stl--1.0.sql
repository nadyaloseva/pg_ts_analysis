-- stl--1.0.sql
-- pg_stl: STL time-series decomposition
-- Requires: CREATE EXTENSION pg_stl;
 

CREATE FUNCTION acf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'acf_array' LANGUAGE C STRICT IMMUTABLE;

CREATE FUNCTION pacf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'pacf_array' LANGUAGE C STRICT IMMUTABLE;


-- -----------------------------------------------------------------------
-- Составной тип результата для stl_decompose
-- -----------------------------------------------------------------------
CREATE TYPE stl_result AS (
    trend    double precision[],
    seasonal double precision[],
    residual double precision[]
);
 
-- -----------------------------------------------------------------------
-- stl_decompose — полная декомпозиция, возвращает composite type
--
-- Параметры:
--   y            double precision[]  — входной ряд
--   period       integer             — период сезонности
--   seasonal     integer  DEFAULT 7  — длина сглаживателя (нечётное, >= 3)
--   robust       boolean  DEFAULT true
--   trend        integer  DEFAULT 0  — 0 = авто (smallest odd > 1.5p/(1-1.5/s))
--   low_pass     integer  DEFAULT 0  — 0 = авто (smallest odd > period)
--   inner_iter   integer  DEFAULT 2
--   outer_iter   integer  DEFAULT 0  — 0 = авто (15 при robust=true)
-- -----------------------------------------------------------------------
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
 
COMMENT ON FUNCTION stl_decompose IS
'STL decomposition (Cleveland 1990). Returns composite (trend[], seasonal[], residual[]).
Example: SELECT (stl_decompose(y, 12)).trend FROM my_table;';
 
-- -----------------------------------------------------------------------
-- Удобные скалярные функции — возвращают только нужную компоненту
-- -----------------------------------------------------------------------
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
 
COMMENT ON FUNCTION stl_trend IS
'STL decomposition: returns trend component array.';
 
-- -----------------------------------------------------------------------
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
 
COMMENT ON FUNCTION stl_seasonal IS
'STL decomposition: returns seasonal component array.';
 
-- -----------------------------------------------------------------------
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
 
COMMENT ON FUNCTION stl_residual IS
'STL decomposition: returns residual component array.';
 
-- -----------------------------------------------------------------------
-- Агрегатный helper: собирает значения колонки в упорядоченный массив
-- (удобно для использования совместно с stl_*)
-- -----------------------------------------------------------------------
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
 
COMMENT ON FUNCTION stl_collect_ordered IS
'Helper: collects a numeric column into an ordered array for STL input.
Example: SELECT stl_trend(stl_collect_ordered(''sales'', ''value'', ''date''), 12);';