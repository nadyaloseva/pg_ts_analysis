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
 
