-- Time Series Analysis Functions
CREATE FUNCTION acf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'acf_array' LANGUAGE C STRICT IMMUTABLE;

CREATE FUNCTION pacf_array(data double precision[], lags integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'pacf_array' LANGUAGE C STRICT IMMUTABLE;

CREATE FUNCTION stl_decompose(data double precision[], period integer)
RETURNS double precision[] AS 'MODULE_PATHNAME', 'stl_decompose' LANGUAGE C STRICT IMMUTABLE;
