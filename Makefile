EXTENSION = stl_ext
OBJS = stl_complete.o acf_array.o pacf_array.o
DATA_built = stl_ext--1.0.sql
DOCS = README.md
REGRESS = ts_analysis_test
EXTRA_CLEAN = pg_ts_analysis.sql

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
