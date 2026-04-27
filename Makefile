# Makefile — pg_stl extension
# Использование:
#   make
#   make install
#   make installcheck   (требует psql и доступ к кластеру)

MODULE_big   = pg_stl
OBJS         = stl.o
EXTENSION    = pg_stl
DATA = pg_stl--1.0.sql
REGRESS      = stl_test

# Флаги компилятора: строгие предупреждения + math
PG_CFLAGS    = -Wall -Wextra -Wno-unused-parameter -O2
SHLIB_LINK   = -lm

PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
