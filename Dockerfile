FROM postgres:16

# Зависимости для сборки C-расширения
RUN apt-get update && apt-get install -y \
    build-essential \
    postgresql-server-dev-16 \
    && rm -rf /var/lib/apt/lists/*

# Копируем исходники расширения
WORKDIR /build/pg_stl
COPY stl.c         .
COPY pg_stl--1.0.sql .
COPY pg_stl.control .
COPY Makefile      .

# Собираем и устанавливаем расширение
RUN make && make install

# Скрипт создания расширения выполняется первым (01_)
COPY docker/01_create_extension.sql /docker-entrypoint-initdb.d/01_create_extension.sql
