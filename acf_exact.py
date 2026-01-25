#!/usr/bin/env python3
import psycopg2
import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import acf
import time

conn = psycopg2.connect(host="localhost", database="demo", user="nadya")
start = time.time()

df = pd.read_sql("""
    SELECT EXTRACT(EPOCH FROM GREATEST('0 minutes'::interval, actual_arrival - scheduled_arrival))/3600.0 as delay_h
    FROM bookings.flights WHERE flight_no LIKE 'PG%' LIMIT 100000
""", conn)

acf_result = acf(df['delay_h'].values, nlags=7, fft=False)
print(f"Python ACF: lag1={acf_result[1]:.4f}, lag7={acf_result[7]:.4f}")
print(f"Время: {time.time()-start:.2f} сек (vs SQL 0.8 сек)")
conn.close()
