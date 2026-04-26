Репозиторий создан для создания расширения PostgreSQL для анализа временных рядов. 

Состав расширение PostgreSQL для анализа временных рядов:
- Автокорреляционная функция (ACF)
- Частичная автокорреляционная функция (PACF)
- STL-декомпозиция

## Установка
```bash
make install
psql -d demo -c "CREATE EXTENSION stl_ext;"
```

Репозиторий содержит следующие файлы:

* Makefile – файл отвечает за сборку и установку расширения PostgreSQL
* stl_complete.c - файл с бизнес-логикой расширения, написанной на языке C (ACF, PACF, STL)
* stl_ext.control – файл с описанием метаданных расширения
* stl_ext.sql - файл, который связывает C-код с SQL

## Использование 

```sql
SELECT acf_array(ARRAY[1,2,3,4,5], 3);
SELECT pacf_array(ARRAY[1,2,3,4,5], 3);
SELECT stl_decompose(ARRAY[1,2,3,4,5], 7);
```


## Описание функций 

### acf_array
Автокорреляционная функция (ACF) используется для оценки зависимости значений временного ряда от его прошлых значений.

**Математическое определение**

Для временного ряда $x_1, x_2, \dots, x_n$ автокорреляция для лага $k$ определяется как:

$$
ACF(k) = \frac{\sum_{i=1}^{n-k}(x_i - \bar{x})(x_{i+k} - \bar{x})}{\sum_{i=1}^{n}(x_i - \bar{x})^2}
$$

где:
- $x_i$ — значение временного ряда;
- $\bar{x}$ — среднее значение ряда;
- $n$ — длина ряда;
- $k$ — лаг.
  
В коде функции
```
return cov_lag / var0;
```

**Описание алгоритма**

Вычисление среднего значения

$$
\bar{x} = \frac{1}{n} \sum_{i=1}^{n} x_i
$$

В коде функции
```
double mean = 0;
for (int i = 0; i < n; i++) mean += data[i];
mean /= n;
```
Центрирование временного ряда

Каждое значение преобразуется следующим образом:

$$
d_i = x_i - \bar{x}
$$

В коде функции
```
double d1 = data[i] - mean;
double d2 = data[i + lag] - mean;
```
Центрирование необходимо для корректного вычисления ковариации.

Вычисление ковариации для лага $k$

$$
cov(k) = \sum_{i=1}^{n-k} d_i \cdot d_{i+k}
$$

В коде функции
```
cov_lag += d1 * d2;
```
Вычисление дисперсии

$$
var = \sum_{i=1}^{n-k} d_i^2
$$

В коде функции
```
var0 += d1 * d1;
```
В реализации используется практическая форма оценки автокорреляции, при которой сумма в числителе и знаменателе вычисляется по одинаковому числу наблюдений (n − k). Это позволяет избежать смещения оценки при больших значениях лага.

Вычисление автокорреляции

$$
ACF(k) = \frac{cov(n)}{var}
$$

В коде функции
```
var0 += d1 * d1;
```
Если дисперсия близка к нулю, значение ACF принимается равным нулю для предотвращения деления на ноль.


### pacf_array
Частичная автокорреляционная функция показывает зависимость между текущим значением временного ряда и его значением на лаге $k$, исключая влияние промежуточных лагов.
В отличие от автокорреляционной функции (ACF), PACF отражает **прямую зависимость** между значениями ряда.

**Математическое определение**

PACF определяется как коэффициент $\phi_{k,k}$ в авторегрессионной модели порядка $k$:

$$
x_t = \phi_{k,1} x_{t-1} + \phi_{k,2} x_{t-2} + \dots + \phi_{k,k} x_{t-k} + \varepsilon_t
$$

где:
- $\phi_{k,k}$ — значение PACF на лаге $k$;
- $\varepsilon_t$ — случайная ошибка.

**Описание алгоритма**

На первом этапе вычисляется автоковариационная функция:

$$
r(k) = \frac{1}{n} \sum_{i=1}^{n-k}(x_i - \bar{x})(x_{i+k} - \bar{x})
$$

где:
- $\bar{x}$ — среднее значение ряда. Вычисляется среднее значение временного ряда, необходимое для центрирования данных.

В коде функции
```
for (int i = 0; i < n - k; i++) {
    sum += (x[i] - mean) * (x[i + k] - mean);
}
r[k] = (n > 0) ? sum / n : 0.0;
```

Для вычисления PACF используется рекурсивная схема Юла–Уокера.

Начальное значение дисперсии:

$$
\sigma_0 = r(0)
$$

Основная формула Юла–Уокера для каждого лага $k$:

Для каждого лага $k$:

$$
\phi_{k,k} = \frac{r(k) - \sum_{j=1}^{k-1} \phi_{k-1,j} \cdot r(k-j)}{\sigma_{k-1}}
$$


Обновление коэффициентов авторегрессионной модели.

$$
\phi_{k,j} = \phi_{k-1,j} - \phi_{k,k} \cdot \phi_{k-1,k-j}, \quad j = 1, \dots, k-1
$$


Обновление дисперсии ошибки

$$
\sigma_k = \sigma_{k-1} \cdot (1 - \phi_{k,k}^2)
$$

В коде функции
```
sigma[0] = r[0];
for (int k = 1; k <= max_lag; k++) {
    double sum = 0.0;
    for (int j = 1; j < k; j++) {sum += phi[k-1][j] * r[k - j];}
    phi[k][k] = (r[k] - sum) / sigma[k-1];
    for (int j = 1; j < k; j++) {phi[k][j] = phi[k-1][j] - phi[k][k] * phi[k-1][k-j];}
    sigma[k] = sigma[k-1] * (1.0 - phi[k][k] * phi[k][k]);
    out[k-1] = phi[k][k];}
}
```
Значение частичной автокорреляции для лага $k$ определяется как:

$$
PACF(k) = \phi_{k,k}
$$

В коде функции
```
out[k-1] = phi[k][k];}
```

Результаты для всех лагов сохраняются в массив и возвращаются пользователю.


### stl_decompose

STL (Seasonal-Trend decomposition using LOESS) — метод разложения временного ряда на три компоненты:

- тренд $T_t$
- сезонность $S_t$
- остаток $R_t$

$$
y_t = T_t + S_t + R_t
$$

Метод основан на локально-взвешенной регрессии (LOESS) и итеративной схеме уточнения компонент.


**1. Общая итеративная схема STL**

Алгоритм выполняется итеративно:

- удаляется тренд
- оценивается сезонность
- пересчитывается тренд
- при необходимости обновляются робастные веса

В коде функции
```
for (int o = 0; o < n_outer; o++) {
    const double *wrob = (o == 0) ? NULL : rw;

    for (int inner = 0; inner < cfg->inner_iter; inner++) {

        for (int i = 0; i < n; i++)
            detr[i] = y[i] - trend[i];

        seasonal_step(detr, n, cfg->period,
                      cfg->seasonal, cfg->low_pass,
                      cfg->seasonal_deg, cfg->low_pass_deg,
                      cfg->seasonal_jump, cfg->low_pass_jump,
                      wrob, season);

        for (int i = 0; i < n; i++)
            desea[i] = y[i] - season[i];

        trend_step(desea, n,
                   cfg->trend, cfg->trend_deg, cfg->trend_jump,
                   wrob, trend);
    }

    if (cfg->robust && o < n_outer - 1)
        compute_robust_weights(y, trend, season, n, rw);
}
```

**LOESS (локальная регрессия)**

Трикубическая функция весов:

$$
w(x)=\left(1-|x|^3\right)^3
$$

В коде функции

```
static double tricube(double x) {
    x = fabs(x);
    if (x >= 1.0) return 0.0;
    double t = 1.0 - x * x * x;
    return t * t * t;
}
```
**Локальная линейная регрессия**

$$
y = a + bx
$$

$$
b = \frac{\sum w x y - \sum w x \sum w y}{\sum w x^2 - (\sum w x)^2}
$$

$$
a = \frac{\sum w y - b \sum w x}{\sum w}
$$

В коде функции
```
double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;

for (int j = 0; j < nd; j++) {
    double w = tricube(dist[j] / dmax);
    if (wrob) w *= wrob[j];

    sw   += w;
    swx  += w * xd[j];
    swy  += w * yd[j];
    swxx += w * xd[j] * xd[j];
    swxy += w * xd[j] * yd[j];
}

double b = (sw * swxy - swx * swy) / (sw * swxx - swx * swx);
double a = (swy - b * swx) / sw;
```

**Jump-интерполяция LOESS**

В коде функции
```
for (int i = 0; i < nq; i += jump)
    sidx[ns++] = i;

if (sidx[ns - 1] != nq - 1)
    sidx[ns++] = nq - 1;
```
**Сезонная компонента**

$$
y^{(p)} = {y_p, y_{p+P}, y_{p+2P}, ...}
$$

В коде функции
```
for (int p = 0; p < period; p++) {

    int m = 0;
    for (int i = p; i < n; i += period) m++;

    double *xsub = malloc(m * sizeof(double));
    double *ysub = malloc(m * sizeof(double));

    int idx = 0;
    for (int i = p; i < n; i += period) {
        xsub[idx] = (double)idx;
        ysub[idx] = detr[i];
        idx++;
    }

    loess_fit(xsub, m, xsub, ysub, m,
              bw, deg_s, wsub, jump_s, yhat);

    idx = 0;
    for (int i = p; i < n; i += period)
        cv[i] = yhat[idx++];
}
```

**Низкочастотная фильтрация**

$$
S_t = C_v - LOESS(C_v)
$$

В коде функции
```
loess_fit(xs, n, xs, cv, n, n_l, deg_l, NULL, jump_l, lp);

for (int i = 0; i < n; i++)
    season[i] = cv[i] - lp[i];
```

**Тренд**

$$
y_t^{(deseasoned)} = y_t - S_t
$$

$$
T_t = LOESS(y_t^{(deseasoned)})
$$

В коде функции
```
for (int i = 0; i < n; i++)
    desea[i] = y[i] - season[i];

loess_fit(xs, n, xs, desea, n,
          n_t, deg_t, wrob, jump_t, trend);
```

**Остаток**

$$
R_t = y_t - T_t - S_t
$$

В коде функции
```
for (int i = 0; i < n; i++)
    residual[i] = y[i] - trend[i] - season[i];
```

**Робастные веса**

$$
MAD = median(|r_t - median(r)|)
$$

$$
w_t = \left(1 - \left(\frac{r_t}{6 \cdot MAD}\right)^2\right)^2
$$

В коде функции
```c
for (int i = 0; i < n; i++)
    r[i] = y[i] - trend[i] - season[i];

double mad = mad_vec(r, n);
double c = 6.0 * mad;

for (int i = 0; i < n; i++)
    rw[i] = bisquare(r[i] / c);
```

**Автовыбор параметров**

Тренд:

$$
n_t = \text{next odd integer} > \frac{1.5P}{1 - 1.5/n_s}
$$


В коде функции
```
double t = 1.5 * c->period / (1.0 - 1.5 / (double)c->seasonal);
c->trend = next_odd_gt(t);
```

**Низкочастотное окно:**

$$
n_l > P
$$

В коде функции
```
c->low_pass = next_odd_gt((double)c->period);
```
