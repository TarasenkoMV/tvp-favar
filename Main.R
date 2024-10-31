# Main
# rm(list = ls())
library(readxl)
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(xlsx)
source('C:\\Users\\Taras\\Desktop\\статья Dean\\tvp_dfm\\FUN_R.R')
sourceCpp('C:\\Users\\Taras\\Desktop\\статья Dean\\tvp_dfm\\tvp-dfm.cpp')

data = fread('C:\\Users\\Taras\\Desktop\\статья Dean\\data favar end IR\\data_saed.csv')[,-1] %>% as.matrix
transform = fread('C:\\Users\\Taras\\Desktop\\статья Dean\\data favar end IR\\transform.csv') %>% pull
Z = fread('C:\\Users\\Taras\\Desktop\\статья Dean\\data favar end IR\\Z.csv')[,-1] %>% as.matrix
vars_to_calc_irf = fread('C:\\Users\\Taras\\Desktop\\статья Dean\\data favar end IR\\vars_to_calc_irf.csv') %>% as.matrix %>% .[,1]
fast_slow = fread('C:\\Users\\Taras\\Desktop\\статья Dean\\data favar end IR\\fast_slow.csv') %>% as.matrix

if (length(transform) != ncol(data)) stop('transform or data not correct')

data_tf = matrix(0, nrow(data) - 1, ncol(data))

for (i in 1:ncol(data)) {
  if (transform[i] == 5) {
    data_tf[,i] = diff(log(data[,i]))
  } else if (transform[i] == 4) {
    data_tf[,i] = diff(data[,i])
  } else {
    data_tf[,i] = data[-1,i]
  }
}

Z = standartise(diff(Z))

data = data_tf
data = standartise(data)
# Очистить данные от ставки
data_new = data
x = cbind(Z, 1)
xTx_xT = solve(t(x) %*% x) %*% t(x)
for (i in 1:ncol(data)) {
  # if (fast_slow[i]) {
  # Инфляцию не чистим
  # if (i != 7) {
  data_new[,i] = data[,i,drop = F] - x %*% xTx_xT %*% data[,i,drop = F]
  # }
}

data = data_new

grid = expand.grid(k_B = c(4),
                   k_A_ = c(4),
                   k_sig_ = c(1),
                   k_Q_ = c(0.1, 0.2),
                   k_S_ = c(0.1, 0.2),
                   k_W_ = c(0.01, 0.1)
)

calc_irf = F

for (i in 8) {
  
  model = estimate_tvp_favar(
    data = data, 
    Z = Z, 
    # Зачем это все тут нагорожено: 
    # 1. Множитель для стартовой дисперсии коэффициентов
    # VAR модели для ФК. Стандартное из литературы
    k_B_ = grid$k_B[i],
    # 2. То же самое, но для внедиагональных элементов
    # разложения Холецкого. Стандартное значение из литературы
    k_A_ = grid$k_A_[i], 
    # 3. То же самое, для диагонали, т.е для дисперсии остатков
    k_sig_ = grid$k_sig_[i], 
    # 4. Множитель для приора на дисперсию коэффициентов.
    # Влияет на то, насколько "узкая" плотность для инноваций
    # коэффициентов VAR
    k_Q_ = grid$k_Q_[i],
    # 5. Множитель для приора на дисперсию внедиагональных элементов
    # разложения Холецкого
    k_S_ = grid$k_S_[i], 
    # 6. То же самое, но для диагональной матрицы, 
    # на которую умножается Холецкий
    k_W_ = grid$k_W_[i], 
    Reps = 40000, Burn = 33000, 
    # Как часто печатать отчет
    print_num = 50, 
    # Включена ли экзогенная переменная. Для которой шок.
    # Если не включена, то зачем вам эта модель
    Z_incl = T, 
    # Количество лагов в VAR
    L = 2, 
    # Сколько наблюдений использовать для стартовых значений для ФК
    # и для приоров. По сути регулирует только стартовые значения
    # Приоры в основном регулируются через параметры 3-6
    tau = 160, 
    # Количество ненаблюдаемых факторов
    K = 5, 
    # Порядок для Холецкого. Указывается, куда ставить Z.
    # Если Z - ставка, то общепринято ставить ее в конец
    ordering = c(0, 0, 0, 0, 0, 1), 
    # Надо ли считать irf (замедляет модель. По сути надо только для ДИ для irf)
    # Можно также на средних значениях после оценки посчитать irf
    calc_irf = calc_irf, 
    # Горизонт для irf. Обычно все затухает в течение года
    horizon = 15, 
    # Берется из файла. Отмечены переменные, для которых
    # надо посчитать irf. Можно отметить все, но тогда
    # ест много памяти
    vars_to_calc_irf = vars_to_calc_irf, 
    # Какие переменные быстрые-медленные.
    # См бернанке 2004 статья "Что-то там про Фавар"
    fast_slow = as.logical(fast_slow), 
    # Надо ли ротировать факторы перед началом оценки 
    # под первые K переменных
    # или оставить ротацию, которая возвращается ГК.
    # По умолчанию все думают, что тут стоит TRUE, т.к.
    # так принято в литературе, хотя это и не обязательно.
    # См Bai "Identification and Bayesian Estimation of Dynamic Factor Models"
    # Там в приложении доказательство. Оно выполняется для произвольной матрицы (K+1) x (K+1)
    rotate_prior_fact = T
  )
  
  # Сохраняет результаты оценки на диск в мои документы
  save_model(model = model, 
             i = i,
             calc_irf = calc_irf)
  
}