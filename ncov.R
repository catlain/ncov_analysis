library(jsonlite)
library(tidyverse)
library(lubridate)
library(deSolve)
library(purrr)

jsondata <- fromJSON('https://interface.sina.cn/news/wap/fymap2020_data.d.json')
# 当前时间
(times <- jsondata$data$times)
# "截至2月5日22时01分"

# 确诊数量
(confirm <- jsondata$data$gntotal)
# "24363"

# 死亡数量
(dead <- jsondata$data$deathtotal)
#> [1] "490"

# 疑似数量
(suspect <- jsondata$data$sustotal)
#> [1] "23260"

# 治愈数量
(cure <- jsondata$data$curetotal)
#> [1] "892"

str(jsondata)

# 历史数据
historylist <- jsondata$data$historylist %>%
  mutate(
    date = ymd(paste0("2020.", date))
  ) %>%
  as_tibble() %>%
  type_convert()
historylist


# 损失函数
MSE <- function(Infected, fit) {
  sum((Infected - fit)^2) / length(Infected)
}
MAE <- function(Infected, fit) {
  sum(abs(Infected - fit)) / length(Infected)
}

R2 <- function(Infected, fit) {
  1 - MSE(Infected, fit)/var(Infected)
}

# -------------基础数据-------------------

# 感染者人数

Type <- "wuhan"

switch (Type,
        "wuhan" = {
          Dat <- historylist %>% pull(wuhan_conNum) %>% rev()
          N <- 15000000
          StartDay <- 6
        },
        "other" = {
          Dat <- (historylist$cn_conNum - historylist$wuhan_conNum) %>% rev()
          N <- 1400000000 - 15000000
          StartDay <- 11
        },
        "total" = {
          Dat <- historylist$cn_conNum %>% rev()
          N <- 1400000000
          StartDay <- 6
        }
)
# Dat <- (historylist$cn_conNum - historylist$wuhan_conNum) %>% rev()

Infected <- Dat[StartDay:length(Dat)]
# 天数
Days <- 1:(length(Infected))
# # 感染最大范围?
# N <- 140000000

# 以多少天窗口计算R0
WatchingWindow <- 12

init <- c(S = N - Infected[1], I = Infected[1], R = 0)

# ---------------计算R0--------------------

# 首先写一个计算 R0 的函数：
R0 <- function(days){
  
  # WatchingWindow = 10
  # WatchingWindow + 1 = 11
  # 11 - 10 = 1
  # 11 - 1 = 10
  # day 1~10
  
  WindowInfected <- Infected[(days - WatchingWindow):(days - 1)]
  Days <- 1:(length(WindowInfected))
  init <- c(S = N - WindowInfected[1], I = WindowInfected[1], R = 0)
  
  LOSS <- function(parameters) {
    names(parameters) <- c("beta", "gamma")
    out <- ode(y = init, times = Days, func = SIR, parms = parameters)
    fit <- out[ , 3]
    -R2(WindowInfected,fit)
    # MSE(WindowInfected, fit)
  }
  
  Opt = optim(c(0.5, 0.5), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, Inf)) # gamma <= 0.5 ~ 感染到恢复 >= 2日
  Opt_par = setNames(Opt$par, c("beta", "gamma"))
  as_tibble(list(days = days, 
                 R0 = Opt_par["beta"] / Opt_par["gamma"],
                 beta = Opt_par["beta"], 
                 gamma = Opt_par["gamma"], 
                 value = Opt$value,
                 date = ymd("2020-01-11") + days(StartDay + days - 1) - 1)) # StartDay - 1 | days - 1
}

RODF <- map_df((WatchingWindow + 1):length(Infected), R0)
RODF

qplot(date, R0, data = filter(RODF, gamma > 0, gamma < 1), 
      geom = "point", 
      main = paste("Trend of", Type, "'s R0")) +
  geom_text(aes(date, R0, label = round(R0, 3)), nudge_x = 0.05, nudge_y = 0.02) +
  scale_x_date(date_breaks = "2 day")  



# --------------作图分析------------------

# 为什么要/N?
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

LOSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Days, func = SIR, parms = parameters)
  fit <- out[ , 3]
  # MSE(Infected,fit)^2/var(Infected)
  -R2(Infected,fit)
}

# 最优化
Opt <- optim(c(0.5, 0.5), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1))
Opt$message
#> "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
#>      beta     gamma
#> 0.6701221 0.3298780

R0 <- setNames(Opt_par["beta"] / Opt_par["gamma"], "R0")
R0


# 时间：1 ～ 70 天
t <- 1:100
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par)) %>%
  as_tibble()
fit

# 感染人数
plot(fit$I)

max(fit$I)

# 检查拟合情况
plot(fit$I[seq_len(length(Infected))], Infected)

fit %>%
  gather(S, I, R, key = "kind", value = "value") %>%
  ggplot(aes(time, value, color = kind)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(
    # breaks = c(0, 5e+08, 1e+9, 1.5e+9),
    # labels = c("0", "5 亿", "10 亿", "15 亿"),
    # expand = c(0, 0.2e+9)
  ) +
  labs(title = "新型冠状病毒肺炎疫情的 SIR 模型",
       x = "天数", y = "人数", color = "")
