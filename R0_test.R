# Necessary libraries
library(R0)
library(nCov2019)
library(tidyverse)

# Initialize random number generator
set.seed(679236)

# ebola <- load_nCov2019("en")$data %>%
#   filter(province == 'Hunan') %>%
#   group_by(time, province) %>%
#   summarise(confirmed = sum(confirmed)) %>%
#   group_by(province) %>%
#   arrange(time) %>%
#   mutate(cum_confirmed = cumsum(confirmed))

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

# 获取数据
ebola <- load_nCov2019("en")$data %>%
  filter(province == 'Hunan') %>%
  group_by(time, province) %>%
  summarise(confirmed = sum(confirmed)) %>%
  arrange(time) %>%
  mutate(cum_confirmed = sum(confirmed))

TD_LOSS <- function(param, plot = FALSE) {
  # print(param)
  
  mGT <- generation.time("gamma", c(param[[1]], param[[2]]))
  TD <- est.R0.TD(ebola$confirmed, mGT, begin = 1L, end = nrow(ebola), nsim = 1e4, correct = TRUE)
  # plot(TD)
  
  res <- list()
  
  res$R2 <- R2(TD$pred, ebola$confirmed)
  res$MSE <- MSE(TD$pred, ebola$confirmed)
  
  if (plot) {
    res_plot <- tibble(time = ebola$time, R0 = TD$R, pred = TD$pred, confirmed = ebola$cum_confirmed) 
    
    # print(res_plot, n = 100)
    res$plot <- ggplot(res_plot, aes(x = time)) +
      geom_line(aes(y = pred), color = "red", stat = "identity") +
      geom_line(aes(y = confirmed), color = "blue", stat = "identity") +
      scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
      labs(title = "pred vs confirmed")
  }
  
  return(res)

}
params <- cross2(seq(2, 10, length.out = 20), seq(0.5, 5, length.out = 20))

res <- lapply(params, TD_LOSS)

temp <- tibble(R2 = unlist(res)[seq(1, 199, by = 2)], 
               MSE = unlist(res)[seq(2, 200, by = 2)], 
               ci_mean = unlist(params)[seq(1, 199, by = 2)], 
               ci_std = unlist(params)[seq(2, 200, by = 2)])

arrange(temp, desc(R2)) %>% print(n = 100)


params_1 <- cross2(seq(1, 2, length.out = 10), seq(0.5, 2, length.out = 10))
res_1 <- lapply(params_1, TD_LOSS)

temp_1 <- tibble(R2 = unlist(res_1)[seq(1, 199, by = 2)], 
               MSE = unlist(res_1)[seq(2, 200, by = 2)], 
               ci_mean = unlist(params)[seq(1, 199, by = 2)], 
               ci_std = unlist(params)[seq(2, 200, by = 2)])

arrange(temp_1, desc(R2)) %>% print(n = 100)
