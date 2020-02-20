# init
rm(list = ls())

# Necessary libraries
# library(chron)
library(R0)
library(nCov2019)
library(tidyverse)

# Initialize random number generator
set.seed(679236)

ebola <- load_nCov2019("en")$data %>%
  filter(province == 'Hunan') %>%
  group_by(time, province) %>%
  summarise(confirmed = sum(confirmed)) %>%
  arrange(time) %>%
  mutate(cum_confirmed = sum(confirmed))

# Estimating case reproduction number
# mGT <- generation.time("gamma", c(5.2,3.9)) # Generation time (serial interval) as reported by WHO Ebola Response Team (2014, NEJM)

mGT <- generation.time("gamma", c(3, 1.5))
TD <- est.R0.TD(ebola$cum_confirmed, mGT, begin = 1L, end = nrow(ebola), nsim = 1e4, correct = TRUE)
res_plot <- tibble(time = ebola$time, R0 = TD$R, pred = TD$pred, confirmed = ebola$cum_confirmed) 

# 实际跟预测的拟合
ggplot(res_plot, aes(x = time)) +
  geom_line(aes(y = pred), color = "red", stat = "identity") +
  geom_line(aes(y = confirmed), color = "blue", stat = "identity") +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
  labs(title = "pred vs confirmed")

# 预测的R0趋势
ggplot(res_plot, aes(x = time)) +
  geom_line(aes(y = R0), color = "red", stat = "identity") +
  geom_text(aes(y = R0, label = round(R0, 2))) +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
  labs(title = "R0 trend")