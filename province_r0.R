# Necessary libraries
library(R0)
library(nCov2019)
library(tidyverse)

# Initialize random number generator
set.seed(679236)

dat <- load_nCov2019("en")$data %>%
  group_by(time, province) %>%
  summarise(confirmed = sum(confirmed)) %>%
  group_by(province) %>%
  arrange(time) %>%
  mutate(cum_confirmed = sum(confirmed))

dat_provinces <- filter(dat, time == '2020-02-17') %>%
  arrange(desc(cum_confirmed)) 
  

dates <- unique(dat$time) %>% sort()

# mGT <- generation.time("gamma", c(8.5, 3.8)) #8.5+-3.8是SARS的，现在应该有新病毒的参数了
mGT <- generation.time("gamma", c(3, 1.5))
get_est <- function(dat, dates){
  names(dat) <- dates
  est <- est.R0.TD(dat, mGT, nsim = 1e4, correct = TRUE, begin = 4L, end = length(dat))
  return(est)
}

get_sen <- function(dat, dates){
  names(dat) <- dates
  sen <- sensitivity.analysis(sa.type = "time", incid = dat, GT = mGT, begin = 4L:14L, end = (length(dat)-10):length(dat), est.method = "EG")
  return(sen)
}

resRaw <- group_by(dat, province) %>%
  summarise(est = list(get_est(confirmed, time)), 
            sen = list(get_sen(confirmed, time)), 
            time = list(time))


filter(resRaw, province == "Shanghai") %>% 
  pull(sen) %>% 
  `[[`(1) %>%  
  plot(what=c("criterion", "heatmap"))

res <- lapply(1:nrow(resRaw), function(x){
  # message(x)
  temp <- resRaw[x,]
  tibble(province = temp$province[[1]], 
         R0 = temp$est[[1]]$R,
         pred = temp$est[[1]]$pred,
         lower = temp$est[[1]]$conf.int$lower,
         upper = temp$est[[1]]$conf.int$upper,
         time = temp$time[[1]]
         )
}) %>%
  bind_rows()
res <- inner_join(res, dat) %>%
  arrange(desc(cum_confirmed))

ggplot(filter(res, province %in% dat_provinces$province[21:31], time >= "2020-01-30", time < "2020-02-19"), 
       aes(x = time, y = R0, color = province, group = province)) +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = "dodge") +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") 



ggplot(filter(res, province %in% c("Shandong", "Hubei", "Guangdong", "Zhejiang"), time >= "2020-01-30", time < "2020-02-19"), 
       aes(x = time, y = R0, color = province, group = province)) +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = "dodge") +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") 


ggplot(filter(res, province %in% c("Hainan", "Hebei", "Zhejiang", "Guangdong"), time >= "2020-01-30", time < "2020-02-19"), 
       aes(x = time, y = R0, color = province, group = province)) +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = "dodge") +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") 

ggplot(filter(res, province %in% c("Inner Mongolia", "Xinjiang", "Tibet", "Zhejiang", "Guangdong"), time >= "2020-01-30", time < "2020-02-19"), 
       aes(x = time, y = R0, color = province, group = province)) +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, position = "dodge") +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") 

  
  
  
  facet_grid(rows = vars(province))


  
  data(Germany.1918)
  mGT<-generation.time("gamma", c(2.6,1))
  ## sensitivity analysis for begin between day 1 and 15, and end between day 16 and 30
  sen = sensitivity.analysis(sa.type="time", incid=Germany.1918, GT=mGT, begin=1:15, end=16:30, 
                             est.method="EG")
  
  plot(sen, what=c("criterion","heatmap"))