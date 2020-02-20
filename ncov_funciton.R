# remotes::install_github("microly/alimap")
library(jsonlite)
library(tidyverse)
library(lubridate)
library(deSolve)
# library(alimap)
# library(sf)
# library(gganimate)
# # library(Cairo)
# library(magick)

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

# 主函数
getTrendR0 <- function(Infected, WatchingWindow) {
  
  Infected <- as.numeric(Infected) # named numeric not work
  Ndays <- length(Infected)
  
  if (Ndays < WatchingWindow) {
    return(NULL)
  }
  
  # # 感染最大范围?需要足够大
  N <- 1400000000
  
  # 以多少天窗口计算R0
  # WatchingWindow <- 12
  
  init <- c(S = N - Infected[1], I = Infected[1], R = 0)
  
  # 建模函数
  SIR <- function(time, state, parameters) {
    par <- as.list(c(state, parameters))
    with(par, {
      dS <- -beta/N * I * S
      dI <- beta/N * I * S - gamma * I
      dR <- gamma * I
      list(c(dS, dI, dR))
    })
  }
  
  # 首先写一个计算 R0 的函数：
  R0 <- function(Day){
    
    WindowInfected <- Infected[(Day - WatchingWindow):Day]
    Days <- 1:(length(WindowInfected))
    init <- c(S = N - WindowInfected[1], I = WindowInfected[1], R = 0)
    
    LOSS <- function(parameters) {
      names(parameters) <- c("beta", "gamma")
      out <- ode(y = init, times = Days, func = SIR, parms = parameters)
      fit <- out[ , 3]
      -R2(WindowInfected,fit)
      # MSE(WindowInfected, fit)
    }
    
    Opt = optim(c(0.5, 0.14), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 0.5))
    Opt_par = setNames(Opt$par, c("beta", "gamma"))
    as_tibble(list(Day = Day, 
                   R0 = Opt_par["beta"] / Opt_par["gamma"],
                   beta = Opt_par["beta"], 
                   gamma = Opt_par["gamma"], 
                   R2 = -Opt$value
                   # date = ymd("2020-01-11") + days(1 + dat - 1) - 1 # 1 - 1 | dat - 1
    )) 
  }
  return(map_df((WatchingWindow+1):Ndays, R0))
}

# ------------获取历史数据---------------
jsondata <- fromJSON('https://interface.sina.cn/news/wap/fymap2020_data.d.json')
# 当前时间
(times <- jsondata$data$times)
# "截至2月5日22时01分"

# 历史数据
historylist <- jsondata$data$historylist %>%
  mutate(
    date = ymd(paste0("2020.", date))
  ) %>%
  as_tibble() %>%
  type_convert()
historylist


Type <- "wuhan"

switch (Type,
        "wuhan" = {
          dat <- historylist %>% pull(wuhan_conNum) %>% rev()
        },
        "total" = {
          dat <- historylist %>% pull(cn_conNum) %>% rev()
        }
)
# -----------------计算R0---------------------

# 窗口期
WatchingWindow <- 12

res <- getTrendR0(dat, WatchingWindow = WatchingWindow)
res$date <- historylist$date[1 : (nrow(historylist) - WatchingWindow)] %>% rev()
print(res, n = 100)


ggplot(filter(res, R2 > 0.7), aes(x = date)) +
  geom_line(aes(y = R0), color = "red", stat = "identity") +
  geom_text(aes(y = R0, label = round(R0, 2))) +
  scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
  labs(title = paste(Type,"R0:", WatchingWindow,"days trend"))


# qplot(date, R0, data = filter(res, R2 > 0.7), 
#       geom = "line", 
#       main = paste(Type,"R0:", WatchingWindow,"days trend"), alpha = 0.1) +
#   geom_text(aes(date, R0, label = round(R0, 2))) +
#   scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
#   theme(legend.position = "none")

# -----------------获取历史数据---------------
# remotes::install_github("GuangchuangYu/nCov2019")

library(nCov2019)
rawDat <- load_nCov2019()$data 

# 窗口期
WatchingWindow <- 12

# NA处理成0
rawDat <- rawDat %>%
  as_tibble() %>%
  select(province, city, time, confirm, dead, heal, cum_confirm, cum_heal, cum_dead) %>%
  mutate_at(vars(confirm, dead, heal, cum_confirm, cum_heal, cum_dead), as.integer) %>%
  mutate_if(is.numeric, list(~ifelse(is.na(.), 0L, .))) 

st <- Sys.time()
resRaw <- filter(rawDat) %>%
  group_by(province, city) %>%
  arrange(time) %>%
  summarise(res = list(tryCatch(getTrendR0(cum_confirm, WatchingWindow = WatchingWindow), error = function(e) NULL)), dates = list(time))
et <- Sys.time()
print(et - st)

res <- lapply(1:nrow(resRaw), function(x){
  # message(x)
  temp <- resRaw[x,]
  if (length(unlist(temp$res)) == 0) {
    return(data.frame())
  } else {
    dat <- temp$res[[1]] %>%
      mutate(province = temp$province,
             city = temp$city,
             date = temp$dates[[1]][(WatchingWindow + 1):length(temp$dates[[1]])]) %>%
      select(date, province, city, R0, beta, gamma, R2)
    return(dat)
  }
}) %>%
  bind_rows()

res_final <- left_join(res, rawDat, by = c("province" = "province", "city" = "city", "date" = "time"))
write_csv(res_final, "res_02_18_12.csv")

res_final <- read_csv("res_02_18_12.csv")

filter(res_final, R0 > 1.5, R2 > 0.8, date %in% as_date(c('2020-02-09')), province != city) %>% 
  arrange(desc(R0), province, city, date) %>%
  print(n = 50)


filter(res1, R0 > 1.5, R2 > 0.8, date %in% as_date(c('2020-02-07')), province == city) %>% 
  arrange(province, city, date) %>%
  print(n = 100)


filter(res1, R0 > 1.5, R2 > 0.8, date >= '2020-02-04', province == city) %>% 
  arrange(province, city, date) %>%
  print(n = 100)




filter(res1, province == "湖北", date == '2020-02-07') %>% 
  arrange(desc(R0)) %>%
  print(n = 100)

filter(res, province == "黑龙江", date >= '2020-02-04') %>% 
  print(n = 100)


filter(res_final, city == "武汉") %>% 
  print(n = 100)

# 
# filter(rawDat, city == "荆门") %>% 
#   print(n = 100)



# -----------------------------------------------------------------------------------------------

rawDat1 <- filter(rawDat, province == "湖北") %>%
  mutate(city = ifelse(city == "武汉", "武汉", "其他")) %>%
  group_by(province, city, time) %>%
  summarise_at(vars(confirm:cum_dead), sum) 

st <- Sys.time()
resRaw1 <- filter(rawDat1) %>%
  group_by(province, city) %>%
  arrange(time) %>%
  summarise(res = list(tryCatch(getTrendR0(cum_confirm, WatchingWindow = WatchingWindow), error = function(e) NULL)), dates = list(time))
et <- Sys.time()
print(et - st)

res1 <- lapply(1:nrow(resRaw1), function(x){
  # message(x)
  temp <- resRaw1[x,]
  if (length(unlist(temp$res)) == 0) {
    return(data.frame())
  } else {
    dat <- temp$res[[1]] %>%
      mutate(province = temp$province,
             city = temp$city,
             date = temp$dates[[1]][(WatchingWindow + 1):length(temp$dates[[1]])]) %>%
      select(date, province, city, R0, beta, gamma, R2)
    return(dat)
  }
}) %>%
  bind_rows()

res_final1 <- left_join(res1, rawDat1, by = c("province" = "province", "city" = "city", "date" = "time"))
write_csv(res_final1, "res_hb_02_09_12.csv")


rawDat2 <- mutate(rawDat, province = "全国", city = "全国") %>%
  group_by(province, city, time) %>%
  summarise_at(vars(confirm:cum_dead), sum) 

st <- Sys.time()
resRaw2 <- group_by(rawDat2, province, city) %>%
  arrange(time) %>%
  summarise(res = list(tryCatch(getTrendR0(cum_confirm, WatchingWindow = WatchingWindow), error = function(e) NULL)), dates = list(time))
et <- Sys.time()
print(et - st)

res2 <- lapply(1:nrow(resRaw2), function(x){
  # message(x)
  temp <- resRaw2[x,]
  if (length(unlist(temp$res)) == 0) {
    return(data.frame())
  } else {
    dat <- temp$res[[1]] %>%
      mutate(province = temp$province,
             city = temp$city,
             date = temp$dates[[1]][(WatchingWindow + 1):length(temp$dates[[1]])]) %>%
      select(date, province, city, R0, beta, gamma, R2)
    return(dat)
  }
}) %>%
  bind_rows()

res_final2 <- left_join(res2, rawDat2, by = c("province" = "province", "city" = "city", "date" = "time"))
write_csv(res_final2, "res_total_02_09_12.csv")


# -----------------计算R0 长短期曲线绘图---------------------


resList <- lapply(c(5,12), function(WatchingWindow){
  res <- getTrendR0(dat, WatchingWindow = WatchingWindow)
  res$date <- historylist$date[1 : (nrow(historylist) - WatchingWindow)] %>% rev()
  res
}) 

res <- left_join(resList[[1]], resList[[2]], by = "date", suffix = c("_long", "_short")) %>%
  select(date, R0_short, R0_long, beta_short, beta_long, gamma_short, gamma_long, R2_short, R2_long)

ggplot(res, aes(x = date)) +
  geom_line(aes(y = R0_short), color = "red") +
  geom_line(aes(y = R0_long), color = "blue") 


lapply(1:10, function(x){
  qplot(date, R0, data = filter(resList[[x]], R2 > 0.7), 
        geom = "line", 
        main = paste(Type,"R0 Trend"), alpha = 0.1) +
    geom_text(aes(date, R0, label = round(R0, 2))) +
    scale_x_date(date_breaks = "1 day", date_labels = "%m-%d") +
    theme(legend.position = "none")
})









resRaw <- filter(rawDat) %>%
  group_by(province, city) %>%
  arrange(time) %>%
  summarise(
    R0 = list(tryCatch(getTrendR0(cum_confirm)$data[["R0"]], error = function(e) NULL)),
    date = list(tryCatch(getTrendR0(cum_confirm)$data[["date"]], error = function(e) NULL)))

res <- lapply(1:nrow(resRaw), function(x){
    temp <- resRaw[x,]
    if (length(unlist(temp$R0)) == 0) {
      return(data.frame())
    } else {
      return(data.frame(
        province = rep(temp$province, length(temp$R0[[1]])),
        city = rep(temp$city, length(temp$R0[[1]])),
        R0 = temp$R0[[1]],
        beta = temp$beta[[1]],
        gamma = temp$gamma[[1]],
        value = temp$value[[1]],
        date = temp$date[[1]]))
    }

  }) %>%
  bind_rows()




lapply(1:nrow(res), function(x){
  temp1 <- res[x,]
  data.frame(city = rep(temp1$city, length(temp1$R0[[1]])),
             R0 = temp1$R0[[1]],
             date = temp1$date[[1]])
}) %>%
  bind_rows()



temp <- filter(rawDat, city == "恩施") %>%
  group_by(province, city) %>%
  arrange(time)
  
new_sum <- function(x) {
  res <- list()
  res$sum <- sum(x)
  res$message <- "ok"
  return(res)
}

temp1 <- summarise(temp, n = list(getTrendR0(cum_confirm)))

# ---------------------------------------------------------------
# 读取R0数据
res_final2 <- read_csv("res_02_09_12.csv")

# 取城市名称前两位做匹配
Chinamap_cities_sf <- map_prefecture_city() %>%
  mutate(c2 = str_sub(name, 1, 2))

# 定义分段显示
plot_breaks <-c(0, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3)
plot_labels <-c("0", "0.6-0.8", "0.8-1", "1-1.2", "1.2-1.4", "1.4-1.6",
             "1.6-1.8", "1.8-2.0", "2.0-2.2", "2.2-2.4", "2.4+")

epidemic_df <- filter(res_final2, province != city, province != "全国") %>%
  mutate(R02 = as.character(cut(R0, breaks = plot_breaks,
                   labels = plot_labels, 
                   include.lowest = TRUE,
                   right = FALSE, 
                   ordered_result = TRUE))) %>%
  mutate(c2 = str_sub(city, 1, 2)) 



# 绘图存储
the_days <- filter(epidemic_df, date >= '2020-01-29') %>% pull(date) %>% unique() %>% sort()

# create temporary document
dir.create(dir1 <- file.path("ncov_images"))

for (i in seq_len(length(the_days))) {
  
  # 绘图
  epidemic_the_date <- filter(epidemic_df, date == the_days[i]) %>%
    right_join(Chinamap_cities_sf, by = "c2") %>%
    mutate(R0 = ifelse(is.na(R0), 0, R0),
           R02 = ifelse(is.na(R02), "0", R02)) %>%
    select(city, R0, R02, adcode, name, level, geometry)
  
  gg_epidemic <- ggplot(epidemic_the_date, color="grey") +
    geom_sf(aes(fill = R02, geometry = geometry)) +
    coord_sf() +
    scale_fill_brewer(palette = "YlOrRd", direction = 1) +
    guides(fill = guide_legend(title = "R0", reverse = T)) +
    labs(title = "2019-ncov疫情数据可视化:12日窗口拟合R0变化趋势",
         subtitle = the_days[i],
         caption = "数据来源:GuangchuangYu/nCov2019") +
    theme(
      # 标题
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 20, color = "red"),
      plot.caption = element_text(face = "bold", hjust = 1, color = "blue"),
      # 图例
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(face = "bold", color = "black"),
      legend.background = element_rect(colour = "black"),
      legend.key = element_rect(fill = NA), # 图例箱体无背景
      legend.position = c(0.85, 0.2),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      # 绘图面板
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", linetype = "solid", size = 1, fill = NA)
    )
  # save picture
  ggsave(filename = paste0(the_days[i], ".png"),
         plot = gg_epidemic, path = dir1,
         width = 20, height = 20, units = "cm")
}

# 制作动图
path_pre <- "./"
animate_epidemic <- image_animate(image = image_read(path = paste0(dir1, "/", the_days, ".png")), fps = 20)
anim_save(filename = "疫情地图可视化_各地12日拟合R0_动态图.gif", animation = animate_epidemic, path = path_pre)

unlink(dir1)











epidemic_plot <- epidemic_df %>%
  right_join(Chinamap_cities_sf, by = "c2") %>%
  mutate(R0 = ifelse(is.na(R0), 0, R0),
         R02 = ifelse(is.na(R02), "0", R02)) %>%
  select(date, city, R0, R02, adcode, name, level, geometry)

gg_epidemic <- ggplot(epidemic_plot, color = "grey") +
  facet_wrap(~date) +
  geom_sf(aes(fill = R02, geometry = geometry)) +
  coord_sf() +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  guides(fill = guide_legend(title = "R0", reverse = T)) +
  labs(title = "2019-ncov疫情数据可视化:12日窗口拟合R0变化趋势",
       subtitle = the_days[i],
       caption = "数据来源:GuangchuangYu/nCov2019") +
  theme(
    # 标题
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 20, color = "red"),
    plot.caption = element_text(face = "bold", hjust = 1, color = "blue"),
    # 图例
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.background = element_rect(colour = "black"),
    legend.key = element_rect(fill = NA), # 图例箱体无背景
    legend.position = c(0.85, 0.2),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    # 绘图面板
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linetype = "solid", size = 1, fill = NA))





the_days[1:6]

ggplot(filter(epidemic_plot, date %in% the_days[7:12]), color = "grey") +
  facet_wrap(~date) +
  geom_sf(aes(fill = R02, geometry = geometry)) +
  coord_sf() +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  guides(fill = guide_legend(title = "R0", reverse = T)) +
  labs(title = "2019-ncov疫情数据可视化:12日窗口拟合R0变化趋势",
       subtitle = the_days[i],
       caption = "数据来源:GuangchuangYu/nCov2019")




library(R0)
data("Germany.1918")

dat1 <- Germany.1918 %>% cumsum()

forged <- Germany.1918
forged[1:60] <- forged[1:60] * seq(0.1, 5, length.out = 60)

dat1_forged <- forged %>% as.integer() %>% cumsum()

res <- getTrendR0(Infected = as.numeric(dat1), WatchingWindow = 12)
res_forged <- getTrendR0(Infected = as.numeric(dat1_forged), WatchingWindow = 12)




