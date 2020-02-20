checkpoint::setSnapshot(Sys.Date() - 2)
remotes::install_github("GuangchuangYu/nCov2019")


library(nCov2019)
library(tidyverse)
# x <- get_nCov2019()
x <- load_nCov2019()

dat <- as_tibble(x$data) %>%
  select(province, city, time, confirm, dead, heal, cum_confirm, cum_heal, cum_dead) %>%
  mutate_at(vars(confirm, dead, heal, cum_confirm, cum_heal, cum_dead), funs(as.integer(.))) %>%
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0L, .))) 


filter(dat, city == "北京") %>% print(n = 100, width = 1000)






test <- x['湖北',] %>%
  as_tibble() 

filter(test, city == "武汉") %>% fill

map(test, function(x){
  x[is.na(x)] <- 0
})



head(test) %>%
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .)))


