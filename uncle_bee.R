library(R0)
library(nCov2019)
library(tidyverse)

mGT <- generation.time("gamma", c(8.5, 3.8))
#8.5+-3.8是SARS的，现在应该有新病毒的参数了
#把现在的参数，就是以上结果，直接放Y叔的Code里面。 发Bees.msu@gmail.com.
#Code是以下；这些code 需要Modify， 变成两个文件。 
places <- x[]$name #store names to places
places

# places <- c("Hunan")
mGT <- generation.time("gamma", c(8.5, 3.8))
dat <- load_nCov2019("en")$data

for (province in places) {
  sick = filter(dat, province == province) %>%
    group_by(time, province) %>%
    summarise(confirmed = sum(confirmed)) %>%
    arrange(time) %>%
    pull(confirmed)
  s.R <- estimate.R(sick, mGT,methods=c("EG"))
  pdf(paste0(province, ".pdf"))
  plotfit(s.R)  #save this to a file place.name.pdf
  save(s.R, file =paste0(province, ".rda"))
  dev.off()
}
