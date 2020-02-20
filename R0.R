require(R0)

data(Germany.1918)
mGT <- generation.time("gamma", c(3, 1.5))
estR0 <- estimate.R(Germany.1918, mGT, begin = 1, end = 27, methods = c("EG", "ML", "TD", "AR", "SB"), pop.size = 100000, nsim = 100)
attributes(estR0)

plot(Germany.1918, type = "line", col = "red")
lines(estR0$estimates$TD$pred, col = "blue")



diff_dt <- diff(dat,lag = 1)


mean_ci <- 5.2
max_ci <- 14
min_ci <- 1

sd_ci <- (max_ci - min_ci - mean_ci)/2



mGT <- generation.time("gamma", c(mean_ci, sd_ci))


estR0 <- estimate.R(diff_dt[5:31], 
                    mGT,
                    begin = 1, 
                    end = 27, 
                    methods = c("EG", "ML", "TD", "AR", "SB"), 
                    pop.size = 100000, 
                    nsim = 100)

estR0$estimates


estR0_sb <- est.R0.SB(diff_dt[5:31], mGT, time.step = 1)
estR0_sb
plot(estR0_sb$R, type = "line")


estR0_td <- est.R0.TD(diff_dt[5:32], mGT, time.step = 1, q = c(0.025, 0.975))
estR0_td$R
plot(estR0_td$R, type = "line")





res <- getTrendR0(dat[5:31], WatchingWindow = 12)
plot(res$R0, type = "line")



plot(diff_dt[5:32], type = "line", col = "red")
lines(estR0_td$pred, col = "blue")


plot(estR0_td$R, res$R0)



rgamma(n, shape, scale=1)
rgamma(3, 0.15)
