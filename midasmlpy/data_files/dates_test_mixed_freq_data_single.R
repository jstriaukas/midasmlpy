rm(list=ls())
setwd("~/Documents/GitHub/midasmlpy/R_test")

require(midasml)

### mixed_freq_data_single ### 
data(us_rgdp)
rgdp <- us_rgdp$rgdp
cfnai <- us_rgdp$cfnai
data.refdate <- rgdp$date
data.x <- cfnai$cfnai
data.xdate <- cfnai$date
est.start <- as.Date("1990-01-01")
est.end <- as.Date("2002-03-01")
x.lag <- 12
horizon <- 1
disp.flag <- FALSE
write.table(data.refdate, file = "input_data_refdate.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(data.x, file = "input_data_x.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(x.lag, file = "input_x_lag.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(horizon, file = "input_horizon.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(est.start, file = "input_est_start.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(est.end, file = "input_est_end.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

out <- mixed_freq_data_single(data.refdate, data.x, data.xdate, x.lag = 12, horizon = 1,
                       est.start, est.end, disp.flag = FALSE)

est.refdate = out$est.refdate
est.x = out$est.x
est.xdate = out$est.xdate
out.refdate = out$out.refdate
out.x = out$out.x
out.xdate = out$out.xdate
x.lag = out$x.lag
min.date = out$min.date
max.date = out$max.date

write.table(est.refdate, file = "output_est_refdate.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(est.x, file = "output_est_x.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(est.xdate, file = "output_est_xdate.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.refdate, file = "output_out_refdate.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.x, file = "output_out_x.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.xdate, file = "output_out_xdate.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(x.lag, file = "output_x_lag.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(min.date, file = "output_min_date.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(max_date, file = "output_max_date.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
