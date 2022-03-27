setwd("~/Documents/GitHub/midasmlpy/R_test")
#install.packages("midasml")
library(midasml)

set.seed(123)
x = matrix(rnorm(100 * 20), 100, 20)
beta = c(5,4,3,2,1,rep(0, times = 15))
y = x%*%beta + rnorm(100)
gindex = sort(rep(1:4,times=5))
gamma = 0.5
standardize = FALSE
intercept = TRUE
seed = 1234
K = 99

write.table(x, file = "input_x.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(y, file = "input_y.txt", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

solution <- tscv.sglfit(x = x, y = y, gindex = gindex, gamma = gamma, standardize = standardize, intercept = intercept, seed = seed, K = K)
save(x, y, gindex, gamma, standardize, intercept, seed, K, file = 'input.RData')
save(solution, file = 'output.RData')


