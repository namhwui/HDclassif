# ART experiment
# real data


source("R/controls.R")
source("R/hdc.R")
source("R/hddc.R")
source("R/ART functions.R")
library(e1071)
library(mclust)

data0 <- na.omit(read.csv("data/bankrupt.csv"))
label <- data0[, 1]
data <- scale(data0[, -c(1, 95)])

set.seed(42)
suppressWarnings({model1 <- hddc(data, K = 1:6, model = 1:6, d_select = "art", threshold = 0.0001, itermax = 500, show = F, kmeans.control = list(nstart = 100, iter.max = 20))}) 
suppressWarnings({model2 <- hddc(data, K = 1:6, model = 1:6, d_select = "cattell", itermax = 500, show = F, kmeans.control = list(nstart = 100, iter.max = 20))})
suppressWarnings({model3 <- hddc(data, K = 1:6, model = 1:6, d_select = "cattell", threshold = 0.001, itermax = 500, show = F, kmeans.control = list(nstart = 100, iter.max = 20))})
suppressWarnings({model4 <- hddc(data, K = 1:6, model = 1:6, d_select = "bic", itermax = 500, show = F, kmeans.control = list(nstart = 100, iter.max = 20))})

saveRDS(list(model1, model2, model3, model4), "ART experiment real data.rds")


