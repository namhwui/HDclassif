# ART experiment 1
# scree vs art threshold sensitivity

source("R/controls.R")
source("R/hdc.R")
source("R/hddc.R")
source("R/ART functions.R")



ev_true <- c(seq(5,4, length = 10), rep(0.5, 40))
n <- seq(100, 1000, by = 100)
reps <- 500

val_art1 <- array(dim = c(reps, 5, 10))
val_art2 <- array(dim = c(reps, 5, 10))
val_scree <- array(dim = c(reps, 5, 10))

set.seed(112358)
for (ii in 1:10) {
	print(paste("n is", ii))
	for (jj in 1:reps) {
		#print(jj)
		data <- Rfast::rmvnorm(n[ii], rep(0, 50), diag(ev_true))
		ev <- matrix(eigen(cov(data), only.values = T)$values, nrow = 1)
		
		val_art1[jj, , ii] <- sapply(c(0.0001, 0.001, 0.01, 0.05, 0.1), function(x) {
			art_main(ev, n, threshold = x, model = "AKJBKQKDK")
		})
		val_art2[jj, , ii] <- sapply(c(0.0001, 0.001, 0.01, 0.05, 0.1), function(x) {
			art_main(ev, n, threshold = x, model = "AKBKQKDK")
		})
		val_scree[jj, , ii] <- sapply(c(0.2, 0.1, 0.05, 0.01, 0.001), function(x) {
			scree_test(ev, x)
		})
	}
}


saveRDS(val_art1, "ART experiment 1 art1.rds")
saveRDS(val_art2, "ART experiment 1 art2.rds")
saveRDS(val_scree, "ART experiment 1 scree.rds")
