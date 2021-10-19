# ART experiment 2
# scree submodel identification

source("R/controls.R")
source("R/hdc.R")
source("R/hddc.R")
source("R/ART functions.R")

set.seed(112358)
reps <- 500
ng <- c(100, 500)
ev_true <- matrix(nrow = 2, ncol = 50)
ev_true[1, ] <- c(seq(5, 4, length.out = 5), rep(0.5, 45))
ev_true[2, ] <- c(seq(10, 8, length.out = 15), rep(1, 35))
Q <- list(pracma::randortho(ncol(ev_true), "orthonormal"), 
		  pracma::randortho(ncol(ev_true), "orthonormal"))
sigma0 <- list(Q[[1]] %*% diag(ev_true[1, ]) %*% t(Q[[1]]),
			   Q[[2]] %*% diag(ev_true[2, ]) %*% t(Q[[2]]))



val1 <- lapply(1:4, function(x) {
	array(dim = c(500, 2, 2))
})
for (ii in 1:2) {
	for (jj in 1:reps) {
		if (jj %in% seq(100, 500, by = 100)) {
			print(jj)
		}
		data <- list(Rfast::rmvnorm(ng[ii], rep(0, 50), sigma0[[1]]),
					 Rfast::rmvnorm(ng[ii], rep(0, 50), sigma0[[2]]))
		ev <- rbind(eigen(cov(data[[1]]), only.values = T)$values,
					eigen(cov(data[[2]]), only.values = T)$values)
		
		val1[[1]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "AKJBKQKDK")
		val1[[2]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "AKBKQKDK")
		val1[[3]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "ABKQKDK")
		val1[[4]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "ABQKDK")
	}
}



ev_true[1, ] <- c(seq(5, 4, length.out = 10), rep(0.5, 40))
ev_true[2, ] <- c(seq(10, 8, length.out = 10), rep(1, 40))
sigma0 <- list(Q[[1]] %*% diag(ev_true[1, ]) %*% t(Q[[1]]),
			   Q[[2]] %*% diag(ev_true[2, ]) %*% t(Q[[2]]))


val2 <- lapply(1:4, function(x) {
	array(dim = c(500, 2, 2))
})
for (ii in 1:2) {
	for (jj in 1:reps) {
		if (jj %in% seq(100, 500, by = 100)) {
			print(jj)
		}
		data <- list(Rfast::rmvnorm(ng[ii], rep(0, 50), sigma0[[1]]),
					 Rfast::rmvnorm(ng[ii], rep(0, 50), sigma0[[2]]))
		ev <- rbind(eigen(cov(data[[1]]), only.values = T)$values,
					eigen(cov(data[[2]]), only.values = T)$values)
		
		val2[[1]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "AKJBKQKDK")
		val2[[2]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "AKBKQKDK")
		val2[[3]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "ABKQKDK")
		val2[[4]][jj, , ii] <- art_main(ev, rep(ng[ii], 2), 0.001, model = "ABQKDK")
	}
}

saveRDS(val1, "ART experiment 2-1.rds")
saveRDS(val2, "ART experiment 2-2.rds")
