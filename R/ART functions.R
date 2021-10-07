# IF components don't interact: test 1 component at a time
# IF components interact:
#    IF common dimension: test all at the same time (harmonic p value)
#    IF diff dimension: start with all 1, find individual components to raise dim (instead of brute force grid search)

art_AkjBk <- function(v, nk, x) {
	p <- length(v)
	a <- v[1:x]
	b <- mean(v[-(1:x)])
	ratio <- v / c(a, rep(b, p - x))
	#print(range(v))
	-p*(nk-1)*log(nk*exp(1)/(nk-1)) - (nk-1)*sum(log(ratio)) + nk*sum(ratio)
}

art_AkBk <- function(v, nk, x) {
	p <- length(v)
	a <- mean(v[1:x])
	b <- mean(v[-(1:x)])
	ratio <- v / c(rep(a, x), rep(b, p - x))
	-p*(nk-1)*log(nk*exp(1)/(nk-1)) - (nk-1)*sum(log(ratio)) + nk*sum(ratio)
}

art_ABk <- function(ev, n) {
	prop <- n / sum(n)
	K <- nrow(ev)
	p <- ncol(ev)
	function(d) {
		a_weighted <- sapply(1:K, function(ii) {
			prop[ii] * sum(ev[ii, 1:d[ii]])
		})
		a <- sum(a_weighted) / sum(prop * d)
		bk <- sapply(1:K, function(ii) {
			mean(ev[ii, -(1:d[ii])])
		})
		
		sapply(1:K, function(ii) {
			ratio <- ev[ii, ] / c(rep(a, d[ii]), rep(bk[ii], p - d[ii]))
			-p*(n[ii]-1)*log(n[ii]*exp(1)/(n[ii]-1)) - (n[ii]-1)*sum(log(ratio)) + n[ii]*sum(ratio)
		})
	}
}

art_AkjB <- function(ev, n) {
	prop <- n / sum(n)
	K <- nrow(ev)
	p <- ncol(ev)
	function(d) {
		b_weighted <- sapply(1:K, function(ii) {
			prop[ii] * sum(ev[ii, -(1:d[ii])])
		})
		b <- sum(b_weighted) / (p - sum(prop * d))
		
		sapply(1:K, function(ii) {
			ratio <- ev[ii, ] / c(ev[ii, 1:(d[ii])], rep(b, p-d[ii]))
			-p*(n[ii]-1)*log(n[ii]*exp(1)/(n[ii]-1)) - (n[ii]-1)*sum(log(ratio)) + n[ii]*sum(ratio)
		})
	}
}


art_AkB <- function(ev, n) {
	prop <- n / sum(n)
	K <- nrow(ev)
	p <- ncol(ev)
	function(d) {
		b_weighted <- sapply(1:K, function(ii) {
			prop[ii] * sum(ev[ii, -(1:d[ii])])
		})
		b <- sum(b_weighted) / (p - sum(prop * d))
		
		sapply(1:K, function(ii) {
			a <- mean(ev[ii, 1:(d[ii])])
			ratio <- ev[ii, ] / c(rep(a, d[ii]), rep(b, p-d[ii]))
			-p*(n[ii]-1)*log(n[ii]*exp(1)/(n[ii]-1)) - (n[ii]-1)*sum(log(ratio)) + n[ii]*sum(ratio)
		})
	}
}

art_AB <- function(ev, n) {
	prop <- n / sum(n)
	K <- nrow(ev)
	p <- ncol(ev)
	function(d) {
		a_weighted <- sapply(1:K, function(ii) {
			prop[ii] * sum(ev[ii, 1:d[ii]])
		})
		a <- sum(a_weighted) / sum(prop * d)
		
		b_weighted <- sapply(1:K, function(ii) {
			prop[ii] * sum(ev[ii, -(1:d[ii])])
		})
		b <- sum(b_weighted) / (p - sum(prop * d))
		
		sapply(1:K, function(ii) {
			ratio <- ev[ii, ] / c(rep(a, d[ii]), rep(b, p-d[ii]))
			-p*(n[ii]-1)*log(n[ii]*exp(1)/(n[ii]-1)) - (n[ii]-1)*sum(log(ratio)) + n[ii]*sum(ratio)
		})
	}
}

art_dim_indep <- function(ev, n, threshold, model, eps = 1e-6) {
	K <- nrow(ev)
	#p <- ncol(ev)
	#d_all <- 1:(p-1)
	val <- rep(1,K)
	
	if (model == "AKJBKQKDK") {
		for (ii in 1:K) {
			ev_temp <- ev[ii, ][ev[ii, ] > eps]
			p <- length(ev_temp)
			d_all <- 1:(p-1)
			df_all <- p*(p+1)/2 - d_all*(p-(d_all+1)/2) - d_all - 1
			#print(df_all)
			test <- sapply(d_all, function(x) {
				art_AkjBk(ev_temp, n[ii], x)
			})
			#print(test)
			val[ii] <- which(test <= qchisq(1 - threshold, df_all))[1]
			if (is.na(val[ii])) {
				val[ii] <- p-1
			}
		}
	} else if (model == "AKBKQKDK") {
		for (ii in 1:K) {
			ev_temp <- ev[ii, ][ev[ii, ] > eps]
			p <- length(ev_temp)
			d_all <- 1:(p-1)
			df_all <- p*(p+1)/2 - d_all*(p-(d_all+1)/2) - 2
			
			#print(df_all)
			test <- sapply(d_all, function(x) {
				art_AkBk(ev_temp, n[ii], x)
			})
			val[ii] <- which(test <= qchisq(1 - threshold, df_all))[1]
			if (is.na(val[ii])) {
				val[ii] <- p-1
			}
		}
	}
	val
}

art_dim_dep <- function(ev, n, threshold, model, eps = 1e-6) {
	K <- nrow(ev)
	#d_all <- 1:(p-1)
	d_init <- rep(1, K)
	continue <- T
	
	# remove excessively small eigenvalues
	# since dividing by them causes NaN
	d_max <- max(apply(ev > eps, 2, prod) * (1:ncol(ev)))
	ev <- ev[, 1:d_max]
	p <- ncol(ev)
	
	if (model == "ABKQKDK") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*d_init/sum(n*d_init)
			test <- art_ABk(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			#print(d_init)
			if (whm_pval < threshold) {
				names(pval) <- 1:K
				pval_min <- suppressWarnings(min(pval[d_init < p-1]))
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- as.numeric(names(which.min(pval[d_init < p-1])))[1]
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "AKJBQKDK") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - d_init - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AkjB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				names(pval) <- 1:K
				pval_min <- suppressWarnings(min(pval[d_init < p-1]))
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- as.numeric(names(which.min(pval[d_init < p-1])))[1]
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "AKBQKDK") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AkB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				names(pval) <- 1:K
				pval_min <- suppressWarnings(min(pval[d_init < p-1]))
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- as.numeric(names(which.min(pval[d_init < p-1])))[1]
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "ABQKDK") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - n*d_init/sum(n*d_init) - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				names(pval) <- 1:K
				pval_min <- suppressWarnings(min(pval[d_init < p-1]))
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- as.numeric(names(which.min(pval[d_init < p-1])))[1]
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
	}
}

art_comdim <- function(ev, n, threshold, model, eps = 1e-6) {
	K <- nrow(ev)
	d_init <- rep(1, K)
	continue <- T
	
	# remove excessively small eigenvalues
	# since dividing by them causes NaN
	d_max <- max(apply(ev > eps, 2, prod) * (1:ncol(ev)))
	ev <- ev[, 1:d_max]
	p <- ncol(ev)
	#d_all <- 1:(p-1)
	
	if (model == "AKJBKQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - d_init
			test <- sapply(1:K, function(ii) {
				art_AkjBk(ev[ii, ], n[ii], d_init[ii])
				})
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "AKBKQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 2
			test <- sapply(1:K, function(ii) {
				art_AkBk(ev[ii, ], n[ii], d_init[ii])
			})
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "ABKQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*d_init/sum(n*d_init)
			test <- art_ABk(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "AKJBQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - d_init - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AkjB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "AKBQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AkB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "ABQKD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - n*d_init/sum(n*d_init) - n*(p-d_init)/sum(n*(p-d_init))
			test <- art_AB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				if (d_init[1] == p-1) {
					continue <- F
					return(d_init)
				} else {
					d_init <- d_init + 1
				}
			} else {
				continue <- F
			}
		}
		return(d_init)
	} else if (model == "AJBQD") {
		
	} else if (model == "ABQD") {
		
	}
}
art_main <- function(ev, n, threshold, graph = F, noise.ctrl = 1e-6, model) {
	
	# n = vector of component-wise sample size
	n <- round(n)
	d <- numeric(nrow(ev))
	if (model %in% c("AKJBKQKDK", "AKBKQKDK")) {
		# component-wise covariance don't interact
		d <- art_dim_indep(ev, n, threshold, model)
	} else if (model %in% c("ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK")) {
		d <- art_dim_dep(ev, n, threshold, model)
	} else if (model %in% c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD")) {
		d <- art_comdim(ev, n, threshold, model)
	} else if (model %in% c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AkBQkD", "ABQKD", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD")) {
		
	} else {
		stop("Unrecognised submodel.")
	}
	d
}


###################################
# Scree Test
###################################

scree_test <- function(ev, threshold) {
	if (is.null(dim(ev))) {
		ev <- matrix(ev, nrow = 1)
	}
	abs_diff <- apply(ev, 1, function(x) {
		x_diff <- abs(diff(x))
		x_diff / max(x_diff)
	})
	apply((abs_diff > threshold) * 1:(ncol(ev) - 1), 2, max)
	
}
