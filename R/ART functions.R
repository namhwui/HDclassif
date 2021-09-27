# IF components don't interact: test 1 component at a time
# IF components interact:
#    IF common dimension: test all at the same time (harmonic p value)
#    IF diff dimension: start with all 1, find individual components to raise dim (instead of brute force grid search)

art_AkjBk <- function(v, nk, x) {
	p <- length(v)
	a <- v[1:x]
	b <- mean(v[-(1:x)])
	ratio <- v / c(a, rep(b, p - x))
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

art_dim_indep <- function(ev, n, threshold, model) {
	K <- nrow(ev)
	p <- ncol(ev)
	d_all <- 1:(p-1)
	val <- rep(1,K)
	
	if (model == "AkjBkQkDk") {
		for (ii in 1:K) {
			df_all <- p*(p+1)/2 - d_all*(p-(d_all+1)/2) - d_all - 1
			print(df_all)
			test <- sapply(d_all, function(x) {
				art_AkjBk(ev[ii, ], n[ii], x)
			})
			val[ii] <- which(test <= qchisq(1 - threshold, df_all))[1]
		}
	} else if (model == "AkBkQkDk") {
		for (ii in 1:K) {
			df_all <- p*(p+1)/2 - d_all*(p-(d_all+1)/2) - 2
			print(df_all)
			test <- sapply(d_all, function(x) {
				art_AkBk(ev[ii, ], n[ii], x)
			})
			val[ii] <- which(test <= qchisq(1 - threshold, df_all))[1]
		}
	}
	val
}

art_dim_dep <- function(ev, n, threshold, model) {
	K <- nrow(ev)
	p <- ncol(ev)
	#d_all <- 1:(p-1)
	d_init <- rep(1, K)
	continue <- T
	
	if (model == "ABkQkDk") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*d_init/sum(n*d_init)
			test <- art_ABk(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				pval_min <- min(pval[d_init < p-1])
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- which(pval == pval_min)
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "AkjBQkDk") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - d_init - n*d_init/sum(n*d_init)
			test <- art_AkjB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				pval_min <- min(pval[d_init < p-1])
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- which(pval == pval_min)
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "AkBQkDk") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*d_init/sum(n*d_init)
			test <- art_AkB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				pval_min <- min(pval[d_init < p-1])
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- which(pval == pval_min)
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
		
	} else if (model == "ABQkDk") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 2*n*d_init/sum(n*d_init)
			test <- art_AB(ev, n)(d_init)
			pval <- pchisq(test, df_K, lower.tail = F)
			whm_pval <- 1 / sum(1 / pval)
			
			if (whm_pval < threshold) {
				pval_min <- min(pval[d_init < p-1])
				if (is.infinite(pval_min)) {
					continue <- F
					return(d_init)
				}
				raise <- which(pval == pval_min)
				d_init[raise] <- d_init[raise] + 1
			} else {
				continue <- F
			}
		}
		return(d_init)
	}
}

art_comdim <- function(ev, n, threshold, model) {
	K <- nrow(ev)
	p <- ncol(ev)
	d_all <- 1:(p-1)
	d_init <- rep(1, K)
	continue <- T
	
	if (model == "AkjBkQkD") {
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
	} else if (model == "AkBkQkD") {
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
	} else if (model == "ABkQkD") {
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
	} else if (model == "AkjBQkD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - d_init - n*d_init/sum(n*d_init)
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
	} else if (model == "AkBQkD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 1 - n*d_init/sum(n*d_init)
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
	} else if (model == "ABQkD") {
		while (continue) {
			df_K <- p*(p+1)/2 - d_init*(p-(d_init+1)/2) - 2*n*d_init/sum(n*d_init)
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
	} else if (model == "AjBQD") {
		
	} else if (model == "ABQD") {
		
	}
}
art_main <- function(ev, n, threshold, graph = F, noise.ctrl = 1e-6, model) {
	
	# n = vector of component-wise sample size
	n <- round(n)
	d <- numeric(nrow(ev))
	if (model %in% c("AkjBkQkDk", "AkBkQkDk")) {
		# component-wise covariance don't interact
		d <- art_dim_indep(ev, n, threshold, model)
	} else if (model %in% c("ABkQkDk", "AkjBQkDk", "AkBQkDk", "ABQkDk")) {
		d <- art_dim_dep(ev, n, threshold, model)
	} else if (model %in% c("AkjBkQkD", "AkBkQkD", "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD")) {
		d <- art_comdim(ev, n, threshold, model)
	} else if (model %in% c("AkjBkQkDk", "AkBkQkDk", "ABkQkDk", "AkjBQkDk", "AkBQkDk", "ABQkDk", "AkjBkQkD", "AkBkQkD", "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD")) {
		
	} else {
		error("Unrecognised submodel.")
	}
	d
}
