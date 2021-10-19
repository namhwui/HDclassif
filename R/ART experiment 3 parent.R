# ART experiment 3
# compare ART against existing dimension selection methods

# intel xeon gold 6150 @ 2.70ghz (36 cores, 30 cores used)

library(parallel)
library(pbapply)
#library(intrinsicDimension)


cl <- makeCluster(30)
reps <- 500
ng <- seq(200, 1000, by = 100)

set.seed(1) # for Q generation
ev_true <- matrix(nrow = 3, ncol = 100)
ev_true[1, ] <- c(seq(10, 9, length.out = 5), rep(1, 95))
ev_true[2, ] <- c(seq(5, 3, length.out = 15), rep(0.5, 85))
ev_true[3, ] <- c(seq(15, 13, length.out = 10), rep(5, 90))
Q <- list(pracma::randortho(100, "orthonormal"),
		  pracma::randortho(100, "orthonormal"),
		  pracma::randortho(100, "orthonormal"))
sigma <- lapply(1:3, function(x) {
	Q[[x]] %*% diag(ev_true[x, ]) %*% t(Q[[x]])
})

clusterExport(cl, c("ng", "sigma"))


d_art <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  source("controls.R", local = T)
	  source("hdc.R", local = T)
	  source("hddc.R", local = T)
	  source("ART functions.R", local = T)
	  
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- hddc(data, K = 3, model = "AkjBkQkDk", 
							d_select = "art", threshold = 0.0001, 
							init = "vector", init.vector = label, 
							itermax = 1, show = F)$d
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_art, "ART experiment 3-1 art.rds")


d_art2 <- lapply(1:length(ng), function(jj) {
  pbsapply(1:reps, function(ii) {
    source("controls.R", local = T)
    source("hdc.R", local = T)
    source("hddc.R", local = T)
    source("ART functions.R", local = T)
    
    set.seed(ii)
    label <- sort(rep(1:3, ng[jj]))
    data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
                  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
                  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
    
    m <- system.time(
      suppressWarnings({
        val <- hddc(data, K = 3, model = "AkjBkQkDk", 
                    d_select = "art", threshold = 0.1, 
                    init = "vector", init.vector = label, 
                    itermax = 1, show = F)$d
      })
    )
    c(m[3], val)
  }, cl = cl)
})
saveRDS(d_art2, "ART experiment 3-1 art2.rds")


d_bic <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  source("controls.R", local = T)
	  source("hdc.R", local = T)
	  source("hddc.R", local = T)
	  
	  
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- hddc(data, K = 3, model = "AkjBkQkDk", 
							d_select = "bic", 
							init = "vector", init.vector = label, 
							itermax = 1, show = F)$d
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_bic, "ART experiment 3-1 bic.rds")

d_scree <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  source("controls.R", local = T)
	  source("hdc.R", local = T)
	  source("hddc.R", local = T)
	  
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- hddc(data, K = 3, model = "AkjBkQkDk", 
							d_select = "cattell", threshold = 0.2, 
							init = "vector", init.vector = label, 
							itermax = 1, show = F)$d
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_scree, "ART experiment 3-1 scree.rds")


d_scree2 <- lapply(1:length(ng), function(jj) {
  pbsapply(1:reps, function(ii) {
    source("controls.R", local = T)
    source("hdc.R", local = T)
    source("hddc.R", local = T)
    
    set.seed(ii)
    label <- sort(rep(1:3, ng[jj]))
    data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
                  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
                  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
    
    m <- system.time(
      suppressWarnings({
        val <- hddc(data, K = 3, model = "AkjBkQkDk", 
                    d_select = "cattell", threshold = 0.001, 
                    init = "vector", init.vector = label, 
                    itermax = 1, show = F)$d
      })
    )
    c(m[3], val)
  }, cl = cl)
})
saveRDS(d_scree2, "ART experiment 3-1 scree2.rds")



d_lpca <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  library(intrinsicDimension)
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- c(pcaLocalDimEst(data[label == 1, ], "fan", verbose = F)$dim.est,
						     pcaLocalDimEst(data[label == 2, ], "fan", verbose = F)$dim.est,
						     pcaLocalDimEst(data[label == 3, ], "fan", verbose = F)$dim.est)
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_lpca, "ART experiment 3-1 lpca.rds")

d_optm <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  library(intrinsicDimension)
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- round(c(median(pcaOtpmPointwiseDimEst(data[label == 1, ], 100)$dim.est),
							   median(pcaOtpmPointwiseDimEst(data[label == 2, ], 100)$dim.est),
							   median(pcaOtpmPointwiseDimEst(data[label == 3, ], 100)$dim.est)))
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_optm, "ART experiment 3-1 optm.rds")

d_ess <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  library(intrinsicDimension)
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- round(c(essLocalDimEst(data[label == 1, ])$dim.est,
							   essLocalDimEst(data[label == 2, ])$dim.est,
							   essLocalDimEst(data[label == 3, ])$dim.est))
			})
		)
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_ess, "ART experiment 3-1 ess.rds")

d_knn <- lapply(1:length(ng), function(jj) {
	pbsapply(1:reps, function(ii) {
	  library(intrinsicDimension)
		set.seed(ii)
		label <- sort(rep(1:3, ng[jj]))
		data <- rbind(Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[1]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[2]]),
					  Rfast::rmvnorm(ng[jj], rep(0, 100), sigma[[3]]))
		
		m <- system.time(
			suppressWarnings({
				val <- round(c(knnDimEst(data[label == 1, ], 2, M = 5, gamma = 2, ps = seq(max(3, round(ng[jj]/2)), ng[jj]-1, by=3))$dim.est,
							   knnDimEst(data[label == 2, ], 2, M = 5, gamma = 2, ps = seq(max(3, round(ng[jj]/2)), ng[jj]-1, by=3))$dim.est,
							   knnDimEst(data[label == 3, ], 2, M = 5, gamma = 2, ps = seq(max(3, round(ng[jj]/2)), ng[jj]-1, by=3))$dim.est))}))
		c(m[3], val)
	}, cl = cl)
})
saveRDS(d_knn, "ART experiment 3-1 knn.rds")