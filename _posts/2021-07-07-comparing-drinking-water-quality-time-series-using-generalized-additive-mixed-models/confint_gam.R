
# critical value for simultaneous interval:

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mgcv::mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

simul_critval <- function(model, N = 10000) {
  
  Vb <- vcov(model)
  
  pred <- predict(model, se.fit = TRUE)
  
  se.fit <- pred$se.fit
  
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  Cg <- predict(model, type = "lpmatrix")
  
  simDev <- Cg %*% t(BUdiff)
  
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  
  masd <- apply(absDev, 2, max)
  
  quantile(masd, prob = 0.95, type = 8)
  
}




