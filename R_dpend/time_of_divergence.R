time_of_divergence(ic10_pol11, ic10_pol113, 11)
time_of_divergence(ic10_pol24, ic10_pol113, 24)
time_of_divergence(ic10_pol53, ic10_pol113, 53)
time_of_divergence(ic10_pol64, ic10_pol113, 64)

time_of_divergence(ic11_pol11, ic11_pol113, 11)
time_of_divergence(ic11_pol24, ic11_pol113, 24)
time_of_divergence(ic11_pol53, ic11_pol113, 53)
time_of_divergence(ic11_pol64, ic11_pol113, 64)

time_of_divergence <- function(sm, lg, sm_nbits) {
  d0 <- (2^(-sm_nbits))
  mag1 <- sqrt((sm$th1)^2 + (sm$w1)^2 + (sm$th2)^2 + (sm$w2)^2)
  mag2 <- sqrt((lg$th1)^2 + (lg$w1)^2 + (lg$th2)^2 + (lg$w2)^2)
  
  mag_diff <- abs(mag2 - mag1)
  
  for (i in 1:nrow(sm)) {
    if (mag_diff[i] >  2*d0) {
      print(sm[i,]) 
      print(lg[i,])
      cat("diff: ", mag_diff[i])
      cat(" mEps: ", d0)
      stop()
    }
  }
}
