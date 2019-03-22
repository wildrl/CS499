time_of_divergence(ic2_mag11, ic2_mag113, 11)
time_of_divergence(ic2_mag24, ic2_mag113, 24)
time_of_divergence(ic2_mag53, ic2_mag113, 53)
time_of_divergence(ic2_mag64, ic2_mag113, 64)


time_of_divergence <- function(sm, lg, sm_nbits) {
  d0 <- (2^(-sm_nbits))
  mag_sm <- sm$magnitude
  mag_lg <- lg$magnitude
  
  mag_diff <- abs(mag_sm - mag_lg)
  
  for (i in 1:nrow(sm)) {
    if (mag_diff[i] >  .000001) {
      print(sm[i,]) 
      print(lg[i,])
      cat("diff: ", mag_diff[i])
      cat(" mEps: ", d0)
      stop()
    }
  }
  
}


