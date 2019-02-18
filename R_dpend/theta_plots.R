#TH IC15
le_plot_ic15 <- function() {
  plot(1,type='n',xlim=c(1,200),ylim=c(0,8),xlab='time step', ylab='exp', main="ic_15")
  
  lines(`ic15_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic15_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic15_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic15_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic15_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(150, 7, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}
