
# IC 20: ic_20,0.0,50.0,180,0,-10,40,10000
le_plot(ic20_lexp11, ic20_lexp24,ic20_lexp53,ic20_lexp64,ic20_lexp113, 
        "ic20", nrow(ic20_lexp11), 8)

# IC 23: ic_23,0.0,50.0,180,0,-10,40.001,10000
le_plot(ic23_lexp11, ic23_lexp24,ic23_lexp53,ic23_lexp64,ic23_lexp113, 
        "ic23", nrow(ic23_lexp11), 8)

# IC 24: ic_24,0.0,50.0,180,0,-10,40.01,10000
le_plot(ic24_lexp11, ic24_lexp24,ic24_lexp53,ic24_lexp64,ic24_lexp113, 
        "ic24", nrow(ic24_lexp11), 8)

# IC 25: ic_25,0.0,50.0,180,0,-10,40.1,10000
le_plot(ic25_lexp11, ic25_lexp24,ic25_lexp53,ic25_lexp64,ic25_lexp113, 
        "ic25", nrow(ic25_lexp11), 8)

#LEXP LINE PLOT
le_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymax) {
  plot(1,type='n',xlim=c(1,nrows),ylim=c(0,ymax),xlab='time step', ylab='exp', main=title)
  
  lines(d11$lexp, type='l', col="red", lwd=1.5)
  lines(d24$lexp, type='l', col="orange", lwd=1.5)
  lines(d53$lexp, type='l', col="green", lwd=1.5)
  lines(d64$lexp, type='l', col="blue", lwd=1.5)
  lines(d113$lexp, type='l', col="purple", lwd=1.5)
  
  legend("topright", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}










