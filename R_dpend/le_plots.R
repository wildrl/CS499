le_plot_ic3()
le_plot_ic4()
le_plot_ic5()
le_plot_ic6()
le_plot_ic7()

le_plot_ic14()
le_plot_ic15()
le_plot_ic16()



#LEXP LINE PLOT: ic_3,0.0,10.0,180.45,0,16.25,0,4000
le_plot_ic3 <- function() {
  plot(1,type='n',xlim=c(1,400),ylim=c(0,5),xlab='time step', ylab='exp', main="ic_3")
  
  lines(`ic3_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic3_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic3_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic3_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic3_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(1, 5, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}

#LEXP LINE PLOT: ic_4,0.0,10.0,180.45,0,16.00,0,4000
le_plot_ic4 <- function() {
  plot(1,type='n',xlim=c(1,400),ylim=c(0,5),xlab='time step', ylab='exp', main="ic_4")
  
  lines(`ic4_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic4_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic4_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic4_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic4_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(1, 5, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}

#LEXP LINE PLOT: ic_5,0.0,10.0,180.45,0,16.001,0,4000
le_plot_ic5 <- function() {
  plot(1,type='n',xlim=c(1,400),ylim=c(0,5),xlab='time step', ylab='exp', main="ic_5")
  
  lines(`ic5_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic5_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic5_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic5_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic5_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(1, 5, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}

#LEXP LINE PLOT: ic_6,0.0,10.0,180.45,0,16.01,0,4000
le_plot_ic6 <- function() {
  plot(1,type='n',xlim=c(1,400),ylim=c(0,5),xlab='time step', ylab='exp', main="ic_6")
  
  lines(`ic6_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic6_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic6_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic6_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic6_lexp113`$lexp, type='l', col="purple", lwd=1.5) 
  
  legend(1, 5, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}

#LEXP LINE PLOT: ic_7,0.0,10.0,180.45,0,16.1,0,4000
le_plot_ic7 <- function() {
  plot(1,type='n',xlim=c(1,400),ylim=c(0,5),xlab='time step', ylab='exp', main="ic_7")
  
  lines(`ic7_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic7_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic7_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic7_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic7_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(1, 5, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}

#LEXP IC14
le_plot_ic14 <- function() {
  plot(1,type='n',xlim=c(1,200),ylim=c(0,8),xlab='time step', ylab='exp', main="ic_14")
  
  lines(`ic14_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic14_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic14_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic14_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic14_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(150, 7, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}


#LEXP IC15
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

#LEXP IC16
le_plot_ic16 <- function() {
  plot(1,type='n',xlim=c(1,200),ylim=c(0,8),xlab='time step', ylab='exp', main="ic_16")
  
  lines(`ic16_lexp11`$lexp, type='l', col="red", lwd=1.5)
  lines(`ic16_lexp24`$lexp, type='l', col="orange", lwd=1.5)
  lines(`ic16_lexp53`$lexp, type='l', col="yellow", lwd=1.5)
  lines(`ic16_lexp64`$lexp, type='l', col="blue", lwd=1.5)
  lines(`ic16_lexp113`$lexp, type='l', col="purple", lwd=1.5)
  
  legend(150, 7, legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}



