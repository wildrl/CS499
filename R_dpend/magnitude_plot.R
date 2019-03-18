
mag_plot(ic10_pol11, ic10_pol24, ic10_pol53, ic10_pol64, ic10_pol113, 
        "ic10 magnitude", nrow(ic10_pol11), 15)

w_plot(ic10_pol11$w1, ic10_pol24$w1, ic10_pol53$w1, ic10_pol64$w1, ic10_pol113$w1, 
      "ic10 w1", nrow(ic10_pol11), -10, 10) 
w_plot(ic10_pol11$w2, ic10_pol24$w2, ic10_pol53$w2, ic10_pol64$w2, ic10_pol113$w2, 
      "ic25 w2", nrow(ic10_pol11), -10, 10)

sum_w_plot(ic10_pol11, ic10_pol24, ic10_pol53, ic10_pol64, ic10_pol113, 
           "ic10 w1+w2", nrow(ic10_pol11), -20, 30)


mag_plot(ic11_pol11, ic11_pol24, ic11_pol53, ic11_pol64, ic11_pol113, 
         "ic11 magnitude", nrow(ic11_pol11), 15)



mag_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymax) {
  plot(1,type='n',xlim=c(0,100000),ylim=c(0,ymax),xlab='time step', ylab='magnitude', main=title)
  
  y1 <- sqrt(d11$th1^2 + d11$w1^2 + d11$th2^2 + d11$w2^2)
  y2 <- sqrt(d24$th1^2 + d24$w1^2 + d24$th2^2 + d24$w2^2)
  y3 <- sqrt(d53$th1^2 + d53$w1^2 + d53$th2^2 + d53$w2^2)
  y4 <- sqrt(d64$th1^2 + d64$w1^2 + d64$th2^2 + d64$w2^2)
  y5 <- sqrt(d113$th1^2 + d113$w1^2 + d113$th2^2 + d113$w2^2)
 
  lines(y1, type='l', col="red", lwd=1.5)
  lines(y2, type='l', col="orange", lwd=1.5)
  lines(y3, type='l', col="green", lwd=1.5)
  lines(y4, type='l', col="blue", lwd=1.5)
  lines(y5, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}


sum_w_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymin, ymax) {
  plot(1,type='n',xlim=c(1,nrows),ylim=c(ymin,ymax),xlab='time step', ylab='exp', main=title)
  
  y1 <- d11$w1 + d11$w2
  y2 <- d24$w1 + d24$w2
  y3 <- d53$w1 + d53$w2
  y4 <- d64$w1 + d64$w2
  y5 <- d113$w1 + d113$w2
  
  lines(y1, type='l', col="red", lwd=1.5)
  lines(y2, type='l', col="orange", lwd=1.5)
  lines(y3, type='l', col="green", lwd=1.5)
  lines(y4, type='l', col="blue", lwd=1.5)
  lines(y5, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

w_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymin, ymax) {
  plot(1,type='n',xlim=c(1,nrows),ylim=c(ymin,ymax),xlab='time step', ylab='exp', main=title)
  
  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  lines(d113, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

