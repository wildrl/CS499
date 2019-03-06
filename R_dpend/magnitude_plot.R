# IC 25: ic_25,0.0,50.0,180,0,-10,40.1,10000
mag_plot(ic25_pol11, ic25_pol24,ic25_pol53,ic25_pol64,ic25_pol113, 
        "ic25 magnitude", nrow(ic25_pol11), 300)

w_plot(ic25_pol11$w1, ic25_pol24$w1,ic25_pol53$w1,ic25_pol64$w1,ic25_pol113$w1, 
         "ic25 w1", nrow(ic25_pol11), -20, 30)
w_plot(ic25_pol11$w2, ic25_pol24$w2,ic25_pol53$w2,ic25_pol64$w2,ic25_pol113$w2, 
       "ic25 w2", nrow(ic25_pol11), -20, 30)

sum_w_plot(ic25_pol11, ic25_pol24,ic25_pol53,ic25_pol64,ic25_pol113, 
           "ic25 w1+w2", nrow(ic25_pol11), -20, 30)


#LEXP LINE PLOT
mag_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymax) {
  plot(1,type='n',xlim=c(1,nrows),ylim=c(0,ymax),xlab='time step', ylab='exp', main=title)
  
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

time_of_divergence(ic25_pol64, ic25_pol113, 64)

time_of_divergence <- function(dn, dm, dn_nbits) {
  d0 <- sqrt(2^(-dn_nbits))
  d <- sqrt((dm$th1-dn$th1)^2 + (dm$w1-dn$w1)^2 + (dm$th2-dn$th2)^2 + (dm$w2-dn$w2)^2)
   
   for (i in 1:nrow(dn)) {
     if (d[i] > 10*d0 + sqrt(2^(-113))) {
         print(dn[i,]) 
         print(dm[i,])
         return()
     }
   }
}
  