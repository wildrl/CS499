# plot parameters and magnitude
fx_plot <- function(d11,d24,d53,d64,m, title, xmax, ymax) {
  plot(1,type='n', xlim=c(-20,xmax), ylim=c(-40,ymax),
       xlab='time step', ylab='L2 norm', main=title)
  
  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  lines(m, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("th1", "w1", "th2", "w2", "L2"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

#Plot relative error
r_error_plot <- function(d11, d24, d53, d64, title, xmax, ymax) {
  plot(1,type='n', xlim=c(0,xmax), ylim=c(0,ymax), 
       xlab='time step', ylab='relative error', main=title)
  
  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64"),
         col=c("red", "orange", "green", "blue"), lty=1:1, cex=0.8)
}

# Plot L2 Norms
mag_plot <- function(d11, d24, d53, d64, d113, title, xmax, ymax) {
  plot(1,type='n', xlim=c(0,xmax), ylim=c(0,ymax),
       xlab='time step', ylab='L2 norm', main=title)

  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  lines(d113, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

# Plot dot product
dot_plot <- function(d11, d24, d53, d64, title, xmax, ymax, xmax) {
  plot(1,type='n', xlim=c(0,xmax), ylim=c(0,ymax),
       xlab='time step', ylab='dot product', main=title)
  
  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64"),
         col=c("red", "orange", "green", "blue"), lty=1:1, cex=0.8)
}



# misc

mag_plot113 <- function(fn, fp_th1, fp_w1, fp_th2, fp_w2, title, nrows, ymax, a, b, c, d, e) {
  
  format(fn,scientific=FALSE)
  plot(1,type='n', xlim=c(0,nrows), ylim=c(0,ymax), 
       xlab='time step', ylab='L2-norm', main=title)
  
  lines(fn, type='l', col="red", lwd=1.5)
  lines(fp_th1, type='l', col="orange", lwd=1.5)
  lines(fp_w1, type='l', col="green", lwd=1.5)
  lines(fp_th2, type='l', col="blue", lwd=1.5)
  lines(fp_w2, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c(a, b, c, d, e),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

r_error <- function(fn, fp_th1, fp_w1, fp_th2, fp_w2, title, xmax, ymax, a, b, c, d) {
  
  format(fn,scientific=FALSE)
  plot(1,type='n', xlim=c(0,xmax), ylim=c(0,ymax), 
       xlab='time step', ylab='relative error', main=title)
  
  lines(fn, type='l', col="red", lwd=1.5)
  lines(fp_th1, type='l', col="orange", lwd=1.5)
  lines(fp_w1, type='l', col="green", lwd=1.5)
  lines(fp_th2, type='l', col="blue", lwd=1.5)
  lines(fp_w2, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c(a, b, c, d,e),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}


dot_plot113(ic6_pol113, ic3_pol113, "dot product of ic3 and ic6 113 bit solns", 
            nrow(ic6_pol113), 50)

dot_plot113 <- function(d1,d2, title, nrows, ymax) {
  plot(1,type='n',xlim=c(0,nrows),ylim=c(0,ymax),xlab='time step', ylab='dot product', main=title)
  
  dp <- sqrt(d1$th1 * d2$th1 + d1$w1 * d2$w1 + d1$th2 * d2$th2 + d1$w2 * d2$w2)

  lines(dp, type='l', col="red", lwd=1.5)
}






