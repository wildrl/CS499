#IC2
mag_plot(ic2_mag11$magnitude, ic2_mag24$magnitude, ic2_mag53$magnitude, ic2_mag64$magnitude, ic2_mag113$magnitude, 
        "ic2 magnitudes", nrow(ic2_mag11), max(ic2_mag113$magnitude))
dot_plot(ic2_mag11$dot_product, ic2_mag24$dot_product, ic2_mag53$dot_product, ic2_mag64$dot_product, 
       "ic2 dot product with 113-bit soln", nrow(ic2_mag11), max(ic2_mag64$dot_product))

#IC3
mag_plot(ic3_mag11$magnitude, ic3_mag24$magnitude, ic3_mag53$magnitude, ic3_mag64$magnitude, ic3_mag113$magnitude, 
         "ic3=(180,0,90,0): magnitude", nrow(ic3_mag11), max(ic3_mag64$magnitude))
dot_plot(ic3_mag11$dot_product, ic3_mag24$dot_product, ic3_mag53$dot_product, ic3_mag64$dot_product, 
         "ic3=(180,0,90,0): dot product with 113-bit soln", nrow(ic3_mag11), max(ic3_mag64$dot_product))

#IC6
mag_plot(ic6_mag11$magnitude, ic6_mag24$magnitude, ic6_mag53$magnitude, ic6_mag64$magnitude, ic6_mag113$magnitude, 
         "ic6=(180,0,91,0): magnitude", nrow(ic6_mag11), max(ic6_mag24$magnitude))
dot_plot(ic6_mag11$dot_product, ic6_mag24$dot_product, ic6_mag53$dot_product, ic6_mag64$dot_product, 
         "ic6=(180,0,91,0): dot product with 113-bit soln", nrow(ic6_mag11), max(ic6_mag64$dot_product))



mag_plot <- function(d11,d24,d53,d64,d113, title, nrows, ymax) {
  plot(1,type='n',xlim=c(0,nrows),ylim=c(0,ymax),xlab='time step', ylab='magnitude', main=title)

  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  lines(d113, type='l', col="purple", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "green", "blue", "purple"), lty=1:1, cex=0.8)
}

dot_plot <- function(d11,d24,d53,d64, title, nrows, ymax, xmax) {
  plot(1,type='n',xlim=c(0,nrows),ylim=c(0,ymax),xlab='time step', ylab='dot product', main=title)
  
  
  lines(d11, type='l', col="red", lwd=1.5)
  lines(d24, type='l', col="orange", lwd=1.5)
  lines(d53, type='l', col="green", lwd=1.5)
  lines(d64, type='l', col="blue", lwd=1.5)
  
  legend("topleft", legend=c("11", "24", "53", "64"),
         col=c("red", "orange", "green", "blue"), lty=1:1, cex=0.8)
}

mag_plot113(ic6_mag113$magnitude, ic3_mag113$magnitude,
            "ic6 & ic3 113-bit soln magnitude", nrow(ic6_mag113),100)

mag_plot113 <- function(d1,d2, title, nrows, ymax) {
  plot(1,type='n',xlim=c(0,nrows),ylim=c(0,ymax),xlab='time step', ylab='magnitude', main=title)
  
  lines(d1, type='l', col="red", lwd=1.5)
  lines(d2, type='l', col="orange", lwd=1.5)
  
  legend("topleft", legend=c("113 ic6", "113 ic3"),
         col=c("red", "orange"), lty=1:1, cex=0.8)
}


dot_plot113(ic6_pol113, ic3_pol113, "dot product of ic3 and ic6 113 bit solns", 
            nrow(ic6_pol113), 50)

dot_plot113 <- function(d1,d2, title, nrows, ymax) {
  plot(1,type='n',xlim=c(0,nrows),ylim=c(0,ymax),xlab='time step', ylab='dot product', main=title)
  
  dp <- sqrt(d1$th1 * d2$th1 + d1$w1 * d2$w1 + d1$th2 * d2$th2 + d1$w2 * d2$w2)

  lines(dp, type='l', col="red", lwd=1.5)
}






