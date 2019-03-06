
th_plot(ic24_pol11, ic24_pol24, ic24_pol53, ic24_pol64, ic24_pol113, "ic24")

th_plot <- function(d11,d24,d53,d64,d113, title) {
  plot(1,type='n',xlim=c(1,10000),ylim=c(-pi,2*pi+1),xlab='time step', ylab='theta (radians)',
       main=title)
  
  n1 <- d11$th1 %/% (2*pi)
  y1 <- (d11$th1) - (n1 * (2 * pi))
  y1 <- ifelse(y1>pi,y1-2*pi,y1)
  lines(y1, type='l', col="red", lwd=1.5)
  
  n2 <- d24$th1 %/% (2*pi)
  y2 <- (d24$th1) - (n2 * (2 * pi))
  y2 <- ifelse(y2>pi,y2-2*pi,y2)
  lines(y2, type='l', col="orange", lwd=1.5)
  
  n3 <- d53$th1 %/% (2*pi)
  y3 <- (d53$th1) - (n3 * (2 * pi))
  y3 <- ifelse(y3>pi,y3-2*pi,y3)
  lines(y3, type='l', col="yellow", lwd=1.5)
  
  n4 <- d64$th1 %/% (2*pi)
  y4 <- (d64$th1) - (n4 * (2 * pi))
  y4 <- ifelse(y4>pi,y4-2*pi,y4)
  lines(y4, type='l', col="blue", lwd=1.5)
  
  n5 <- d113$th1 %/% (2*pi)
  y5 <- (d113$th1) - (n5 * (2 * pi))
  y5 <- ifelse(y5>pi,y5-2*pi,y5)
  lines(y5, type='l', col="purple", lwd=1.5)
  
  legend("topright",legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}
