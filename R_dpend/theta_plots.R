
th_plot(ic10_pol11$th1, ic10_pol24$th1, ic10_pol53$th1, ic10_pol64$th1, ic10_pol113$th1, "ic10 th1 plot")
th_plot(ic10_pol11$th1, ic10_pol24$th1, ic10_pol53$th1, ic10_pol64$th1, ic10_pol113$th1, "ic10 th2 plot")

# This function plots the value of theta at each time step
# for all parameters on a single graph.
th_plot <- function(d11,d24,d53,d64,d113, title) {
  plot(1,type='n',xlim=c(1,100000),ylim=c(-pi,pi),xlab='time step', ylab='theta (radians)',
       main=title)
  
  n1 <- d11 %/% (2*pi)
  y1 <- (d11) - (n1 * (2 * pi))
  y1 <- ifelse(y1>pi,y1-2*pi,y1)
  lines(y1, type='l', col="red", lwd=1.5)
  
  n2 <- d24 %/% (2*pi)
  y2 <- (d24) - (n2 * (2 * pi))
  y2 <- ifelse(y2>pi,y2-2*pi,y2)
  lines(y2, type='l', col="orange", lwd=1.5)
  
  n3 <- d53 %/% (2*pi)
  y3 <- (d53) - (n3 * (2 * pi))
  y3 <- ifelse(y3>pi,y3-2*pi,y3)
  lines(y3, type='l', col="yellow", lwd=1.5)
  
  n4 <- d64 %/% (2*pi)
  y4 <- (d64) - (n4 * (2 * pi))
  y4 <- ifelse(y4>pi,y4-2*pi,y4)
  lines(y4, type='l', col="blue", lwd=1.5)
  
  n5 <- d113 %/% (2*pi)
  y5 <- (d113) - (n5 * (2 * pi))
  y5 <- ifelse(y5>pi,y5-2*pi,y5)
  lines(y5, type='l', col="purple", lwd=1.5)
  
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  legend("topright", inset=c(-0.4,0),legend=c("11", "24", "53", "64", "113"),
         col=c("red", "orange", "yellow", "blue", "purple"), lty=1:1, cex=0.8)
}
