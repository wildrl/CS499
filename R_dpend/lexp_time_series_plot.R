install.packages("plot3D")
library(plot3D)
library(datasets)

plot(1,type='n',xlim=c(1,200),ylim=c(-1.5,1.5),xlab='time step', ylab='exp')

lines(`ly_exp16`$exp, type='l', col="red", lwd=1.5)
lines(`ly_exp23`$exp, type='l', col="orange", lwd=1.5)
lines(`ly_exp32`$exp, type='l', col="yellow", lwd=1.5)
lines(`ly_exp48`$exp, type='l', col="blue", lwd=1.5)
lines(`ly_exp52`$exp, type='l', col="purple", lwd=1.5)
lines(`ly_exp64`$exp, type='l', col="black", lwd=1.5)