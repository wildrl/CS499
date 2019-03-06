install.packages("plot3D")
library(plot3D)
library(datasets)

plot(1,type='n',xlim=c(1,200),ylim=c(-1.5,1.5),xlab='time step', ylab='exp')

lines(`ly_exp16`$lexp, type='l', col="red", lwd=1.5)
lines(`ly_exp23`$lexp, type='l', col="orange", lwd=1.5)
lines(`ly_exp32`$lexp, type='l', col="yellow", lwd=1.5)
lines(`ly_exp48`$lexp, type='l', col="blue", lwd=1.5)
lines(`ly_exp52`$lexp, type='l', col="purple", lwd=1.5)
lines(`ly_exp64`$lexp, type='l', col="black", lwd=1.5)