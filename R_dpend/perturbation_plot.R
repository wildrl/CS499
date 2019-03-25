# type ="l" for lines only
scatter3D(x, y, z, phi = 0, bty = "g", type = "l", 
          ticktype = "detailed", lwd = 4)

arrows3D(x0, y0, z0, x1, y1, z1, colvar = x1^2, col = cols,
         lwd = 2, d = 3, clab = c("Quality", "score"), 
         main = "Arrows 3D", bty ="g", ticktype = "detailed")
# Add starting point of arrow
points3D(x0, y0, z0, add = TRUE, col="darkred", 
         colkey = FALSE, pch = 19, cex = 1)
# Add labels to the arrows
text3D(x1, y1, z1, c("Sepal.L", "Sepal.W", "Petal.L", "Petal.W"),
       colvar = x1^2, col = cols, add=TRUE, colkey = FALSE)



liapunov_plot <- function(f_n, fp_th1, fp_w1, fp_th2, fp_w2) {
  #x=th1, y=th2, z=w2
  
  
  scatter3D(f_n$th1, f_n$th2, f_n$w2, phi = 0, main = "Effects of Perterbations", 
            xlab = "th1", ylab ="th2", zlab = "w2",
            bty = "g", type = "l", 
            ticktype = "detailed", lwd = 2)
  
  for(i in seq(1, nrow(f_n), 1000)) {
    x0 <- f_n$th1[i]
    y0 <- f_n$th2[i]
    z0 <- f_n$w2[i]
    
    fp_th1_x1 <- fp_th1$th1[i]
    fp_th1_y1 <- fp_th1$th2[i]
    fp_th1_z1 <- fp_th1$w2[i]
    
    arrows3D(x0, y0, z0, fp_th1_x1, fp_th1_y1, fp_th1_z1, add=TRUE, col = "black",
             lwd = 2, d = 3, 
             main = "Arrows 3D", bty ="g", ticktype = "detailed")
    
    
   # fp_th2_x1 <- fp_th2$th1[i]
   # fp_th2_y1 <- fp_th2$th2[i]
   # fp_th2_z1 <- fp_th2$w2[i]
    
    
    
    
   # fp_th3_x1 <- fp_th3$th1[i]
   # fp_th3_y1 <- fp_th3$th2[i]
   # fp_th3_z1 <- fp_th3$w2[i]
    
    
  }
  
  
}