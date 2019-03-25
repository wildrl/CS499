import_mag_data()

import_mag_data <- function() {
  #IC2 normal
  ic2_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/mag11.csv")
  ic2_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/mag24.csv")
  ic2_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/mag53.csv")
  ic2_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/mag64.csv")
  ic2_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/mag113.csv")
  
  #IC3 preturb th1
  ic3_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/mag11.csv")
  ic3_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/mag24.csv")
  ic3_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/mag53.csv")
  ic3_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/mag64.csv")
  ic3_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/mag113.csv")
  
  #IC4 preturb w1
  ic4_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/mag11.csv")
  ic4_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/mag24.csv")
  ic4_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/mag53.csv")
  ic4_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/mag64.csv")
  ic4_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/mag113.csv")
  
  #IC5 preturb th2
  ic5_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/mag11.csv")
  ic5_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/mag24.csv")
  ic5_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/mag53.csv")
  ic5_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/mag64.csv")
  ic5_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/mag113.csv")
  
  #IC6 preturb w2
  ic6_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/mag11.csv")
  ic6_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/mag24.csv")
  ic6_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/mag53.csv")
  ic6_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/mag64.csv")
  ic6_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/mag113.csv")
  
  
  
  #IC7 no chaos normal
  ic7_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/mag11.csv")
  ic7_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/mag24.csv")
  ic7_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/mag53.csv")
  ic7_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/mag64.csv")
  ic7_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/mag113.csv")
  
  #IC7 no chaos preturb th1
  ic8_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/mag11.csv")
  ic8_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/mag24.csv")
  ic8_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/mag53.csv")
  ic8_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/mag64.csv")
  ic8_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/mag113.csv")
  
  #IC7 no chaos preturb w1
  ic9_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/mag11.csv")
  ic9_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/mag24.csv")
  ic9_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/mag53.csv")
  ic9_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/mag64.csv")
  ic9_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/mag113.csv")
  
  #IC7 no chaos preturb th2
  ic10_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/mag11.csv")
  ic10_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/mag24.csv")
  ic10_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/mag53.csv")
  ic10_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/mag64.csv")
  ic10_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/mag113.csv")
  
  #IC7 no chaos preturb w2
  ic11_mag11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/mag11.csv")
  ic11_mag24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/mag24.csv")
  ic11_mag53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/mag53.csv")
  ic11_mag64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/mag64.csv")
  ic11_mag113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/mag113.csv")
}
