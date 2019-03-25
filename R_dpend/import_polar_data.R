import_polar_data()

import_polar_data <- function() {
  #IC2 normal
  ic2_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/polar11.csv")
  ic2_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/polar24.csv")
  ic2_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/polar53.csv")
  ic2_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/polar64.csv")
  ic2_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_2/polar113.csv")
  
  #IC3 preturb th1
  ic3_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/polar11.csv")
  ic3_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/polar24.csv")
  ic3_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/polar53.csv")
  ic3_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/polar64.csv")
  ic3_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_3/polar113.csv")
  
  #IC4 preturb w1
  ic4_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/polar11.csv")
  ic4_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/polar24.csv")
  ic4_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/polar53.csv")
  ic4_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/polar64.csv")
  ic4_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_4/polar113.csv")
  
  #IC5 preturb th2
  ic5_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/polar11.csv")
  ic5_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/polar24.csv")
  ic5_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/polar53.csv")
  ic5_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/polar64.csv")
  ic5_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_5/polar113.csv")
  
  #IC6 preturb w2
  ic6_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/polar11.csv")
  ic6_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/polar24.csv")
  ic6_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/polar53.csv")
  ic6_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/polar64.csv")
  ic6_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_6/polar113.csv")
  
  
  #IC7 non chaotic normal
  ic7_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/polar11.csv")
  ic7_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/polar24.csv")
  ic7_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/polar53.csv")
  ic7_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/polar64.csv")
  ic7_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_7/polar113.csv")
  
  #IC8 non chaotic preturb th1
  ic8_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/polar11.csv")
  ic8_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/polar24.csv")
  ic8_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/polar53.csv")
  ic8_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/polar64.csv")
  ic8_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_8/polar113.csv")
  
  #IC9 non chaotic preturb w1
  ic9_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/polar11.csv")
  ic9_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/polar24.csv")
  ic9_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/polar53.csv")
  ic9_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/polar64.csv")
  ic9_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_9/polar113.csv")
  
  #IC7 non chaotic preturb th2
  ic10_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/polar11.csv")
  ic10_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/polar24.csv")
  ic10_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/polar53.csv")
  ic10_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/polar64.csv")
  ic10_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_10/polar113.csv")
  
  #IC7 non chaotic preturb w2
  ic11_pol11 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/polar11.csv")
  ic11_pol24 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/polar24.csv")
  ic11_pol53 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/polar53.csv")
  ic11_pol64 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/polar64.csv")
  ic11_pol113 <- read.csv(file="~/Desktop/mpfr_data/mpfr_data/ic_11/polar113.csv")
}
