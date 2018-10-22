import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

#path = "/MPFRoutput"

NSTEPS=1000
step=0

mpfr8_file = open("./mpfr_data/mpfr8.txt", "r")
mpfr16_file = open("./mpfr_data/mpfr16.txt", "r")
mpfr32_file = open("./mpfr_data/mpfr32.txt", "r")
mpfr64_file = open("./mpfr_data/mpfr64.txt", "r")

t=[NSTEPS]
x1=[NSTEPS]
y1=[NSTEPS]
x2=[NSTEPS]
y2=[NSTEPS]

for line in mpfr8_file.readlines():
  step_data = line.split(' ', 5)
  t.append(float(step_data[0]))
  x1.append(float(step_data[1]))
  y1.append(float(step_data[2]))
  x2.append(float(step_data[3]))
  y2.append(float(step_data[4]))
  step+=1

plt.plot(4,7)



#make_movie = TRUE
#params = {'backend': 'ps',
#          'font.size': 20,
#          'font.family':'serif',
#          'font.serif':['Times'],
#          'ps.usedistiller': 'xpdf',
#          'text.usetex': True,
#          }
#plt.rcParams.update(params)
#plt.ion()
#fig = plt.figure(figsize=(9.6,5.4),dpi=100) # 1920x1080
#fig.subplots_adjust(left=0, right=1, top=1, bottom=0,hspace=0.02,wspace=0.02)
#ax = fig.add_subplot(111)
#ax.axis('off') # no frame

