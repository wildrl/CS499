import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from matplotlib.patches import Circle

NSTEPS=0

mpfr8_file = open("./mpfr_data/cartesian8.txt", "r")
mpfr16_file = open("./mpfr_data/cartesian16.txt", "r")
mpfr32_file = open("./mpfr_data/cartesian32.txt", "r")
mpfr64_file = open("./mpfr_data/cartesian64.txt", "r")
mpfr128_file = open("./mpfr_data/cartesian128.txt", "r")

t_8, t_16, t_32, t_64, t_128 = ([] for i in range(5))
x1_8, x1_16, x1_32, x1_64, x1_128 = ([] for i in range(5))
y1_8, y1_16, y1_32, y1_64, y1_128 = ([] for i in range(5))
x2_8, x2_16, x2_32, x2_64, x2_128 = ([] for i in range(5))
y2_8, y2_16, y2_32, y2_64, y2_128 = ([] for i in range(5))

for line in mpfr8_file.readlines():
  step_data_8 = line.split(' ', 5)
  t_8.append(float(step_data_8[0]))
  x1_8.append(float(step_data_8[1]))
  y1_8.append(float(step_data_8[2]))
  x2_8.append(float(step_data_8[3]))
  y2_8.append(float(step_data_8[4]))
  NSTEPS+=1

for line in mpfr16_file.readlines():
  step_data_16 = line.split(' ', 5)
  t_16.append(float(step_data_16[0]))
  x1_16.append(float(step_data_16[1]))
  y1_16.append(float(step_data_16[2]))
  x2_16.append(float(step_data_16[3]))
  y2_16.append(float(step_data_16[4]))

for line in mpfr32_file.readlines():
  step_data_32 = line.split(' ', 5)
  t_32.append(float(step_data_32[0]))
  x1_32.append(float(step_data_32[1]))
  y1_32.append(float(step_data_32[2]))
  x2_32.append(float(step_data_32[3]))
  y2_32.append(float(step_data_32[4]))

for line in mpfr64_file.readlines():
  step_data_64 = line.split(' ', 5)
  t_64.append(float(step_data_64[0]))
  x1_64.append(float(step_data_64[1]))
  y1_64.append(float(step_data_64[2]))
  x2_64.append(float(step_data_64[3]))
  y2_64.append(float(step_data_64[4]))

for line in mpfr128_file.readlines():
  step_data_128 = line.split(' ', 5)
  t_128.append(float(step_data_128[0]))
  x1_128.append(float(step_data_128[1]))
  y1_128.append(float(step_data_128[2]))
  x2_128.append(float(step_data_128[3]))
  y2_128.append(float(step_data_128[4]))

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.grid()

line8, = ax.plot([], [], 'r-', lw=2)
line16, = ax.plot([], [], 'y-', lw=2)
line32, = ax.plot([], [], 'g-', lw=2)
line64, = ax.plot([], [], 'b-', lw=2)
line128, = ax.plot([], [], 'm-', lw=2)
time_template = 't = %.6fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

step_template = 'step = %fs'
step_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line8.set_data([], [])
    line16.set_data([], [])
    line32.set_data([], [])
    line64.set_data([], [])
    line128.set_data([], [])
    #time_text.set_text('')
    step_text.set_text('')
    return line8, line16, line32, line64, line128, step_text #time_text


def animate(i):
    thisx_8 = [0, x1_8[i], x2_8[i]]
    thisy_8 = [0, y1_8[i], y2_8[i]]
    thisx_16 = [0, x1_16[i], x2_16[i]]
    thisy_16 = [0, y1_16[i], y2_16[i]]
    thisx_32 = [0, x1_32[i], x2_32[i]]
    thisy_32 = [0, y1_32[i], y2_32[i]]
    thisx_64 = [0, x1_64[i], x2_64[i]]
    thisy_64 = [0, y1_64[i], y2_64[i]]
    thisx_128 = [0, x1_128[i], x2_128[i]]
    thisy_128 = [0, y1_128[i], y2_128[i]]

    line8.set_data(thisx_8, thisy_8)
    line16.set_data(thisx_16, thisy_16)
    line32.set_data(thisx_32, thisy_32)
    line64.set_data(thisx_64, thisy_64)
    line128.set_data(thisx_128, thisy_128)
    #time_text.set_text(time_template % t_8[i])
    step_text.set_text(step_template % i)
    return line8, line16, line32, line64, line128, step_text #time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, NSTEPS),
                              interval=25, blit=True, init_func=init)

plt.show()
