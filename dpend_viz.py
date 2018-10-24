import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from matplotlib.patches import Circle

#https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

#path = "/MPFRoutput"

NSTEPS=0

mpfr8_file = open("./mpfr_data/mpfr8.txt", "r")
mpfr16_file = open("./mpfr_data/mpfr16.txt", "r")
mpfr32_file = open("./mpfr_data/mpfr32.txt", "r")
mpfr64_file = open("./mpfr_data/mpfr64.txt", "r")

t_8, t_16, t_32, t_64=([] for i in range(4))
x1_8, x1_16, x1_32, x1_64=([] for i in range(4))
y1_8, y1_16, y1_32, y1_64=([] for i in range(4))
x2_8, x2_16, x2_32, x2_64=([] for i in range(4))
y2_8, y2_16, y2_32, y2_64=([] for i in range(4))

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

for line in mpfr64_file.readlines():
  step_data_64 = line.split(' ', 5)
  t_64.append(float(step_data_64[0]))
  x1_64.append(float(step_data_64[1]))
  y1_64.append(float(step_data_64[2]))
  x2_64.append(float(step_data_64[3]))
  y2_64.append(float(step_data_64[4]))

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'b-', lw=2)
line2, = ax.plot([], [], 'g-', lw=2)
time_template = 't = %.6fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx_8 = [0, x1_8[i], x2_8[i]]
    thisy_8 = [0, y1_8[i], y2_8[i]]
    thisx_64 = [0, x1_64[i], x2_64[i]]
    thisy_64 = [0, y1_64[i], y2_64[i]]

    line.set_data(thisx_8, thisy_8)
    line2.set_data(thisx_64, thisy_64)
    time_text.set_text(time_template % t_8[i])
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, NSTEPS),
                              interval=25, blit=True, init_func=init)

# ani.save('double_pendulum.mp4', fps=15)
plt.show()
