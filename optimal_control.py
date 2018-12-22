import numpy as np
from pylab import *
import itertools
import timeit
import sys
import sympy
from sympy import *
from matplotlib import *
from step5 import *


acc_opt = []
spe_opt = []
pos_opt = []
s = K_values[init]
t0 = int(init)
te = int(t_opt)

x0 = initx
v0 = initv
for t in range (t0,te+2):
    x = ( ((t)**3*(t0*v0 - te*v0 + t0*ve - te*ve - 2.0*x0 + 2.0*xe))/(t0 - te)**3 -((t)**2*(t0**2*v0 + t0*te*v0 - 2.0*te**2*v0 + 2.0*t0**2*ve - t0*te*ve - te**2*ve - 3.0*t0*x0 - 3.0*te*x0 + 3.0*t0*xe + 3.0*te*xe))/(t0 - te)**3 -((t)*(-2.0*t0**2*te*v0 + t0*te**2*v0 + te**3*v0 - t0**3*ve - t0**2*te*ve + 2.0*t0*te**2*ve + 6.0*t0*te*x0 - 6.0*t0*te*xe))/(t0 - te)**3 -(t0**2*te**2*v0 - t0*te**3*v0 + t0**3*te*ve - t0**2*te**2*ve - 3.0*t0*te**2*x0 + te**3*x0 - t0**3*xe + 3.0*t0**2*te*xe)/(t0 - te)**3)
    pos_opt.append(x)
for t in range (t0+1,te+2):
    u = ( (6.0*(t)*(t0*v0 - te*v0 + t0*ve - te*ve - 2.0*x0 + 2*xe))/(t0 - te)**3 -(2.0*(t0**2*v0 + t0*te*v0 - 2.0*te**2*v0 + 2.0*t0**2*ve - t0*te*ve - te**2*ve - 3.0*t0*x0 - 3.0*te*x0 + 3.0*t0*xe + 3.0*te*xe))/(t0 - te)**3 )
    v = ((3*t**2*(t0*v0 - te*v0 + t0*ve - te*ve - 2*x0 + 2*xe))/(t0 - te)**3 -(2*t*(t0**2*v0 + t0*te*v0 - 2*te**2*v0 + 2*t0**2*ve - t0*te*ve - te**2*ve - 3*t0*x0 - 3*te*x0 + 3*t0*xe + 3*te*xe))/(t0 - te)**3 -(-2*t0**2*te*v0 + t0*te**2*v0 + te**3*v0 - t0**3*ve - t0**2*te*ve + 2*t0*te**2*ve + 6*t0*te*x0 - 6*t0*te*xe)/(t0 - te)**3)
    #x = ( ((t)**3*(t0*v0 - te*v0 + t0*ve - te*ve - 2.0*x0 + 2.0*xe))/(t0 - te)**3 -((t)**2*(t0**2*v0 + t0*te*v0 - 2.0*te**2*v0 + 2.0*t0**2*ve - t0*te*ve - te**2*ve - 3.0*t0*x0 - 3.0*te*x0 + 3.0*t0*xe + 3.0*te*xe))/(t0 - te)**3 -((t)*(-2.0*t0**2*te*v0 + t0*te**2*v0 + te**3*v0 - t0**3*ve - t0**2*te*ve + 2.0*t0*te**2*ve + 6.0*t0*te*x0 - 6.0*t0*te*xe))/(t0 - te)**3 -(t0**2*te**2*v0 - t0*te**3*v0 + t0**3*te*ve - t0**2*te**2*ve - 3.0*t0*te**2*x0 + te**3*x0 - t0**3*xe + 3.0*t0**2*te*xe)/(t0 - te)**3)
    #pos_opt.append(x)
    acc_opt.append(u)
    spe_opt.append(v)

total_acc = u_opt + acc_opt
total_spe= v_opt+ spe_opt
total_pos = x_opt + pos_opt

print total_pos
host = host_subplot(111, axes_class=AA.Axes)
plt.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()

offset = 60
new_fixed_axis = par2.get_grid_helper().new_fixed_axis
par2.axis["right"] = new_fixed_axis(loc="right", axes=par2,
                                        offset=(offset, 0))

par2.axis["right"].toggle(all=True)

host.set_xlim(0, te+1)
host.set_ylim(0, xe)

host.set_xlabel("Stage")
host.set_ylabel("Position'")
par2.set_ylabel("Speed")
par1.set_ylabel("Acceleration")



p1, = host.plot(total_pos, label="Position",color='black')
p2, = par2.plot(total_spe, label="Speed")
p3, = par1.plot(total_acc, drawstyle='steps',label="Acceleration")
p4 = plt.axvline(x=T_max,color='red', linestyle=':')

par2.set_ylim(0,20,1.0)
par1.set_ylim(-5,5,1.0)

host.legend()

host.axis["left"].label.set_color(p1.get_color())
par2.axis["right"].label.set_color(p2.get_color())
par1.axis["right"].label.set_color(p3.get_color())

grid(True)
plt.draw()
plt.show()
