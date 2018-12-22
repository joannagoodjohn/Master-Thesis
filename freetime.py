import numpy as np
from pylab import *
import itertools
import timeit
import sys
import sympy
from sympy import *
from matplotlib import *
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
import random
import sys
import time
import pickle
start_time = time.time()

te= Symbol('te',real=True, positive=True)


w = 0.1

T_max = 15.0 ; T_min = 5.0 ; T_tl = 0.0 ;

ve = 11.0 ; xe = 350.0 ; X0 = 0 ; V0 = 0.0

dt = 1.0 ;
dx = 0.25;
dv = 0.5;
du = dv

if (len(sys.argv) > 1):
    dx = sys.argv[1]
if (len(sys.argv) > 2):
    dv = sys.argv[2]
if (len(sys.argv) > 3):
    du = sys.argv[3]


x_min = 0.0 ; x_max  = 300.0
v_min = 0.0 ; v_max = 11.0
u_min = -3.0; u_max = 3.0



x_values = np.arange(x_min,x_max+1,float(dx))
v_values = np.arange(v_min,v_max+1,float(dv))
u_values = np.arange(u_min,u_max+1,float(du))
K_values = np.arange(T_tl,T_max+1, dt)


X0j = np.abs(x_values-X0).argmin()
if (x_values[X0j] != X0):
    x_values = np.insert(x_values, X0j+1, X0)

V0j = np.abs(v_values-V0).argmin()
if (v_values[V0j] != V0):
    v_values = np.insert(v_values, V0j+1, V0)

#KX = int(((x_max-x_min)/float(dx))+1)
#KV = int(((v_max-v_min)/float(dv))+1)
KU = int(((u_max-u_min)/float(du))+1)
KMAX = int(((T_max-T_tl)/dt)+1)
KV = len(v_values)
KX = len(x_values)

tf = np.zeros([KMAX+1, KX, KV])

t_min_temp = np.where( K_values == T_min )
t_min = t_min_temp[0]

t_max_temp = np.where( K_values == T_max )
t_max = t_max_temp[0]
c = (t_min+t_max)/2.0


def tfinal(x,v,ve,xe,t0):
    x0 = x
    v0 = v
    eqn =-4*t0**4*v0**2 - 10*t0**4*v0*ve - 10*t0**4*ve**2 + t0**6*w - 6*t0*te**5*w + te**6*w + te**4*(-10*v0**2 - 10*v0*ve - 4*ve**2 + 15*t0**2*w) + 24*t0**3*v0*x0 + 36*t0**3*ve*x0 - 36*t0**2*x0**2 - 24*t0**3*v0*xe - 36*t0**3*ve*xe + 72*t0**2*x0*xe - 36*t0**2*xe**2 + te**3*(16*t0*v0**2 + 4*t0*v0*ve + 4*t0*ve**2 - 20*t0**3*w - 36*v0*x0 - 24*ve*x0 + 36*v0*xe + 24*ve*xe) + te*(4*t0**3*v0**2 + 4*t0**3*v0*ve + 16*t0**3*ve**2 - 6*t0**5*w - 12*t0**2*v0*x0 - 24*t0**2*ve*x0 + 12*t0**2*v0*xe + 24*t0**2*ve*xe) + te**2*(-6*t0**2*v0**2 + 12*t0**2*v0*ve - 6*t0**2*ve**2 + 15*t0**4*w + 24*t0*v0*x0 + 12*t0*ve*x0 - 36*x0**2 - 24*t0*v0*xe - 12*t0*ve*xe + 72*x0*xe - 36*xe**2)
    sol = solveset(eqn, te, domain = S.Reals)
    result = [s for s in sol if s > 0] or None
    ts = min(result)
    return  ts




for k in reversed(range(KMAX)):
    for i in range (KX):
        for j in range (KV):
            x = x_values[i]
            v = v_values[j]
            t0 = K_values[k]
            if t_min<=k<=t_max:
                if k == KMAX-1:
                    res = 0
                    tf[k,i,j] = tfinal(x,v,ve,xe,t0)
                    print x,v,tf[k,i,j]

                else:

                    tf[k,i,j] = tf[k+1,i,j]+1.0
output = open('data.pkl', 'wb')
pickle.dump(tf, output)
output.close()







expcost = 0
for l in range(t_min,KMAX):
    usofar = u_opt[:l]
    t0 = l
    x = X0
    v = V0
    phisum = 0
    for u in usofar:
        (x, v) = xvnext(x, v, u)
        phisum = phisum + phi(u)
        j_opt_here = 0
        (te) = tfinal(x,v,ve,xe,t0)
    for t in range(int(l), int(te)):
        j_opt_here = j_opt_here + ( (6.0*t*(t0*v - te*v + t0*ve - te*ve - 2.0*x + 2.0*xe))/(t0 - te)**3.0 -(2.0*(t0**2.0 * v + t0*te*v - 2.0 * te**2.0 * v + 2.0 * t0**2.0 * ve - t0*te*ve - te**2.0 * ve - 3.0*t0*x - 3.0*te*x + 3.0*t0*xe +3.0*te*xe))/(t0 - te)**3.0 )**2

    p = 1.0/((T_max-T_min)+1.0)
    if (chkpoint >= t_min):
        if (l == chkpoint):
            p = 1.0
        else:
            p = 0.0

    expcost = expcost + (p)*( phisum +j_opt_here)

    #print "l %s u %s cost %s jopt %s expcost %s x %s v %s" % (l, usofar, phisum, j_opt_here, expcost,x ,v)
    (x,v) = xvnext(x,v,u)


print x_opt
if (chkpoint >= t_min):
    print "controls to checkpoint %s: %s" % (chkpoint, u_opt[:chkpoint])
print "dx %s dv %s du %s time %s cost %s" % (dx, dv, du,time, expcost)
