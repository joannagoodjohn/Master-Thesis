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
import pprint, pickle
import pprint, pickle

pkl_file = open('data.pkl', 'rb') #txt file with final times stored

data1 = pickle.load(pkl_file)
#pprint.pprint(data1)

w = 0.1

T_max = 15.0 ; T_min = 5.0 ; T_tl = 0.0 ;

ve = 11.0 ; xe = 120.0 ; X0 = 100 ; V0 = 9.25
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
v_min = 0.0 ; v_max = 10.0
u_min = -3.0; u_max = 3.0

chkpoint = 15
use_chkpoint = True

if (use_chkpoint == True):
    T_max = chkpoint

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

J_opt = np.zeros([KMAX+2,KX,KV])
J_star = np.full([KMAX+2,KX,KV],np.inf)
U = np.zeros([KMAX+1,KX,KV])
s = np.zeros([KMAX+1])
p = np.zeros([KMAX+1])



t_min_temp = np.where( K_values == T_min )
t_min = t_min_temp[0]

t_max_temp = np.where( K_values == T_max )
t_max = t_max_temp[0]
c= 15.0


def xvnext(x0,v0,u):
    global dt

    v1 = v0+u*dt
    x1 = x0+v0*dt + 0.5*u*(dt**2)
    return (x1,v1)


def ijnext(x0, v0, u):
    global x_values, v_values, x_max, v_max, dt

    (x1, v1) = xvnext(x0, v0, u)

    if (x1 < 0 or x1 > x_max or v1 < 0 or v1 > v_max ):
        return (-1, -1)

    inext = np.abs(x_values-x1).argmin()
    jnext = np.abs(v_values-v1).argmin()

    return (inext, jnext)

def phi(u):
    return 0.5*u**2

start_time = time.time()
for k in reversed(range(KMAX)):
    print k
    t0 = K_values[k]
    for i in range (KX):
        x0 = x_values[i]
        for j in range (KV):
            v0 = v_values[j]
            res = 0
            te  = data1[k,i,j]

            for t in range(int(t0), int(te)):
                res = res + ( (6.0*t*(t0*v0 - te*v0 + t0*ve - te*ve - 2.0*x0 + 2.0*xe))/(t0 - te)**3.0 -(2.0*(t0**2.0 * v0 + t0*te*v0 - 2.0 * te**2.0 * v0 + 2.0 * t0**2.0 * ve - t0*te*ve - te**2.0 * ve - 3.0*t0*x0 - 3.0*te*x0 + 3.0*t0*xe +3.0*te*xe))/(t0 - te)**3.0 )**2
                J_opt [k,i,j] = res

            if t_min<=k<=t_max:

                #if t_min<=k<=c:
                    #p[k]=(2*(k-t_min))/((t_max-t_min)*(c-t_min))

                #elif k==c:
                    #p[k]=2/(t_max-t_min)

                #elif c<=k<=t_max:
                    #p[k]=(2*(t_max-k))/((t_max-t_min)*(t_max-c))
                p[k] = 1.0/((t_max-k)+1.0)

                if (use_chkpoint):

                    if (k == chkpoint):
                        p[k] = 1.0
                    else:
                        p[k] = 0
            else:
                p[k] = 0

            if (k == KMAX-1):
                J_star[k,i,j] = J_opt[k,i,j]
                continue

            #J_star[k,i,j] =float("inf")
            for u in u_values:

                (inext, jnext) = ijnext(x0, v0, u)

                if (inext < 0 or jnext < 0 ):
                    continue

                cost = (1-p[k])*(phi(u)+J_star[k+1,inext,jnext])+(p[k])*J_opt[k,i,j]

                if cost < J_star[k,i,j]:
                    J_star[k,i,j] = cost
                    control = u

            U[k,i,j] = control







(time.time() - start_time)
cost = 0
x_opt = []
v_opt = []
u_opt = []
Cost = 0
real = 0
i_max_temp = np.where( x_values == x_max )
i_max = i_max_temp[0]
x0 = X0
v0 = V0
for k in range (KMAX):
    i_temp = np.where( x_values == x0 )
    i = i_temp[0]
    j_temp = np.where( v_values == v0 )
    j = j_temp[0]
    x_opt.append(x0)
    v_opt.append(v0)
    u_opt.append(U[k,i,j][0])
    #trajectory cost
    real +=0.5*U[k,i,j]**2
    if k < t_min:
        Cost +=0.5*U[k,i,j]**2
    else:
        Cost +=s[k]*J_opt[k,i,j]+(1-s[k])*(0.5*U[k,i,j]**2)
    if k==chkpoint:
        real+=J_opt[k,i,j]
    print('stage:',k,'position:',x0,'speed:',v0,'control:',U[k,i,j][0], 'cost: ',J_star[k,i,j][0], 'jopt:' ,J_opt[k,i,j][0], 'te:', data1[k,i,j],'trajectory cost:',real[0], 'Hybrid:', Cost)
    v = (v0 + U[k,i,j]*dt)
    x = (x0 + v0*dt + 0.5*U[k,i,j]*(dt**2))

    if (x < 0 or x > x_max or v < 0 or v > v_max):
        i = i
        j = j
        v0 = v_values[j]
        x0 = x_values[i]

    else:
        j =np.abs(v_values-v).argmin()
        v0 = v_values[j]
        i =np.abs(x_values-x).argmin()
        x0 = x_values[i]
    if k == chkpoint:
        initx = x_values[i]
        initv = v_values[j]
        init = k
        t_opt = data1[k,i,j]
        break

print x_opt
s = (KMAX*KU*KX*KV)*10**(-6)
print("--- %s seconds ---" % (time.time() - start_time))
print("Complexity:", s)
