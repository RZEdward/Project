import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from IPython.display import display
from scipy.optimize import curve_fit

kB = 1.38*10**-23

r = np.arange(0.8,3,0.01)
model_F = []
V = []
ys = [-3,-2.5,-2,-1.5,-1,-0.5,0,0.5]
eps = 1
sig = 1.7818

R = 1
D = 2*R

yhigh = 10
ylow = -10
ymid = (yhigh+ylow)/2
yupperq = (yhigh+ymid)/2
txtpoint = (9*ylow + ymid)/10
txtshift = 0.95*txtpoint

max_F = 0
max_F_r = 0.8
for i in range(len(r)):
    element = -24*eps*(   ((2*sig**12)/(r[i]**13) - ((sig**6)/(r[i]**7)))  )
    model_F.append(element)
    if element > max_F:
        max_F = element
        max_F_r = r[i]
for i in range(len(r)):
    element = 4*eps*(   (sig/(r[i]))**12 - (sig/(r[i]))**6      )
    V.append(element)

label_max_F = 'Max Force = ' + str(round(max_F,2))

plt.figure()
plt.gca().set_ylim(-10,10)
plt.gca().set_xlim(0,5)
plt.grid()
plt.plot(r, model_F, "k--", label = 'LJ Force')
plt.plot(r, V, "r--", label = 'LJ Potential')
plt.scatter(D, txtpoint, color = 'b')
plt.text(D, txtshift, '$r = R_{i} + R_{j}$')
plt.scatter(max_F_r, max_F, color = 'b')
plt.text(max_F_r, 1.1*max_F, label_max_F)
plt.legend(loc = 'upper left')
plt.xlabel('$r$')
plt.ylabel('Force / Potential')
plt.show()