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

for i in range(len(r)):
    element = -24*eps*(   ((2*sig**12)/(r[i]**13) - ((sig**6)/(r[i]**7)))  )
    model_F.append(element)
for i in range(len(r)):
    element = 4*eps*(   (sig/(r[i]))**12 - (sig/(r[i]))**6      )
    V.append(element)

plt.figure()
plt.gca().set_ylim(-10,10)
plt.grid()
plt.plot(r, model_F, "k--")
plt.plot(r, V, "r--")
plt.show()