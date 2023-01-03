import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from IPython.display import display

w = 10
h = 10
boxlims = 10

R = 1
r = 2

mrk1 = 2*R*h*500/(9*boxlims)
mrk2 = 2*r*h*500/(9*boxlims)

x = [1,3,5,7,9]
y = [7,7,7,7,7]

a = [3,7]
b = [3,3]

plt.figure(figsize=(w,h), dpi = 80)
plt.gca().set_xlim(0,boxlims)
plt.gca().set_ylim(0,boxlims)
plt.xticks([0,boxlims])
plt.yticks([0,boxlims])
plt.plot(x, y, marker = 'o', color = 'white', markeredgecolor = 'black', linestyle = 'none', markersize = mrk1)
plt.plot(a, b, marker = 'o', color = 'white', markeredgecolor = 'red', linestyle = 'none', markersize = mrk2)

plt.savefig("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/test.png")