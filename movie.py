import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from IPython.display import display

properties = pd.read_csv('properties.csv')
R1_particles = pd.read_csv('R1.csv')
R2_particles = pd.read_csv('R2.csv')

dpi = properties.dpi[0]
w = properties.w[0]
h = properties.h[0]
timesteps = properties.timesteps[0]
dt = properties.dt[0]
boxlims = properties.boxlims[0]
markersize1 = properties.markersize1[0]
markersize2 = properties.markersize2[0]
numR1 = properties.numR1[0]
numR2 = properties.numR2[0]

rowsR1 = []
with open('R1.csv') as file1:
    reader_obj = csv.reader(file1)
    for row in reader_obj:
        rowsR1.append(row)

rowsR2 = []
with open('R2.csv') as file2:
    reader_object = csv.reader(file2)
    for row in reader_object:
        rowsR2.append(row)

for i in range(1,timesteps+1):

    x_posR1 = []
    y_posR1 = []
    x_posR2 = []
    y_posR2 = []
    for j in range(numR1):
        if j % 2 == 0:
            x_posR1.append(rowsR1[i][j])
        else:
            y_posR1.append(rowsR1[i][j])
    for j in range(numR2):
        if j % 2 == 0:
            x_posR2.append(rowsR2[i][j])
        else:
            y_posR2.append(rowsR2[i][j])

    x = x_posR1
    y = y_posR1
    print(x)
    print(y)
    plt.figure(1)
    plt.plot(x,y,"ko")
    #plt.plot(x_posR2,y_posR2,"ko")
    iteration = str(i)
    png = ".png"
    loc = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/slide"
    savepoint = loc + iteration + png
    plt.savefig(savepoint)


    




