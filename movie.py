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

rowsR1 = []
with open('R1.csv') as file1:
    reader_obj = csv.reader(file1)
    for row in reader_obj:
        rowsR1.append(row)

rowsR2 = []
with open('R1.csv') as file2:
    reader_obj = csv.reader(file2)
    for row in reader_obj:
        rowsR1.append(row)


for i in range(timesteps):

    x_posR1 = []
    y_posR1 = []
    x_posR2 = []
    y_posR2 = []
    for j in range(len(rowsR1[0])):
        if j % 2 == 0:
            x_posR1.append(rowsR1[i+1][j])
        else:
            y_posR1.append(rowsR1[i+1][j])
    for j in range(len(rowsR2[0])):
        if j % 2 == 0:
            x_posR1.append(rowsR2[i+1][j])
        else:
            y_posR1.append(rowsR2[i+1][j])

    




