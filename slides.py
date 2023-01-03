import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from IPython.display import display

properties = pd.read_csv('properties.csv')

info = pd.read_csv('info.csv')


timesteps = properties.timesteps[0]
dt = properties.dt[0]
boxlims = properties.boxlims[0]
numR1 = properties.numR1[0]
numR2 = properties.numR2[0]
R1 = properties.R1[0]
R2 = properties.R2[0]
phi = properties.phi[0]

steps = info.timesteps[0]
delta_t = info.dt[0]
runs = info.runs[0]
num_particles = info.num_particles[0]

#hello there

w = 10
h = 10
mrk1 = 2*R1*h*500/(9*boxlims)
mrk2 = 2*R2*h*500/(9*boxlims)

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

packingfractions = []
with open('packingfractions.csv') as file3:
    read_object = csv.reader(file3)
    for row in read_object:
        packingfractions.append(row)

v_rms_data = []
with open('rms_v.csv') as file4:
    read_obj = csv.reader(file4)
    for row in read_obj:
        v_rms_data.append(row)

for i in range(1,timesteps+1):
    for j in range(numR1):
        rowsR1[i][j] = float(rowsR1[i][j])

for i in range(1,timesteps+1):
    for j in range(numR2):
        rowsR2[i][j] = float(rowsR2[i][j])

for i in range(runs):
    packingfractions[0][i] = float(packingfractions[0][i])

for i in range(1,runs+1):
    for j in range(steps):
        v_rms_data[i][j] = float(v_rms_data[i][j])

time_data = []
for i in range(steps):
    time_data.append(i*delta_t)

#------------------------------------------------------------
jericho = 1
if jericho == 1:

    plt.figure(1, figsize=(w,h), dpi = 100)
    string1 = str(num_particles)
    string2 = " Particles"
    Title = "Root Mean Squared Velocity Over Time for an Increasing $\\phi$\n" + string1 + string2
    plt.title(Title)
    plt.yscale('log')
    plt.xscale('log')
    colours = ['#D8BFD8', '#DDA0DD', '#EE82EE', '#DA70D6', '#FF00FF', '#BA55D3', '#9370DB', '#8A2BE2', '#9400D3', '#9932CC', '#800080', '#4B0082']
    for i in range(1,runs+1):
        string = str(packingfractions[0][i-1])
        plt.plot(time_data, v_rms_data[i], label = string, color = colours[i-1])
    plt.xlabel("Time ($s$)")
    plt.ylabel("$\\sqrt{\\langle v^{2} \\rangle}$")
    plt.legend(loc = 'upper right', title = '$\\phi$')
    plt.savefig("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/rootmeansquarevelocity.png")
#------------------------------------------------------------
echo = 0
if echo == 1:

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

        x1 = x_posR1
        y1 = y_posR1
        x2 = x_posR2
        y2 = y_posR2

        plt.figure(figsize=(w,h), dpi = 100)
        plt.gca().set_xlim(0,boxlims)
        plt.gca().set_ylim(0,boxlims)
        plt.xticks([])
        plt.yticks([])
        plt.title("Particle Simulation - Repulsion Upon Overlap")
        plt.plot(x1, y1, marker = 'o', color = 'white', markeredgecolor = 'black', linestyle = 'none', markersize = mrk1)
        plt.plot(x2, y2, marker = 'o', color = 'white', markeredgecolor = 'red', linestyle = 'none', markersize = mrk2)
        label1 = "Box Dimensions: "
        label2 = str(boxlims)
        label3 = "\n"
        label4 = "Packing Fraction: "
        label5 = str(phi)
        label6 = "\n"
        label7 = "R1 = 1.0, R2 = 1.4\n"
        label8 = "Slide "
        iteration = str(i)
        label = label1 + label2 + label3 + label4 + label5 + label6 + label7 + label8 + iteration
        plt.xlabel(label)
        png = ".png"
        loc = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/slide"
        savepoint = loc + iteration + png
        plt.savefig(savepoint)






