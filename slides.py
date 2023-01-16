import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from IPython.display import display

properties = pd.read_csv('properties.csv')

info = pd.read_csv('info.csv')

brownianinfo = pd.read_csv('brownianinfo.csv')

browniandisplacementinfo = pd.read_csv('browniandisplacementinfo.csv')

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

timesteps_brownian = brownianinfo.timesteps[0]
dt_brownian = brownianinfo.dt[0]
boxlims_brownian = brownianinfo.boxlims[0]
radius_brownian = brownianinfo.R1[0]
N_movie = brownianinfo.N[0]
N_msd = browniandisplacementinfo.N[0]
timesteps_msd = browniandisplacementinfo.timesteps[0]
dt_msd = browniandisplacementinfo.dt[0]
std_dev = browniandisplacementinfo.std_dev[0]

w = 10
h = 10
mrk1 = 2*R1*h*500/(9*boxlims)
mrk2 = 2*R2*h*500/(9*boxlims)

markersz = 2*radius_brownian*h*500/(9*boxlims_brownian)

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

mean_square_data = []
with open('mean_square_displacement.csv') as file6:
    read_obj_msd = csv.reader(file6)
    for row in read_obj_msd:
        mean_square_data.append(row)

rows_brownian_data = []
with open('brownian_positions.csv') as file5:
    reader_obj_brownian = csv.reader(file5)
    for row in reader_obj_brownian:
        rows_brownian_data.append(row)


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

for i in range(1):
    for j in range(timesteps_msd):
        mean_square_data[i][j] = float(mean_square_data[i][j])

for i in range(1,timesteps_brownian+1):
    for j in range(N_movie):
        rows_brownian_data[i][j] = float(rows_brownian_data[i][j])

time_data = []
for i in range(steps):
    time_data.append(i*delta_t)

time_data_msd = []
for i in range(timesteps_msd):
    time_data_msd.append(i*dt_msd)

#------------------------------------------------------------
# Plot Root Mean Squared Speed of Soft Repulsive Particles

jericho = 1
if jericho == 1:

    #from 10^1 onwards, find average rms_v across all packing fractions at regular intervals
    #use power law decay fit to find \beta
    #

    plt.figure(1, figsize=(w,h), dpi = 100)
    string1 = str(num_particles)
    string2 = " Particles"
    Title = "Root Mean Squared Speed Over Time for an Increasing $\\phi$\n" + string1 + string2
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
# Plot Soft Repulsion Movie

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
        loc = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/soft_repulsion_slide"
        savepoint = loc + iteration + png
        plt.savefig(savepoint)

#---------------------------------------------------------------------------
# Plot Mean Square Displacement of Brownian Particles

yotta = 0
if yotta == 1:

    model_msd = []
    for i in range(timesteps_msd):
        model = 2*time_data_msd[i]*std_dev
        model_msd.append(model)
    
    plt.figure(1, figsize=(w,h), dpi = 100)
    string1 = str(N_msd)
    string2 = " Particles"
    Title = "Root Mean Squared Displacement Over Time for " + string1 + string2
    plt.title(Title)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(time_data_msd, mean_square_data[0], 'k')
    plt.plot(time_data_msd, model_msd, 'b--', label = 'Model Data')
    plt.legend(loc = 'upper left')
    plt.xlabel("Time ($s$)")
    plt.ylabel("$\\sqrt{\\langle x^{2} \\rangle}$")
    plt.savefig("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/rootmeansquaredisplacement.png")

#---------------------------------------------------------------------------
# Plot Brownian Motion Movie

alpha = 0
if alpha == 1:

    for i in range(1,timesteps_brownian+1):

            x_pos = []
            y_pos = []
            for j in range(N_movie):
                if j % 2 == 0:
                    x_pos.append(rows_brownian_data[i][j])
                else:
                    y_pos.append(rows_brownian_data[i][j])

            x = x_pos
            y = y_pos

            plt.figure(figsize=(w,h), dpi = 100)
            plt.gca().set_xlim(0,boxlims_brownian)
            plt.gca().set_ylim(0,boxlims_brownian)
            plt.xticks([])
            plt.yticks([])
            plt.title("Particle Simulation - Brownian")
            plt.plot(x, y, marker = 'o', color = 'white', markeredgecolor = 'blue', linestyle = 'none', markersize = markersz)
            label1 = "Box Dimensions: "
            label2 = str(boxlims_brownian)
            label6 = "\n"
            label7 = "R = 1.0\n"
            label8 = "Slide "
            iteration = str(i)
            label = label1 + label2 + label6 + label7 + label8 + iteration
            plt.xlabel(label)
            png = ".png"
            loc = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/brownian_slide"
            savepoint = loc + iteration + png
            plt.savefig(savepoint)






