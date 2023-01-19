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
    element = (24*eps*pow(sig,6)*(2*pow(sig,6)-pow(r[i],6)))/pow(r[i],13)
    model_F.append(element)
    if element > max_F:
        max_F = element
        max_F_r = r[i]
for i in range(len(r)):
    element = 4*eps*(   (sig/(r[i]))**12 - (sig/(r[i]))**6      )
    V.append(element)


# V_h (r) = eps_h [  (sig_h/r-2R)^v - ( 1 - (r-2R)(v/sig_h) ) ]

# short ranged steep repulsive potential to prevent particle overlaps

# potential holds for 2R <= r <= (2R + sig_h), V=0 otherwise

# sig_h = (2^{1/6} - 1)2R         such that the onset of this steep
# repulsive interaction at particle separations r = 2R + sig_h coincides
# with the minimum in the Lennard-Jones potential at r = 2^{1/6} sig_LJ =
# at r = 2^{7/6}R




label_max_F = 'Max Force = ' + str(round(max_F,2))


ydata = [ylow,yhigh]
xdata = [D,D]

plt.figure()
plt.gca().set_ylim(ylow,yhigh)
plt.gca().set_xlim(0,3)
plt.grid()
plt.plot(xdata,ydata,'b--')
plt.plot(r, model_F, "k--", label = 'LJ Force')
plt.plot(r, V, "r--", label = 'LJ Potential')

plt.text(D, txtshift, '$r = R_{i} + R_{j}$')
plt.scatter(max_F_r, max_F, color = 'b')
plt.text(max_F_r, 1.1*max_F, label_max_F)
plt.legend(loc = 'upper left')
plt.xlabel('$r$')
plt.ylabel('Force / Potential')
plt.savefig("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/LJpotential.png")


#---------------------------------------------------------------------------------------------

LJproperties = pd.read_csv('LJproperties.csv')

LJtimesteps = LJproperties.timesteps[0]
LJboxlims = LJproperties.boxlims[0]
LJnumR1 = LJproperties.numR1[0]
LJnumR2 = LJproperties.numR2[0]
LJR1 = LJproperties.R1[0]
LJR2 = LJproperties.R2[0]
LJphi = LJproperties.phi[0]

dts_LJ = []
with open('dts_LJ.csv') as file1:
    read_object_dts_LJ = csv.reader(file1)
    for row in read_object_dts_LJ:
        dts_LJ.append(row)

time_data_LJ = [float(dts_LJ[0][0])]
for i in range(1,LJtimesteps):
    prev_time_element = float(time_data_LJ[i-1])
    del_t = float(dts_LJ[0][i])
    element = prev_time_element + del_t
    time_data_LJ.append(element)

rowsR1 = []
with open('LJR1.csv') as file2:
    reader_obj = csv.reader(file2)
    for row in reader_obj:
        rowsR1.append(row)

rowsR2 = []
with open('LJR2.csv') as file3:
    reader_object = csv.reader(file3)
    for row in reader_object:
        rowsR2.append(row)

for i in range(1,LJtimesteps+1):
    for j in range(LJnumR1):
        rowsR1[i][j] = float(rowsR1[i][j])

for i in range(1,LJtimesteps+1):
    for j in range(LJnumR2):
        rowsR2[i][j] = float(rowsR2[i][j])

#-----------------
#start plotting

w = 10
h = 10
mrk1 = 2*LJR1*h*500/(9*LJboxlims)
mrk2 = 2*LJR2*h*500/(9*LJboxlims)


for i in range(750,1000): # 1, timesteps + 1

        LJ_dt = str(dts_LJ[0][i-1])

        x_posR1 = []
        y_posR1 = []
        x_posR2 = []
        y_posR2 = []
        for j in range(LJnumR1):
            if j % 2 == 0:
                x_posR1.append(rowsR1[i][j])
            else:
                y_posR1.append(rowsR1[i][j])
        for j in range(LJnumR2):
            if j % 2 == 0:
                x_posR2.append(rowsR2[i][j])
            else:
                y_posR2.append(rowsR2[i][j])

        x1 = x_posR1
        y1 = y_posR1
        x2 = x_posR2
        y2 = y_posR2

        plt.figure(figsize=(w,h), dpi = 100)
        plt.gca().set_xlim(0,LJboxlims)
        plt.gca().set_ylim(0,LJboxlims)
        plt.xticks([])
        plt.yticks([])
        plt.title("Particle Simulation - Lennard-Jones Potential")
        plt.plot(x1, y1, marker = 'o', color = 'white', markeredgecolor = 'red', linestyle = 'none', markersize = mrk1)
        plt.plot(x2, y2, marker = 'o', color = 'white', markeredgecolor = 'red', linestyle = 'none', markersize = mrk2)
        label1 = "Box Dimensions: "
        label2 = str(LJboxlims)
        label3 = "\n"
        label31 = "Timestep (s) Between Now & Next Slide = "
        label32 = LJ_dt
        label33 = "\n"
        label4 = "Packing Fraction: "
        label5 = str(LJphi)
        label6 = "\n"
        Rlabel = "R1 = " + str(LJR1) + ", R2 = " + str(LJR2) + "\n"
        label7 = Rlabel
        label8 = "Slide "
        iteration = str(i)
        label = label1 + label2 + label3 + label31 + label32 + label33 + label4 + label5 + label6 + label7 + label8 + iteration
        plt.xlabel(label)
        png = ".png"
        loc = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/Lennard-Jones/LJ_slide"
        savepoint = loc + iteration + png
        plt.savefig(savepoint)