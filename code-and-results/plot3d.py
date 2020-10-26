import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

N = int(input("Insert number of objects N: "))

#Reading in posiion values from textfiles.
infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

#Assigning variable names to file input.
x = np.loadtxt(infile1)
y = np.loadtxt(infile2)
z = np.loadtxt(infile3)
planets = np.genfromtxt("./results/planet_names.txt",dtype='str')  #Reading planet names from file.

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


## chose colour of axis
ax.tick_params(axis='x',colors='grey')
ax.tick_params(axis='y',colors='grey')
ax.tick_params(axis='z',colors='grey')

ax.w_xaxis.line.set_color("grey")
ax.w_yaxis.line.set_color("grey")
ax.w_zaxis.line.set_color("grey")

#Checking and plotting for the Sun.
if x[0,0] == 0.0:
    ax.plot(x[:,0], y[:,0], ".", color = "orange", label = "Sun")
    for i in range(1, N):
        ax.plot(x[:,i], y[:,i], z[:,i] , label = planets[i])
        ax.set_xlabel('x [AU]',color = "grey", fontsize = 15)
        ax.set_ylabel('y [AU]',color = "grey", fontsize = 15)
        ax.set_zlabel('z [AU]',color = "grey", fontsize = 15)
        plt.legend(loc = "upper left", fontsize = 10)
    plt.show()

#Plotting the rest of the planets.
else:
    for i in range(0, N):
        ax.plot(x[:,i], y[:,i], z[:,i], label = planets[i])
        ax.set_xlabel('x [AU]',color = "grey", fontsize = 15)
        ax.set_ylabel('y [AU]',color = "grey", fontsize = 15)
        ax.set_zlabel('z [AU]',color = "grey",fontsize = 15)
        plt.legend(loc = "upper left",fontsize = 10)
    plt.show()
