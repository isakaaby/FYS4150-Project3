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

#Checking and plotting for the Sun.
if x[0,0] == 0.0:
    plt.plot(x[:,0], y[:,0], ".", label = "Sun")
    for i in range(1, N):
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.plot(x[:, i], y[:, i], z[:, i])
        plt.plot(x[:, i], y[:, i], label = planets[i])
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.legend()
    plt.show()

#Plotting the rest of the planets.
else:
    for i in range(0, N):
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.plot(x[:, i], y[:, i], z[:, i])
        plt.plot(x[:, i], y[:, i], label = planets[i])
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.legend()
    plt.show()
