import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

N = int(input("Insert number of objects N: "))

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

x = np.loadtxt(infile1)
y = np.loadtxt(infile2)
z = np.loadtxt(infile3)

planets = ["Sun", "Earth", "Jupiter", "Mars", "Venus", "Saturn", "Mercury", "Uranus", "Neptune", "Pluto"]

for i in range(0, N):
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax.plot(x[:, i], y[:, i], z[:, i])
    plt.plot(x[:,i], y[:, i], label = planets[i])
    plt.xlabel('x (AU)')
    plt.ylabel('y (AU)')
    plt.legend()
plt.show()
