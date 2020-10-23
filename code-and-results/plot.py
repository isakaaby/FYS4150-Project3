import matplotlib.pyplot as plt
import numpy as np

N = int(input("Insert number of objects N: "))

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

x = np.loadtxt(infile1)
y = np.loadtxt(infile2)
z = np.loadtxt(infile3)
planets = np.genfromtxt("./results/planet_names.txt",dtype='str')

if x[0,0] == 0.0:
    plt.plot(x[:,0], y[:,0], ".", color = "orange", label = "Sun")
    for i in range(1, N):
        plt.plot(x[:, i], y[:, i], label = planets[i])
        plt.xlabel('x (AU)',fontsize = 13)
        plt.ylabel('y (AU)',fontsize = 13)
        plt.legend(loc = "upper right",fontsize = 15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    plt.show()

else:
    for i in range(0, N):
        plt.plot(x[:, i], y[:, i], label = planets[i])
        plt.xlabel('x (AU)',fontsize = 13)
        plt.ylabel('y (AU)',fontsize = 13)
        plt.legend(loc = "upper right",fontsize = 15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    plt.show()
