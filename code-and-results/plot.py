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

#plot_3D = string(input("Plot in 3D? (yes/no):"))


if x[0,0] == 0.0:
    plt.plot(x[:,0], y[:,0], ".", label = "Sun")
    for i in range(1, N):
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.plot(x[:, i], y[:, i], z[:, i])
        plt.plot(x[99500:, i], y[99500:, i], label = planets[i])
        plt.plot(x[:350, i], y[:350, i], label = 'merkur etter 1 år')
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.legend()
    plt.show()

else:
    for i in range(0, N):
        plt.plot(x[:, i], y[:, i], label = planets[i])
        plt.xlabel('x (AU)')
        plt.ylabel('y (AU)')
        plt.legend()
    plt.show()
