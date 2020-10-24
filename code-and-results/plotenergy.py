import matplotlib.pyplot as plt
import numpy as np


infile1 = open("./results/potential.txt", "r")
infile2 = open("./results/kinetic.txt", "r")
infile3 = open("./results/total.txt", "r")

x = np.loadtxt(infile1.readlines()[:-1])
y = np.loadtxt(infile2.readlines()[:-1])
z = np.loadtxt(infile3.readlines()[:-1])

t = np.linspace(1, len(x), len(x))


#plt.plot(t, x[:, 1])
#plt.plot(t, y[:, 1])
#plt.plot(t, z[:, 0])

for i in range(0, len(x[1, :])):
    #plt.plot(t, x[:, i])
    #plt.plot(t, y[:, i])
    plt.plot(t, z[:, i])
plt.show()
