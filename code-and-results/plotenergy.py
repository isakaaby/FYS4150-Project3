import matplotlib.pyplot as plt
import numpy as np

#Reading in energy from textfiles.
infile1 = open("./results/potential.txt", "r")
infile2 = open("./results/kinetic.txt", "r")
infile3 = open("./results/total.txt", "r")

#Assigning variable names to file input.
x = np.loadtxt(infile1.readlines(), skiprows=1)
y = np.loadtxt(infile2.readlines(), skiprows=1)
z = np.loadtxt(infile3.readlines(), skiprows=1)

t = np.linspace(1, len(x), len(x))

#Plotting potential, kinetic and total energy as x, y, z.
for i in range(0, len(x[1, :])):
    plt.plot(t, x[:, i])
    plt.plot(t, y[:, i])
    plt.plot(t, z[:, i])
    plt.xlabel("Timestep (dt)")
    plt.ylabel("Energy (J)")
    plt.title("Earth's energy")
plt.show()
