import numpy as np
import matplotlib.pyplot as plt

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")


k = 10000
N = 2
"""
x1 = np.zeros(k)
y1 = np.zeros(k)

x2 = np.zeros(k)
y2 = np.zeros(k)


xlines = infile1.readlines()
ylines = infile2.readlines()
zlines = infile3.readlines()

x = []
y = []
z = []

L = len(xlines)
for i in range(L):
    xvals = xlines[i].split()
    yvals = ylines[i].split()
    zvals = zlines[i].split()

    x.append(float(xvals[0]))
    x.append(float(xvals[1]))

    y.append(float(yvals[0]))
    y.append(float(yvals[1]))

    z.append(float(zvals[0]))
    z.append(float(zvals[1]))
"""
x_sun = np.zeros(k)
x_earth = np.zeros(k)

i = 0
for line in infile1:
    numbers1 = line.split()
    x_sun[i] = float(numbers1[0])
    x_earth[i] = float(numbers1[1])
    i = i + 1

y_sun = np.zeros(k)
y_earth = np.zeros(k)

j = 0
for line in infile2:
    numbers2 = line.split()
    y_sun[j] = float(numbers2[0])
    y_earth[j] = float(numbers2[1])
    j = j + 1

plt.plot(x_earth,y_earth)
plt.plot(x_sun,y_sun)
plt.show()
