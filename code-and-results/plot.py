import numpy as np
import matplotlib.pyplot as plt

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

k = 10000
N = 2

x_1 = np.zeros(k)
x_2 = np.zeros(k)
x_3 = np.zeros(k)
x_4 = np.zeros(k)
x_5 = np.zeros(k)
x_6 = np.zeros(k)
x_7 = np.zeros(k)
x_8 = np.zeros(k)
x_9 = np.zeros(k)
x_10 = np.zeros(k)

y_1 = np.zeros(k)
y_2 = np.zeros(k)
y_3 = np.zeros(k)
y_4 = np.zeros(k)
y_5 = np.zeros(k)
y_6 = np.zeros(k)
y_7 = np.zeros(k)
y_8 = np.zeros(k)
y_9 = np.zeros(k)
y_10 = np.zeros(k)

i = 0
for line in infile1:
    numbers1 = line.split()
    x_1[i] = float(numbers1[0])
    x_2[i] = float(numbers1[1])
    x_3[i] = float(numbers1[2])
    x_4[i] = float(numbers1[3])
    x_5[i] = float(numbers1[4])
    x_6[i] = float(numbers1[5])
    x_7[i] = float(numbers1[6])
    x_8[i] = float(numbers1[7])
    x_9[i] = float(numbers1[8])
    x_10[i] = float(numbers1[9])
    i = i + 1

i = 0
for line in infile2:
    numbers2 = line.split()
    y_1[i] = float(numbers2[0])
    y_2[i] = float(numbers2[1])
    y_3[i] = float(numbers2[2])
    y_4[i] = float(numbers2[3])
    y_5[i] = float(numbers2[4])
    y_6[i] = float(numbers2[5])
    y_7[i] = float(numbers2[6])
    y_8[i] = float(numbers2[7])
    y_9[i] = float(numbers2[8])
    y_10[i] = float(numbers2[9])
    i = i + 1


plt.plot(x_1, y_1)
plt.plot(x_2, y_2)
plt.plot(x_3, y_3)
plt.plot(x_4, y_4)
plt.plot(x_5, y_5)
plt.plot(x_6, y_6)
plt.plot(x_7, y_7)
plt.plot(x_8, y_8)
plt.plot(x_9, y_9)
plt.plot(x_10, y_10)
plt.show()
