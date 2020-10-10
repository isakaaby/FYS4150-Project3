import numpy as np
import matplotlib.pyplot as plt

infile = open("./results/Earth_sun.txt", "r")

k = 100
N = 2

x1 = np.zeros(k)
y1 = np.zeros(k)

x2 = np.zeros(k)
y2 = np.zeros(k)


for line in infile:
    numbers = line.split()
    for j in range(k):
        x1[j] = float(numbers[j*N + 0])
