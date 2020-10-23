import numpy as np

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")

x = np.loadtxt(infile1)
y = np.loadtxt(infile2)

tol = 1e-14

x_mercury = x[:,1]
y_mercury = y[:,1]


r2 = np.sqrt(x_mercury**2 + y_mercury**2)


idx = np.where(abs(r2-np.min(r2)) < tol)


print(idx[0][-1])

I = idx[0][-1]





theta_rad = np.arctan(y_mercury[I]/x_mercury[I])
print(theta_rad)


theta_arc = (4.848e6)*theta_rad

print(theta_arc)
