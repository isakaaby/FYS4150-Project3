import numpy as np

#Reading in posiion values from textfiles.
infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")

#Assigning variable names to file input.
x = np.loadtxt(infile1)
y = np.loadtxt(infile2)

tol = 1e-14

#Extracting Mercury's position.
x_mercury = x[:,1]
y_mercury = y[:,1]

#Computing the r-term.
r2 = np.sqrt(x_mercury**2 + y_mercury**2)

#Finding position closest to the sun.
idx = np.where(abs(r2-np.min(r2)) < tol)


print(idx[0][-1])

I = idx[0][-1]




#Computing the perihelion precession.
theta_rad = np.arctan(y_mercury[I]/x_mercury[I])
print(theta_rad)


theta_arc = (4.848e6)*theta_rad

print(theta_arc)
