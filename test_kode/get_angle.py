### Programme to investigate perihelion precession ###
import matplotlib.pyplot as plt
import numpy as np

#N = int(input("Insert number of objects N: "))

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

x = np.loadtxt(infile1)
y = np.loadtxt(infile2)
z = np.loadtxt(infile3)
planets = np.genfromtxt("./results/planet_names.txt",dtype='str')

# Get radii of Mercury from sun
r = np.sqrt(x[:,1]*x[:,1] + y[:,1]*y[:,1] + z[:,1]*z[:,1]);
#print(r)

# calculate differentials assuming evenly spacing in time
dr = np.gradient(r)
ddr = np.gradient(dr)

print(len(dr))
print(len(ddr))


#get indices where ddr positive (>0) and dr =(approx) 0.
perihelion_points = np.logical_and(ddr > 0, np.absolute(dr) < 0.43e-4)
perihelion_indices = np.where(perihelion_points == True)
print(len(perihelion_indices[0]))


indices = perihelion_indices[0]
r_mins =  r[indices]
x_p = x[indices,1]; y_p = y[indices,1];

angle_arcsec = np.arctan(y_p/x_p)*(4.848e6); #*206264.8062471

axis =  np.linspace(0,len(r_mins),len(r_mins))

anglefit = np.polyfit(axis,angle_arcsec,1)

print(anglefit)

def fit(x):
    return anglefit[0]*x + anglefit[1];


#print(fit(axis[1]))
plt.plot(axis,r_mins, label = 'r_mins')
plt.legend()
plt.show()
#plt.plot(axis,angle_arcsec)
#plt.show()
plt.plot(axis,fit(axis),label = 'angles-arcsec');
plt.legend()
plt.show()
