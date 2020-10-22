from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation

# N = int(input("Insert number of objects N: ")) Could choose to use this
# Maybe we can make animation for mercury sun and jupyter earth sun

infile1 = open("./results/position_x.txt", "r")
infile2 = open("./results/position_y.txt", "r")
infile3 = open("./results/position_z.txt", "r")

x = np.loadtxt(infile1)
y = np.loadtxt(infile2)
z = np.loadtxt(infile3)
planets = np.genfromtxt("./results/planet_names.txt",dtype='str') # order given in main

pos =  np.array([x[:,0],y[:,0],z[:,0]]); pos1 = np.array([x[:,1],y[:,1],z[:,1]])
#pos2 = np.array([x[:,2],y[:,2],z[:,2]]); Why is jupyter number two here?
pos3 = np.array([x[:,3],y[:,3],z[:,3]])
pos4 = np.array([x[:,4],y[:,4],z[:,4]]);

#pos5 = np.array([x[:,5],y[:,5],z[:,5]])

pos6 = np.array([x[:,6],y[:,6],z[:,6]]);
"""
pos7 = np.array([x[:,7],y[:,7],z[:,7]])
pos8 = np.array([x[:,8],y[:,8],z[:,8]]); pos9 = np.array([x[:,9],y[:,9],z[:,9]])
"""

# initiate figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Setting the axes properties so they are constant for each frame
#ax.set_xlim3d([-5.0, 40.0])
ax.set_xlim3d([-1.5, 1.5])
ax.set_xlabel('X [AU]')

#ax.set_ylim3d([-30.0, 40.0])
ax.set_ylim3d([-1.5, 1.5])
ax.set_ylabel('Y [AU]')

#ax.set_zlim3d([-15.0, 10.0])
ax.set_zlim3d([-1.5, 1.5])
ax.set_zlabel('Z [AU]')

"""
#line.set_data(data[:2, :num])
#line.set_3d_properties(data[2, :num])
"""

# This is the first frame, and it is a point in 3D for each planet
line, = ax.plot(pos[0,0:1],pos[1,0:1],pos[2,0:1],'.', color = 'orange', label = planets[0])
line1, = ax.plot(pos1[0,0:1],pos1[1,0:1],pos1[2,0:1],'.',label = planets[1])
#line2, = ax.plot(pos2[0,0:1],pos2[1,0:1],pos2[2,0:1],label = planets[2])
line3, = ax.plot(pos3[0,0:1],pos3[1,0:1],pos3[2,0:1],'.',label = planets[3])
line4, = ax.plot(pos4[0,0:1],pos4[1,0:1],pos4[2,0:1],'.', label = planets[4])
#line5, = ax.plot(pos5[0,0:1],pos5[1,0:1],pos5[2,0:1],label = planets[5])
line6, = ax.plot(pos6[0,0:1],pos6[1,0:1],pos6[2,0:1],'.',label = planets[6])
"""
line7, = ax.plot(pos7[0,0:1],pos7[1,0:1],pos7[2,0:1],label = planets[7])
line8, = ax.plot(pos8[0,0:1],pos8[1,0:1],pos8[2,0:1],label = planets[8])
line9, = ax.plot(pos9[0,0:1],pos9[1,0:1],pos9[2,0:1],label = planets[9])
"""
plt.legend()

def animation_frame(num):
    # Update lines  with x,y and z data
    line.set_data(pos[:2, :num])         # x,y axis
    line.set_3d_properties(pos[2, :num]) # z-axis
    line1.set_data(pos1[:2, :num]); line1.set_3d_properties(pos1[2, :num])
    #line2.set_data(pos2[:2, :num]); line2.set_3d_properties(pos2[2, :num])
    line3.set_data(pos3[:2, :num]); line3.set_3d_properties(pos3[2, :num])
    line4.set_data(pos4[:2, :num]); line4.set_3d_properties(pos4[2, :num])
    #line5.set_data(pos5[:2, :num]); line5.set_3d_properties(pos5[2, :num])
    line6.set_data(pos6[:2, :num]); line6.set_3d_properties(pos6[2, :num])
    """
    line7.set_data(pos7[:2, :num]); line7.set_3d_properties(pos7[2, :num])
    line8.set_data(pos8[:2, :num]); line8.set_3d_properties(pos8[2, :num])
    line9.set_data(pos9[:2, :num]); line9.set_3d_properties(pos9[2, :num])
    """
    return [line,line1,line3,line4, line6]#,line4]#,line5,line6,line7,line8,line9]

N_points = len(x[:,1]);

#ani = animation.FuncAnimation(fig, func = animation_frame, N, fargs=(data, line), interval=10000/N, blit=False)

animation = FuncAnimation(fig, func = animation_frame, frames = range(N_points),interval = 150)#, #blit=True)
#animation.save('matplot003.gif', writer='imagemagick')
plt.show()
