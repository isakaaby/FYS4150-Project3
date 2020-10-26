# FY4150-Project3 - Exploring the Solar System: Planetary Orbits and ODESolvers
Git for Project 3 in Computational Physics (FYS4150) 

### Main overview
In this project we have built a model for different solar systems using differential equations: forward Euler, Euler-Cromer and the Verlet method. The project description can be found [here](https://github.com/isakaaby/FYS4150-Project3/blob/master/Report/Project3.pdf).


### Code: Link and description of programs:
- [main.cpp](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/main.cpp) is the file containing the main function used to run the program, fitted with a menu to terminal for input to solve the different problems.

- [makefile](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/makefile) is the file through which all programs are run.

- [mercurysunsolver.cpp](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/mercurysunsolver.cpp) is used to model a solar system consisting of the Sun and Mercury, using a relativistic correction to the newtonian force to find perihelion precession.

- [particlesolver.cpp](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/particlesolver.cpp) is the "main" file consisting all the ODEs along with methods for finding the gravitational force and the potential, kinetic and total energy of the system.

- [planets.cpp](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/planets.cpp) sets up planet objects for the desired solar system giving the objects an initial position and velocity.

- [planetsolver.cpp](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/planetsolver.cpp) creates a solar system and adjusts it according to the center of mass and velocity, and tracks the trajectory while writing out position and velocity to separate files.

- [plot.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plot.py) makes a 2D plot for a chosen number of planets. 

- [plot3d.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plot3d.py) makes a 3D plot for a chosen number of planets. 

- [plotenergy.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plotenergy.py) makes a plot of the kinetic, potential and total energy. 

- [animate.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/animate.py) makes an animation of the planets with data from NASA.

- [get_angle.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/get_angle.py) calculates the derivative and double derivatives in order to find minima in position from the Sun. This is applied to find the perihelion precession of Mercury.

The files can be compiled with make all. 

How to run the programmes to reproduce the results discussed in the article: The menu gives 5 alternatives, where Forward Euler, Euler Cromer and Velocity Verlet (Störmer) can be chosen by changing input in main.cpp.
  1.  Run for Earth-Sun system: This method was used with inital conditions x = 1, y = 0, z = 0, vx = 0, vy = 6.28, vz = 0 for circular orbits. The same inital conditions are used for elliptical orbits, but with vy = 5. This method was also used to test a non-square inverse gravitational law (beta -> 3). It also tests if energy is constant, constant angular momentum and if the orbit is circular.
  2. Run for Earth-Jupiter-Sun system: This method was used to investigate the effects of Jupiter on the Earth with masses 1,100 and 1000 times of Jupiter's original mass. The planets were initialized with x = 1, y = 0, z = 0, vx = 0, vy = 6.28, vz = 0 for the Earth and x = 6.2, y = 0, z = 0, vx = 0, vy = 2.76 and vz = 0. It also tests if energy is constant, constant angular momentum and if the orbit is circular.
  3. Run for all planets: Solves the entire solar system with data from NASA. 
  4. Run for the Mercury Sun system: Solves for the case with Mercury around the Sun. Used together with get_angle.py, this is used to get the perihelion precession of Mercury.
  5. Calculates a convergence test for the Eart-Sun system with orbit of 1 AU (circular) as exact solution. Choose number of experiments as terminal input and wanted method in main.




### Links and packages
- The NASA data used to initialize some of the systems was pulled from [here.](https://ssd.jpl.nasa.gov/horizons.cgi#top)

- Documentation for Matplotlib from python from [here](https://matplotlib.org/)

