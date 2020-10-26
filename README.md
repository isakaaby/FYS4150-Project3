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

- [plot3d.py] (https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plot3d.py) makes a 3D plot for a chosen number of planets.

- [plotLE.py] (https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plotLE.py) makes 

- [plotenergy.py] (https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/plotenergy.py)

- [animate.py] (https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/animate.py)

- [get_angle.py](https://github.com/isakaaby/FYS4150-Project3/blob/master/code-and-results/get_angle.py)
