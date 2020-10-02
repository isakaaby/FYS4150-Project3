#ifndef SolarSystem_HPP
#define SolarSystem_HPP

#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

// Setting up solarsystem as a superclass
// Work in units: time [years], distance [AU], mass [solar masses]

class SolarSystem {

private:

protected:
  string m_name
  vec m_x,_y,m_m_z
  double m_force
  int m_N                  // number of planets
  int m_k                  // number of time steps
  double m_G = 39.478;     //AU^3 yr^-2 M_sol^-1




public:
  void initialize(m_N,m_k)      // Use keys for each planet
  void grav_force()         // calculate force between stellar objects
  void Verlet()             // Verlet solver
  void EulerChromer()       // EulerChromer solver
  void RungeKutta4()        // RungeKutta4 solver


};



#endif
