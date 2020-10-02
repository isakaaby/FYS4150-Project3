#ifndef SolarSystem_HPP
#define SolarSystem_HPP

#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

// Setting up solarsystem as a superclass

class SolarSystem {

private:


protected:
  string m_name
  vec m_x,_y,m_m_z
  double m_force
  int m_N           // number of planets



public:
  void initialize(m_N)      // Use keys for each planet
  void grav_force()         // calculate force between stellar objects
  void Verlet()             // Verlet solver
  void EulerChromer()       // EulerChromer solver
  void RungeKutta4()        // RungeKutta4 solver


};



#endif
