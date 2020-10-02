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
  m_force



public:
  void initialize(m_name)   // Use keys for each planet
  void grav_force(m_name1,m_)         // calculate force between stellar objects
  void Verlet()             // Verlet solver
  void EulerChromer()       // EulerChromer solver
  void RungeKutta4()        // RungeKutta4 solver


};



#endif
