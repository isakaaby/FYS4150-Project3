#ifndef PARTICLESOLVER_HPP
#define PARTICLESOLVER_HPP

#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

// Setting up solarsystem as a superclass
// Work in units: time [years], distance [AU], mass [solar masses]

class ParticleSolver {

private:
  int K;

protected:
  string m_name;
  vec m_X,m_Y,m_m_Z;
  vec m_Vx,m_Vy, m_Vx_prev, m_Vy_prev;
  vec m_ax, m_ax_prev, m_ay,m_ay_prev;

  double m_force;
  int m_N;         // number of planets
  int m_k;          // number of time steps

public:                          // general solver
  void initialize(int m_N, int k);      // Use keys for each planet
  void Verlet();                // Verlet solver
  void EulerChromer();          // EulerChromer solver
  void RungeKutta4();           // RungeKutta4 solver
};

//subclass to solve planet case
class PlanetSolver : public ParticleSolver {
private:

public:
  void init();            //init special solver for planet case
  void solvesystem();      //  solve for planet system

};



#endif
