#ifndef ParticleSolver_HPP
#define ParticleSolver_HPP

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
  vector<string> m_names;
  vec m_X,m_Y,m_m_Z;
  vec m_Vx,m_Vy, m_Vx_prev, m_Vy_prev;
  vec m_ax, m_ax_prev, m_ay,m_ay_prev;

  int m_N;          // number of planets
  int m_k;          // number of time steps
  double m_T, m_T0;

public:                          // general solver
  void initialize(int m_N, int k);      // Use keys for each planet
  void verlet(double f(double x, double y));                // Verlet solver
  void eulerchromer();          // EulerChromer solver
  void rungekutta4();           // RungeKutta4 solver
};

//subclass to solve planet case
class PlanetSolver : public ParticleSolver {
private:

public:
  void init(vector <string> m_name);            //init special solver for planet case
  void solvesystem();      //  solve for planet system

};



#endif
