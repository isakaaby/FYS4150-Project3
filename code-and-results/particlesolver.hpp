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
  vec m_X,m_Y,m_Z;
  vec m_Vx,m_Vy, m_Vz;
  vec m_ax,m_ay,m_az;
  vec m_xupdate, m_yupdate; // params for RK4
  double m_beta;
  vec m_masses;

  int m_N;          // number of planets
  int m_k;          // number of time steps
  double m_T, m_h;
  double m_T0 = 0.0;

public:                          // general solver
  void initialize(double m_beta,int m_N, int k, int m_T);      // Use keys for each planet
  double force_a(vec pos, int l, int j);
  void verlet(double f(double s, double x, double y));                // Verlet solver
  void eulerchromer();          // EulerChromer solver
  void RK4(double f1(double t, double x, double y, double v),double force(double t, double x, double y, double v));           // RungeKutta4 solver
  void RK4_xupdate(double t, double x, double y, double v, double f1(double t, double x,  double y, double v), double f2(double t, double x, double y, double v));
  void RK4_yupdate(double t, double x, double y, double v, double f1(double t, double x,  double y, double v), double f2(double t, double x, double y, double v));
};

//subclass to solve planet case
class PlanetSolver : public ParticleSolver {
private:

public:
  void init(vector <string> m_name);            //init special solver for planet case
  void solvesystem();      //  solve for planet system

};



#endif
