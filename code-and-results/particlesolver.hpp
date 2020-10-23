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

protected:
  vec m_X,m_Y,m_Z;
  vec m_Vx,m_Vy, m_Vz;
  vec m_ax,m_ay,m_az;
  vec m_xupdate, m_yupdate; // params for RK4
  double m_beta;
  vec m_masses;
  vec m_x0,m_y0,mz0;
  vec m_vx0,m_vy0,m_vz0;
  int m_N;          // number of planets
  int m_k;          // number of time steps
  double m_T, m_h, hh, m_T0;
  double M;
  vec m_Etot, m_Lx, m_Ly, m_Lz;
  vector<string> m_names;

public:                          // general solver
  void initialize(double beta, int N, int k, double T);      // Use keys for each planet
  void force_a(int l, int j);
  void verlet_pos(int l, int j);                // Verlet solver
  void verlet_vel(int l, int j);                // Verlet solver

  void eulerchromer(int l, int j);          // EulerChromer solver
  void forwardeuler(int l, int j);          // Forward euler solver
  double kinetic_energy(int i, int j);
  double potential_energy(double r, int l, int i, int j);
  double angular_momentum(double pos1, double v1, double pos2, double v2);
  void get_angular_momentum(int j);

};

//subclass to solve planet case
class PlanetSolver : public ParticleSolver {
private:

public:
  void init(vector<string> names, double beta, int N, int k, double T); //init special solver for planet case
  void init_sun_center(vector<string> names, double beta, int N, int k, double T, bool check);
  void solvesystem(bool check,int method);                         //  solve for planet system
  void write_pos_to_file();
  void write_vel_to_file();
  void test_constant_energy(double tol);
  void test_constant_angular(double tol);
  void test_circular_orbit(double tol);
  int random_index_generator(int min, int max);
  void test_convergence(vector<string> names,double beta, int N,int k, double T, int N_experiments, int method);
};


class MercurySunSolver : public ParticleSolver {

public:
  void init(vector<string> names, double beta, int N, int k, double T);           //init special solver for planet case
  void force_mercury_rel(int l, int j);
  void solve_mercury_sun_verlet();                          //  solve for planet system
  void solve_mercury_sun_eulerchromer();
  void write_pos_to_file();
};


#endif
