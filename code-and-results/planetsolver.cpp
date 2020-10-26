#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

void PlanetSolver::init(vector<string> names, double beta, int N, int k, double T){
  initialize(beta, N, k, T);

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);

  m_masses = zeros<vec>(m_N);
  M = 0;
  double posMx = 0;
  double posMy = 0;
  double posMz = 0;
  double velMx = 0;
  double velMy = 0;
  double velMz = 0;

  double velx_sun = 0;
  double vely_sun = 0;
  double velz_sun = 0;

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
    m_X(i*m_k) = params(1); m_Y(i*m_k) = params(2); m_Z(i*m_k) = params(3);
    m_Vx(i*m_k) = 365*params(4); m_Vy(i*m_k) = 365*params(5); m_Vz(i*m_k) = 365*params(6);
    m_ax(i*m_k) = m_ay(i*m_k) = m_az(i*m_k) = 0.0;
    M += m_masses[i];
    posMx += m_masses[i]*m_X(i*m_k);
    posMy += m_masses[i]*m_Y(i*m_k);
    posMz += m_masses[i]*m_Z(i*m_k);
    velMx += m_masses[i]*m_Vx(i*m_k);
    velMy += m_masses[i]*m_Vy(i*m_k);
    velMz += m_masses[i]*m_Vz(i*m_k);
    if (i != 0){
      velx_sun -= m_masses[i]*m_Vx(i*m_k);
      vely_sun -= m_masses[i]*m_Vy(i*m_k);
      velz_sun -= m_masses[i]*m_Vz(i*m_k);
    }
  }
  m_Vx(0) = velx_sun;
  m_Vy(0) = vely_sun;
  m_Vz(0) = velz_sun;
  for (int i = 0; i < m_N; i++){
    m_X(i*m_k) -= posMx/M;
    m_Y(i*m_k) -= posMy/M;
    m_Z(i*m_k) -= posMz/M;
    m_Vx(i*m_k) -= velMx/M;
    m_Vy(i*m_k) -= velMy/M;
    m_Vz(i*m_k) -= velMz/M;
  }
};

void PlanetSolver::init_sun_center(vector<string> names, double beta, int N, int k, double T, bool test_convergence = false){
  initialize(beta, N, k, T);

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);

  m_masses = zeros<vec>(m_N);

  double posx0, posy0, posz0, velx0, vely0, velz0;
  if (test_convergence == false){
  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
    if (i != 0){
      cout << "x0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> posx0;
      cout << "y0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> posy0;
      cout << "z0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> posz0;
      cout << "vx0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> velx0;
      cout << "vy0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> vely0;
      cout << "vz0 for" << " " << m_names[i] << " " << "=" << " ";
      cin >> velz0;
      m_X(i*m_k) = posx0; m_Y(i*m_k) = posy0; m_Z(i*m_k) = posz0;
      m_Vx(i*m_k) = velx0; m_Vy(i*m_k) = vely0; m_Vz(i*m_k) = velz0;
      m_ax(i*m_k) = m_ay(i) = m_az(i*m_k) = 0.0;
      }
    }
  } else {
    for (int i = 0; i < m_N; i++){
      params = Planet.initialize(m_names[i]);
      m_masses(i) = params(0);
      }
      posx0 = 1.0 ; posy0 = 0.0; posz0 = 0.0;
      velx0 = 0.0; vely0 = 2*M_PI; velz0 = 0.0;

      int i = 1;
      m_X(i*m_k) = posx0; m_Y(i*m_k) = posy0; m_Z(i*m_k) = posz0;
      m_Vx(i*m_k) = velx0; m_Vy(i*m_k) = vely0; m_Vz(i*m_k) = velz0;
      m_ax(i*m_k) = m_ay(i) = m_az(i*m_k) = 0.0;
  }
};


void PlanetSolver::solvesystem(bool check, int method){
 int s;
 if(check == true) {
   s = 1;
 } else {
   s = 0;
 }

  if (method == 1){ //use velocity verlet
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        verlet_pos(i,j);
      }
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j+1);
        verlet_vel(i,j);
      }
    }
  }
  if (method == 2){
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j);
        forwardeuler(i,j);
        }
      }
    }

  if (method == 3){
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j);
        eulerchromer(i,j);
      }
    }
  }
}


int PlanetSolver::random_index_generator(int min, int max){
  // Using random generator from namespace std
  random_device seed;
  mt19937 rng(seed());
  uniform_int_distribution<int> uni(min,max);
  int random_integer = uni(rng);
  return random_integer;
}

void PlanetSolver::test_constant_energy(double tol){
  int j = random_index_generator(0,m_k);
  double diff;
  for (int i = 1; i < m_N; i++){
    diff = abs(m_Etot(i*m_k + 1) - m_Etot(i*m_k + j));
    if (diff < tol) {
      continue;
    } else {
    cout << "Error: Energy not conserved for celestial bodies with tolerance:" << " " << tol << "\n";
    break;
    }
  }
  cout << "Absolute error total energy:" << " " <<  diff << "\n";
}

void PlanetSolver::test_constant_angular(double tol){
  double diffL;
  // get angular momentum for all times
  for (int j = 0; j < m_k; j++){
    get_angular_momentum(j);
  }
  int j = random_index_generator(0,m_k);
  for (int i = 1; i < m_N; i++){
    diffL = abs((m_Lx(i*m_k + 1)*m_Lx(i*m_k + 1) + m_Ly(i*m_k + 1)*m_Ly(i*m_k + 1))  \
        - (m_Lx(i*m_k + j)*m_Lx(i*m_k + j) + m_Ly(i*m_k + j)*m_Ly(i*m_k + j)));
    if (diffL < tol) {
      continue;
    } else {
    cout << "Error: Angular momentum not conserved for celestial bodies with tolerance:" << " " << tol << "\n";
    break;
    }
  }
  cout << "Absolute error total angular momentum:" << " " << diffL << "\n";
}

void PlanetSolver::test_circular_orbit(double tol){
  int i = 1;
  int j = random_index_generator(0,m_k);
  double diffr = abs(1 - sqrt(m_X(i*m_k+j)*m_X(i*m_k+j) \
              + m_Y(i*m_k+j)*m_Y(i*m_k+j) + m_Z(i*m_k+j)*m_Z(i*m_k+j)));
  if (diffr > tol){
    cout << "Orbit is not circular with tolerance:" << " " << tol << "\n";
    }
    cout << "Relative error from a circular orbit: " << " " << diffr << "\n";
}

void PlanetSolver::test_convergence(vector<string> names,double beta, int N,int k, double T, int N_experiments, int method){
  // must solve for several dt's and check stability
  // --> need to call the solver in a loop
  // Check for conservation of energy og angular momentum?
  vec numpoints = zeros<vec>(N_experiments);
  vec E = zeros<vec>(N_experiments);
  vec h = zeros<vec>(N_experiments);
  vec r = zeros<vec>(N_experiments);
  r(N_experiments-1)= 0.0;

  bool sun_center = true;
  int steps = 0;
  double next_error;
  int i = 1;
  double factor = 5;
  bool test_convergence = true;
  while (steps < N_experiments){
    double error = 0;
    init_sun_center(names,beta,N,k,T,test_convergence);
    solvesystem(sun_center,method);
    // Calculate error and convergence rate
    for (int j = 0; j < m_k; j++){
      next_error =  1 - sqrt(m_X(i*m_k+j)*m_X(i*m_k+j) \
      + m_Y(i*m_k+j)*m_Y(i*m_k+j) + m_Z(i*m_k+j)*m_Z(i*m_k+j));
      error = error + next_error*next_error;
      }

    E(steps) = sqrt(m_h*error);
    h(steps) = m_h;
    numpoints(steps) = m_k;
    k = m_k*factor;
    steps += 1;
    }
    for (int j = 1; j < N_experiments; j++){
      r(j-1) = log(E(j)/E(j-1))/log(1/factor);
    }
    write_error_to_file(N_experiments,E,r,h,numpoints);
    cout << "Convergence rates:\n"  << r <<"\n";
    cout << "relative error:\n"  << E << "\n";
    cout << "step size:\n" << h << "\n";
}

void PlanetSolver::write_pos_to_file(){
  ofstream x;
  ofstream y;
  ofstream z;
  ofstream planet_names;

  string filename_1("./results/position_x.txt");
  string filename_2("./results/position_y.txt");
  string filename_3("./results/position_z.txt");
  string filename_4("./results/planet_names.txt");
  x.open(filename_1);
  y.open(filename_2);
  z.open(filename_3);
  planet_names.open(filename_4);

  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      x << setprecision(15) << m_X(i*m_k+j) << " ";
      y << setprecision(15) << m_Y(i*m_k+j) << " ";
      z << setprecision(15) << m_Z(i*m_k+j) << " ";
    }
    x << "\n";
    y << "\n";
    z << "\n";
  }
  for (int i = 0; i < m_N; i++){
    planet_names << m_names[i] << "\n";
  }
  x.close();
  y.close();
  z.close();
  planet_names.close();
}

void PlanetSolver::write_vel_to_file(){
  ofstream x;
  ofstream y;
  ofstream z;

  string filename_1("./results/velocity_x.txt");
  string filename_2("./results/velocity_y.txt");
  string filename_3("./results/velocity_z.txt");
  x.open(filename_1);
  y.open(filename_2);
  z.open(filename_3);
  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      x << m_X(i*m_k+j) << " ";
      y << m_Y(i*m_k+j) << " ";
      z << m_Z(i*m_k+j) << " ";
    }
    x << "\n";
    y << "\n";
    z << "\n";
  }
  x.close();
  y.close();
  z.close();
}

void PlanetSolver::write_error_to_file(int N_experiments, vec E, vec r, vec h, vec nump){
  ofstream error;

  string filename_1("./results/error_params.txt");
  error.open(filename_1);

  error << "E" << " "; error << "h" << " ";
  error << "r" << " "; error << "num_points" << " ";
  error << "\n";
  for (int j = 0; j < N_experiments; j++){
    error << E(j) << " ";
    error << h(j) << " ";
    error << r(j) << " ";
    error << nump(j) << " ";
    error << "\n";
    }
  error.close();

}
