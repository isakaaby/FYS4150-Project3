#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

//Method that intializes a solver, returns nothing.
void PlanetSolver::init(vector<string> names, double beta, int N, int k, double T){
  initialize(beta, N, k, T); //Method that sets up initializer for the solver, returns nothing.

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();    //Reads in position and velocity from NASA data.
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
    params = Planet.initialize(m_names[i]);   //Creating planet objects.

    //Extracting planet information.
    m_masses(i) = params(0);
    m_X(i*m_k) = params(1); m_Y(i*m_k) = params(2); m_Z(i*m_k) = params(3);
    m_Vx(i*m_k) = 365*params(4); m_Vy(i*m_k) = 365*params(5); m_Vz(i*m_k) = 365*params(6);
    m_ax(i*m_k) = m_ay(i*m_k) = m_az(i*m_k) = 0.0;
    M += m_masses[i];

    //Updating position and velocity.
    posMx += m_masses[i]*m_X(i*m_k);
    posMy += m_masses[i]*m_Y(i*m_k);
    posMz += m_masses[i]*m_Z(i*m_k);
    velMx += m_masses[i]*m_Vx(i*m_k);
    velMy += m_masses[i]*m_Vy(i*m_k);
    velMz += m_masses[i]*m_Vz(i*m_k);

    if (i != 0){    //Checking if object is sun or not.
      velx_sun -= m_masses[i]*m_Vx(i*m_k);
      vely_sun -= m_masses[i]*m_Vy(i*m_k);
      velz_sun -= m_masses[i]*m_Vz(i*m_k);
    }
  }

  //Updating the sun's velocity.
  m_Vx(0) = velx_sun;
  m_Vy(0) = vely_sun;
  m_Vz(0) = velz_sun;

  //Computing the center of mass and velocity center.
  for (int i = 0; i < m_N; i++){
    m_X(i*m_k) -= posMx/M;
    m_Y(i*m_k) -= posMy/M;
    m_Z(i*m_k) -= posMz/M;
    m_Vx(i*m_k) -= velMx/M;
    m_Vy(i*m_k) -= velMy/M;
    m_Vz(i*m_k) -= velMz/M;
  }
};

//Method that sets up an initializer with the Sun in the center, returns nothing.
void PlanetSolver::init_sun_center(vector<string> names, double beta, int N, int k, double T, bool test_convergence = false){
  initialize(beta, N, k, T);    //Method that sets up initializer for the solver, returns nothing.

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();    //Reads in position and velocity from NASA data.
  vec params = vec(7);

  m_masses = zeros<vec>(m_N);

  double posx0, posy0, posz0, velx0, vely0, velz0;
  if (test_convergence == false){
    for (int i = 0; i < m_N; i++){
      params = Planet.initialize(m_names[i]);
      m_masses(i) = params(0);

      //Getting values from terminal.
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

        //Updating the values from terminal input.
        m_X(i*m_k) = posx0; m_Y(i*m_k) = posy0; m_Z(i*m_k) = posz0;
        m_Vx(i*m_k) = velx0; m_Vy(i*m_k) = vely0; m_Vz(i*m_k) = velz0;
        m_ax(i*m_k) = m_ay(i) = m_az(i*m_k) = 0.0;
        }
      }
  } else {
    for (int i = 0; i < m_N; i++){
      params = Planet.initialize(m_names[i]);   //Creating planet objects.
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

//Method that solves the problem, using different methods taken in as parameter, returns nothing.
void PlanetSolver::solvesystem(bool check, int method){
 int s;
 if(check == true) {
   s = 1;
 } else {
   s = 0;
 }

  if (method == 1){ //use velocity verlet
    auto start = chrono::high_resolution_clock::now(); //Start timing velocity verlet method on the whole system.
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        verlet_pos(i,j);
      }
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j+1);
        verlet_vel(i,j);
        kin(i*m_k + j+1) = kinetic_energy(i, j);
        tot(i*m_k + j+1) = kin(i*m_k + j+1) + pot(i*m_k + j+1); //get total energy E = K + V for each planet
      }
    }
    auto finish = chrono::high_resolution_clock::now(); //End timer
    auto cpu_time_verlet = duration_cast<milliseconds>(finish-start).count();
    cout << "Time it took to compute the verlet method for the system: " << cpu_time_verlet << "ms"<< endl;
  }
  if (method == 2){ //use velocity verlet
    auto start = chrono::high_resolution_clock::now(); //Start timing forward euler method on the whole system.
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j);
        forwardeuler(i,j);
        kin(i*m_k + j+1) = kinetic_energy(i, j);
        tot(i*m_k + j+1) = kin(i*m_k + j+1) + pot(i*m_k + j+1); //get total energy E = K + V for each planet
        }
      }
      auto finish = chrono::high_resolution_clock::now(); //End timer
      auto cpu_time_fe = duration_cast<milliseconds>(finish-start).count();
      cout << "Time it took to compute the Forward-Euler method for the system: " << cpu_time_fe << "ms"<< endl;
    }

  if (method == 3){ //Use Euler-Cromer
    auto start = chrono::high_resolution_clock::now(); //Start timing euler-cromer method on the whole system.
    for (int j = 0; j < m_k-1; j++){ // for time
      for (int i = s; i < m_N; i++){ //for planets
        force_a(i,j);
        eulerchromer(i,j);
        kin(i*m_k + j+1) = kinetic_energy(i, j);
        tot(i*m_k + j+1) = kin(i*m_k + j+1) + pot(i*m_k + j+1); //get total energy E = K + V for each planet
      }
    }
    auto finish = chrono::high_resolution_clock::now(); //End timer
    auto cpu_time_ec = duration_cast<milliseconds>(finish-start).count();
    cout << "Time it took to compute the Euler-Cromer method for the system: " << cpu_time_ec << "ms"<< endl;
  }
}

//Method that gives a random, uniformly distributed index, returns the index.
int PlanetSolver::random_index_generator(int min, int max){
  // Using random generator from namespace std
  random_device seed;
  mt19937 rng(seed());
  uniform_int_distribution<int> uni(min,max);
  int random_integer = uni(rng);
  return random_integer;
}

//Method that tests if the energy is conserved, returns nothing.
void PlanetSolver::test_constant_energy(double tol){
  int j = random_index_generator(0,m_k);    //Picking out random index.
  for (int i = 1; i < m_N; i++){
    if (tot(i*m_k + 1) - tot(i*m_k + j) < tol) { //Test against tolerance.
      continue;
    } else {
    cout << "Error: Energy not conserved for celestial bodies with tolerance:" << " " << tol << "\n";
    break;
    }
  }
}

//Method that tests if the angular momentum is constant, returns nothing.
void PlanetSolver::test_constant_angular(double tol){
  double diffL;
  // get angular momentum for all times
  get_angular_Momentum();
  int j = random_index_generator(0,m_k); //Random index.
  for (int i = 1; i < m_N; i++){
    diffL = (m_Lx(i*m_k + 1)*m_Lx(i*m_k + 1) + m_Ly(i*m_k + 1)*m_Ly(i*m_k + 1))  \
        - (m_Lx(i*m_k + j)*m_Lx(i*m_k + j) + m_Ly(i*m_k + j)*m_Ly(i*m_k + j));    //Computing the difference in angular momentum at different time steps.
    if (diffL < tol) {
      continue;
    } else {
    cout << "Error: Angular momentum not conserved for celestial bodies with tolerance:" << " " << tol << "\n";
    break;
    }
  }
}

//Method that tests if an orbit is cicular or not, returns nothing.
void PlanetSolver::test_circular_orbit(double tol){
  int i = 1;
  int j = random_index_generator(0,m_k);      //Random index method.
  double diffr = abs(1 - sqrt(m_X(i*m_k+j)*m_X(i*m_k+j) \
              + m_Y(i*m_k+j)*m_Y(i*m_k+j) + m_Z(i*m_k+j)*m_Z(i*m_k+j)));    //Checking the difference in distance.
  if (diffr > tol){
    cout << "Orbit is not circular with tolerance:" << " " << tol << "\n";
    }
}

//Method that tests the convergence, returns nothing.
void PlanetSolver::test_convergence(vector<string> names,double beta, int N,int k, double T, int N_experiments, int method){
  // must solve for several dt's and check stability
  // --> need to call the solver in a loop

  //Initializing arrays.
  vec E = zeros<vec>(N_experiments);
  vec rel_error = zeros<vec>(N_experiments);
  vec r = zeros<vec>(N_experiments-1);

  bool sun_center = true;   //Set Sun as center of system.
  int steps = 0;
  double next_error;
  int i = 1;
  double factor = 2;
  bool test_convergence = true;
  while (steps < N_experiments){
    double error = 0;
    init_sun_center(names,beta,N,k,T,test_convergence);   //Initializing a system with Sun as center and testing convergence.
    solvesystem(sun_center,method);
    // Calculate error and convergence rate
    for (int j = 0; j < m_k; j++){    //For planets
      next_error =  1 - sqrt(m_X(i*m_k+j)*m_X(i*m_k+j) \
      + m_Y(i*m_k+j)*m_Y(i*m_k+j) + m_Z(i*m_k+j)*m_Z(i*m_k+j));   //Calulating the error.
      error = error + next_error*next_error;
      }

    rel_error(steps) = next_error; //last point
    E(steps) = sqrt(m_h*error);
    k = m_k*factor;
    steps += 1;
    }
    for (int j = 1; j < N_experiments; j++){
      r(j-1) = log(E(j)/E(j-1))/log(1/factor);
    }
    cout << "Convergence rates:" << " " << r <<"\n";
    cout << "relative error:" << " " << rel_error << "\n";
    // get relative error for last step
}

//Method that writes position to file, returns nothing.
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

//Method that writes velocity to file, returns nothing.
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
