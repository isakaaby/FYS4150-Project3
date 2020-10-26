#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;



//Method that sets up initializer for the solver, returns nothing.
void ParticleSolver::initialize(double beta, int N, int k, double T){
  m_N = N;      // N number of planets
  m_k = k;      //m_k number of time steps
  m_T = T;
  m_T0 = 0.0;
  m_h = (m_T-m_T0)/(m_k - 1); // time step
  //cout << m_h << "\n";
  m_beta = beta;
  hh = m_h*m_h;

  //initialize position vectors
  m_X = zeros<vec>(m_N*m_k);
  m_Y = zeros<vec>(m_N*m_k);
  m_Z = zeros<vec>(m_N*m_k);

  //initialize velocity vectors
  m_Vx = zeros<vec>(m_N*m_k);
  m_Vy = zeros<vec>(m_N*m_k);
  m_Vz = zeros<vec>(m_N*m_k);

  //initialize acceleration vectors
  m_ax = zeros<vec>(m_N*m_k);
  m_ay = zeros<vec>(m_N*m_k);
  m_az = zeros<vec>(m_N*m_k);

  //initialize energy vectors
  m_Etot = zeros<vec>(m_N*m_k);
  pot = zeros<vec>(m_N*m_k);
  kin = zeros<vec>(m_N*m_k);
  tot = zeros<vec>(m_N*m_k);

  //initialize angular momentum vectors
  m_Lx = zeros<vec>(m_N*m_k);
  m_Ly = zeros<vec>(m_N*m_k);
  m_Lz = zeros<vec>(m_N*m_k);


};

//Method that computes the gravitational force, returns nothing.
void ParticleSolver::force_a(int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double a_x = 0;
  double a_y = 0;
  double a_z = 0;
  double V = 0;
  double diffx,diffy,diffz,diffr,r;
  for (int i = 0; i < m_N; i++){ //for planets
    if (i != l){                 //l is the index of the specific planet we are looking at

      //Calculating the differences.
      diffx = m_X(l*m_k+j) - m_X(i*m_k+j);
      diffy = m_Y(l*m_k+j) - m_Y(i*m_k+j);
      diffz = m_Z(l*m_k+j) - m_Z(i*m_k+j);
      diffr = diffx*diffx + diffy*diffy + diffz*diffz;
      r = pow(diffr,0.5);
      double r_term = pow(r,(m_beta+1));

      // calculate potential energy for a planet
      V += potential_energy(r,l,i,j);

      //calculate gravitational acceleration
      a_x += ((m_X(l*m_k+j)-m_X(i*m_k+j))*G*m_masses(i))/r_term;
      a_y += ((m_Y(l*m_k+j)-m_Y(i*m_k+j))*G*m_masses(i))/r_term;
      a_z += ((m_Z(l*m_k+j)-m_Z(i*m_k+j))*G*m_masses(i))/r_term;

    }
  }

  //Updating the accelaration and energies.
  m_ax(l*m_k+j) = -a_x;
  m_ay(l*m_k+j) = -a_y;
  m_az(l*m_k+j) = -a_z;

  pot(l*m_k + j) = V;
  kin(l*m_k + j) = kinetic_energy(l, j);
  tot(l*m_k + j) = kin(l*m_k + j) + pot(l*m_k + j); //get total energy E = K + V for each planet
}

//Method that computes the positions using the Verlet method, returns nothing.
void ParticleSolver::verlet_pos(int l, int j){
  m_X(l*m_k+j+1) = m_X(l*m_k+j) + m_h*m_Vx(l*m_k+j) + (1./2)*hh*m_ax(l*m_k+j);
  m_Y(l*m_k+j+1) = m_Y(l*m_k+j) + m_h*m_Vy(l*m_k+j) + (1./2)*hh*m_ay(l*m_k+j);
  m_Z(l*m_k+j+1) = m_Z(l*m_k+j) + m_h*m_Vz(l*m_k+j) + (1./2)*hh*m_az(l*m_k+j);
}

//Method that computes the velocities using the Verlet method, returns nothing.
void ParticleSolver::verlet_vel(int l, int j){
  m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + (1./2)*m_h*(m_ax(l*m_k+j+1) + m_ax(l*m_k+j));
  m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + (1./2)*m_h*(m_ay(l*m_k+j+1) + m_ay(l*m_k+j));
  m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + (1./2)*m_h*(m_az(l*m_k+j+1) + m_az(l*m_k+j));
}

//Method that calculates the position, velocity and acceleration using the Eueler-Cromer method, returns nothing.
void ParticleSolver::eulerchromer(int l, int j){
  m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + m_h*m_ax(l*m_k+j);
  m_X(l*m_k+j+1) = m_X(l*m_k+j) +   m_h*m_Vx(l*m_k+j+1);

  m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + m_h*m_ay(l*m_k+j);
  m_Y(l*m_k+j+1) = m_Y(l*m_k+j) +   m_h*m_Vy(l*m_k+j+1);

  m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + m_h*m_az(l*m_k+j);
  m_Z(l*m_k+j+1) = m_Z(l*m_k+j) +   m_h*m_Vz(l*m_k+j+1);
};

//Method that calculates the position, velocity and acceleration using the Forward Eurler method, returns nothing.
void ParticleSolver::forwardeuler(int l, int j){
    m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + m_h*m_ax(l*m_k+j);
    m_X(l*m_k+j+1) = m_X(l*m_k+j) +   m_h*m_Vx(l*m_k+j);

    m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + m_h*m_ay(l*m_k+j);
    m_Y(l*m_k+j+1) = m_Y(l*m_k+j) +   m_h*m_Vy(l*m_k+j);

    m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + m_h*m_az(l*m_k+j);
    m_Z(l*m_k+j+1) = m_Z(l*m_k+j) +   m_h*m_Vz(l*m_k+j);
};

//Method that calculates the kinetic energy, returns the energy.
double ParticleSolver::kinetic_energy(int l, int j){ // kinetic energy for one planet
  double vnorm2 = m_Vx(l*m_k+j)*m_Vx(l*m_k + j) + m_Vy(l*m_k+j)*m_Vy(l*m_k + j)  +m_Vy(l*m_k+j)*m_Vy(l*m_k + j);
  double kinetic = 0.5*m_masses(l)*vnorm2;
  return kinetic;

}

//Method that calculates the potential energy, returns the energy.
double ParticleSolver::potential_energy(double r, int l, int i, int j){ // potential energy between two objects /general central potential
  double potential = 4*(1/(-m_beta + 1))*M_PI*m_masses(l)*m_masses(i)/pow(r,m_beta-1); // check sign!!
  return potential;
}

//Method that calculates the angular momentum, returns the momentum.
double ParticleSolver::angular_Momentum(double pos1, double v1, double pos2, double v2){
  return pos1*v2 - pos2*v1;
}

//Method that updates the value for the angular momentum in three directions, returns nothing.
void ParticleSolver::get_angular_Momentum(){
  for (int i = 0; i < m_N; i++){ //for planets
    for(int j = 0; j < m_k; j++) {
    m_Lx(i*m_k+j) = angular_Momentum(m_Y(i*m_k+j), m_Vy(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Ly(i*m_k+j) = angular_Momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Lz(i*m_k+j) = angular_Momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Y(i*m_k+j), m_Vy(i*m_k+j));
    }
  }
}

//Method that writes energy to file, returns nothing.
void ParticleSolver::write_energy_to_file() {
  ofstream p;
  ofstream k;
  ofstream t;

  string filename_1("./results/kinetic.txt");
  string filename_2("./results/potential.txt");
  string filename_3("./results/total.txt");

  p.open(filename_1);
  k.open(filename_2);
  t.open(filename_3);

  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      p << pot(i*m_k+j) << " ";
      k << kin(i*m_k+j) << " ";
      t << tot(i*m_k+j) << " ";
    }
    p << "\n";
    k << "\n";
    t << "\n";
  }
  p.close();
  k.close();
  t.close();
}
