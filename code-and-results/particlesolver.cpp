#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(double beta, int N, int k, double T){
  m_N = N;
  m_k = k;
  m_T = T;
  m_T0 = 0.0;
  m_h = (m_T-m_T0)/(m_k - 1); // time step
  //cout << m_h << "\n";
  m_beta = beta;
  hh = m_h*m_h;

  //initialize vectors
  m_X = zeros<vec>(m_N*m_k);  // N number of planets, m_k number of time steps
  m_Y = zeros<vec>(m_N*m_k);
  m_Z = zeros<vec>(m_N*m_k);

  m_Vx = zeros<vec>(m_N*m_k);
  m_Vy = zeros<vec>(m_N*m_k);
  m_Vz = zeros<vec>(m_N*m_k);

  m_ax = zeros<vec>(m_N*m_k);
  m_ay = zeros<vec>(m_N*m_k);
  m_az = zeros<vec>(m_N*m_k);

  m_Etot = zeros<vec>(m_N*m_k);
  pot = zeros<vec>(m_N*m_k);
  kin = zeros<vec>(m_N*m_k);
  tot = zeros<vec>(m_N*m_k);
  m_Lx = zeros<vec>(m_N*m_k);
  m_Ly = zeros<vec>(m_N*m_k);
  m_Lz = zeros<vec>(m_N*m_k);

  //for (int i = 0; i < ) fill initial conditions for position, velcocity and acceleration

};

void ParticleSolver::force_a(int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double a_x = 0;
  double a_y = 0;
  double a_z = 0;
  double V = 0;
  double diffx,diffy,diffz,diffr,r;
  for (int i = 0; i < m_N; i++){ //for planets
    if (i != l){                 //l is the index of the specific planet we are looking at
      diffx = m_X(l*m_k+j) - m_X(i*m_k+j);   //j is the given timestep put planet distance in own method?
      diffy = m_Y(l*m_k+j) - m_Y(i*m_k+j);
      diffz = m_Z(l*m_k+j) - m_Z(i*m_k+j);
      diffr = diffx*diffx + diffy*diffy + diffz*diffz;
      r = pow(diffr,0.5);
      double r_term = pow(r,(m_beta+1));

      // calculate potential energies for a planet
      V += potential_energy(r,l,i,j);
      //calculate gravitational acceleration
      a_x += ((m_X(l*m_k+j)-m_X(i*m_k+j))*G*m_masses(i))/r_term;
      a_y += ((m_Y(l*m_k+j)-m_Y(i*m_k+j))*G*m_masses(i))/r_term;
      a_z += ((m_Z(l*m_k+j)-m_Z(i*m_k+j))*G*m_masses(i))/r_term;

    }
  }
  m_ax(l*m_k+j) = -a_x;
  m_ay(l*m_k+j) = -a_y;
  m_az(l*m_k+j) = -a_z;

  pot(l*m_k + j) = V;
  kin(l*m_k + j) = kinetic_energy(l, j);
  tot(l*m_k + j) = kin(l*m_k + j) + pot(l*m_k + j); //get total energy E = K + V for each planet
}


void ParticleSolver::verlet_pos(int l, int j){
  m_X(l*m_k+j+1) = m_X(l*m_k+j) + m_h*m_Vx(l*m_k+j) + (1./2)*hh*m_ax(l*m_k+j);
  m_Y(l*m_k+j+1) = m_Y(l*m_k+j) + m_h*m_Vy(l*m_k+j) + (1./2)*hh*m_ay(l*m_k+j);
  m_Z(l*m_k+j+1) = m_Z(l*m_k+j) + m_h*m_Vz(l*m_k+j) + (1./2)*hh*m_az(l*m_k+j);
}

void ParticleSolver::verlet_vel(int l, int j){
  m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + (1./2)*m_h*(m_ax(l*m_k+j+1) + m_ax(l*m_k+j));
  m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + (1./2)*m_h*(m_ay(l*m_k+j+1) + m_ay(l*m_k+j));
  m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + (1./2)*m_h*(m_az(l*m_k+j+1) + m_az(l*m_k+j));
}


void ParticleSolver::eulerchromer(int l, int j){
  m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + m_h*m_ax(l*m_k+j);
  m_X(l*m_k+j+1) = m_X(l*m_k+j) +   m_h*m_Vx(l*m_k+j+1);

  m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + m_h*m_ay(l*m_k+j);
  m_Y(l*m_k+j+1) = m_Y(l*m_k+j) +   m_h*m_Vy(l*m_k+j+1);

  m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + m_h*m_az(l*m_k+j);
  m_Z(l*m_k+j+1) = m_Z(l*m_k+j) +   m_h*m_Vz(l*m_k+j+1);
};

void ParticleSolver::forwardeuler(int l, int j){
    m_Vx(l*m_k+j+1) = m_Vx(l*m_k+j) + m_h*m_ax(l*m_k+j);
    m_X(l*m_k+j+1) = m_X(l*m_k+j) +   m_h*m_Vx(l*m_k+j);

    m_Vy(l*m_k+j+1) = m_Vy(l*m_k+j) + m_h*m_ay(l*m_k+j);
    m_Y(l*m_k+j+1) = m_Y(l*m_k+j) +   m_h*m_Vy(l*m_k+j);

    m_Vz(l*m_k+j+1) = m_Vz(l*m_k+j) + m_h*m_az(l*m_k+j);
    m_Z(l*m_k+j+1) = m_Z(l*m_k+j) +   m_h*m_Vz(l*m_k+j);
};


double ParticleSolver::kinetic_energy(int l, int j){ // kinetic energy for one planet
  double vnorm2 = m_Vx(l*m_k+j)*m_Vx(l*m_k + j) + m_Vy(l*m_k+j)*m_Vy(l*m_k + j)  +m_Vy(l*m_k+j)*m_Vy(l*m_k + j);
  double kinetic = 0.5*m_masses(l)*vnorm2;
  return kinetic;

}


double ParticleSolver::potential_energy(double r, int l, int i, int j){ // potential energy between two objects /general central potential
  double potential = 4*(1/(-m_beta + 1))*M_PI*m_masses(l)*m_masses(i)/pow(r,m_beta-1); // check sign!!
  return potential;
}

double ParticleSolver::angular_momentum(double pos1, double v1, double pos2, double v2){
  return pos1*v2 - pos2*v1;
}

void ParticleSolver::get_angular_momentum(int j){
  for (int i = 0; i < m_N; i++){ //for planets
    m_Lx(i*m_k+j) = angular_momentum(m_Y(i*m_k+j), m_Vy(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Ly(i*m_k+j) = angular_momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Lz(i*m_k+j) = angular_momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Y(i*m_k+j), m_Vy(i*m_k+j));
  }
}

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

void ParticleSolver::write_angular_momentum_to_file() {
  ofstream x;
  ofstream y;
  ofstream z;

  string filename_1("./results/ang_x.txt");
  string filename_2("./results/ang_y.txt");
  string filename_3("./results/ang_z.txt");

  x.open(filename_1);
  y.open(filename_2);
  z.open(filename_3);

  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      x << m_Lx(i*m_k) << " ";
      y << m_Ly(i*m_k) << " ";
      z << m_Lz(i*m_k) << " ";
    }
    x << "\n";
    y << "\n";
    z << "\n";
  }
  x.close();
  y.close();
  z.close();
}
