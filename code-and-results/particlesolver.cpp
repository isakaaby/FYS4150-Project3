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
  m_Lx = zeros<vec>(m_N*m_k);
  m_Ly = zeros<vec>(m_N*m_k);
  m_Lz = zeros<vec>(m_N*m_k);

  //for (int i = 0; i < ) fill initial conditions for position, velcocity and acceleration

};

double ParticleSolver::force_a(vec pos, int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double a = 0;
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
      a += ((pos(l*m_k+j)-pos(i*m_k+j))*G*m_masses(i))/r_term;
    }
  }
  m_Etot(l*m_k+j) = kinetic_energy(l,j) + V; // get total energy E = K + V for each planet
  return -a;
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
