#include "particlesolver.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(int N, int k){
  m_N = N;
  //to be used in all derived classes

  //initialize vector
  m_k = k;
  vec X = zeros<vec>(N*m_k);
  vec Y = zeros<vec>(N*m_k);
  vec Z = zeros<vec>(N*m_k);

};

/*void ParticleSolver::Verlet(double force(double x, double y, double z)){

};

void ParticleSolver::EulerChromer(){

};

void ParticleSolver::RungeKutta4(){

};*/
