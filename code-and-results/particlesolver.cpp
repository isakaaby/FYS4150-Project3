#include "particlesolver.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

mat ParticleSolver::initialize(int N, int k,){
  m_N = N;
  //to be used in all derived classes

  //initialize vectors
  vec X = zeros<vec>(N*K)
  vec Y = zeros<vec>(N*K)
  vec Z = zeros<vec>(N*K)

};

void ParticleSolver::Verlet(double force(double x, double y, double z)){

};

void ParticleSolver::EulerChromer(){

};

void ParticleSolver::RungeKutta4(){

};
