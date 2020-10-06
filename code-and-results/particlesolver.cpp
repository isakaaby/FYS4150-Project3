#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(int N, int k){
  m_N = N;

  //to be used in all derived classes
  vec m = vec(N);
  Planets Planet;
  m = Planet.initialize("Earth");

  /*m(0) = Planet.get_Mass("earth");
  m(1) = Planet.get_Mass("jupiter");
  m(2) = Planet.get_Mass("mars");
  m(3) = Planet.get_Mass("venus");
  m(4) = Planet.get_Mass("saturn");
  m(5) = Planet.get_Mass("mercury");
  m(6) = Planet.get_Mass("uranus");
  m(7) = Planet.get_Mass("neptune");
  m(8) = Planet.get_Mass("pluto");

  vec r = vec(N);
  r(0) = Planet.get_Distance("earth");
  r(1) = Planet.get_Distance("jupiter");
  r(2) = Planet.get_Distance("mars");
  r(3) = Planet.get_Distance("venus");
  r(4) = Planet.get_Distance("saturn");
  r(5) = Planet.get_Distance("mercury");
  r(6) = Planet.get_Distance("uranus");
  r(7) = Planet.get_Distance("neptune");
  r(8) = Planet.get_Distance("pluto");*/

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

double get_force() {


}

void ParticleSolver::RungeKutta4(){

};*/
