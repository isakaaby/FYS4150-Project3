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

  //initialize vectors
  m_k = k;
  vec m_X = zeros<vec>(N*m_k);  // N number of planets, m_k number of time step
  vec m_Y = zeros<vec>(N*m_k);

  vec m_Vx = zeros<vec>(N);
  vec m_Vy = zeros<vec>(N);
  vec m_Vx_prev = zeros<vec>(N);
  vec m_Vy_prev = zeros<vec>(N);

  vec m_ax zeros<vec>(N);
  vec m_ax_prev zeros<vec>(N);
  vec m_ay zeros<vec>(N);
  vec m_ay_prev zeros<vec>(N);
};

void ParticleSolver::Verlet(double force(double x,y){
  h = (m_T-T0)/(m_k - 1); // time step
  for (int j = 0; int j < m_k){ // for time
    for (int i = 0; int i < m_N){ //for planets
      m_X(i*m_k+j+1) = m_X(i*m_k+j) + h*m_Vx_prev(i) + (1./2)*h*h*m_ax_prev(i);
      m_ax(i) = force(m_X(i*m_k+j+1),m_Y(i*m_k+j+1));
      m_Vx(i) = m_Vx(i) + (1./2)*h*h*(m_ax(i) + m_ax_prev(i));

      m_Y(i*m_k+j+1) = m_Y(i*m_k+j) + h*m_Vy_prev(i) + (1./2)*h*h*m_ay_prev(i);
      m_ay(i) = force(m_X(i*m_k+j+1),m_Y(i*m_k+j+1));
      m_Vy(i) = m_Vy(i) + (1./2)*h*h*(m_ay(i) + m_ay_prev(i));

      // Update to prev
      m_ax_prev(i) = m_ax(i);
      m_ay_prev(i) = m_ay(i);
      m_Vx_prev(i) = m_Vx(i)
      m_Vy_prev(i) = m_Vy(i)
    }
  }
};

/*void ParticleSolver::EulerChromer(){

};

double get_force() {


}

void ParticleSolver::RungeKutta4(){
  Write Runge Kutta here
};*/
