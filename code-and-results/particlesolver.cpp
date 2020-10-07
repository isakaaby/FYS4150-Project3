#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(int N, int k, int m_T){

  m_N = N;
  m_k = k;
  double step = (m_T-m_T0)/(m_k - 1); // time step
  m_h = step;

  //initialize vectors
  vec m_X = zeros<vec>(N*m_k);  // N number of planets, m_k number of time step
  vec m_Y = zeros<vec>(N*m_k);

  vec m_Vx = zeros<vec>(N);
  vec m_Vy = zeros<vec>(N);
  vec m_Vx_prev = zeros<vec>(N);
  vec m_Vy_prev = zeros<vec>(N);

  vec m_ax = zeros<vec>(N);
  vec m_ax_prev = zeros<vec>(N);
  vec m_ay = zeros<vec>(N);
  vec m_ay_prev = zeros<vec>(N);
};

void ParticleSolver::verlet(double force(double x, double y)){
  double h = m_h;
  for (int j = 0; j < m_k; ++j){ // for time
    for (int i = 0; i < m_N; ++i){ //for planets
      m_X(i*m_k+j+1) = m_X(i*m_k+j) + h*m_Vx_prev(i) + (1./2)*h*h*m_ax_prev(i);
      m_ax(i) = force(m_X(i*m_k+j+1),m_Y(i*m_k+j+1));
      m_Vx(i) = m_Vx(i) + (1./2)*h*h*(m_ax(i) + m_ax_prev(i));

      m_Y(i*m_k+j+1) = m_Y(i*m_k+j) + h*m_Vy_prev(i) + (1./2)*h*h*m_ay_prev(i);
      m_ay(i) = force(m_X(i*m_k+j+1),m_Y(i*m_k+j+1));
      m_Vy(i) = m_Vy(i) + (1./2)*h*h*(m_ay(i) + m_ay_prev(i));

      // Update to prev
      m_ax_prev(i) = m_ax(i);
      m_ay_prev(i) = m_ay(i);
      m_Vx_prev(i) = m_Vx(i);
      m_Vy_prev(i) = m_Vy(i);
    };
  };
};

void ParticleSolver::RK4_xupdate(double t, double x,double y, double v, double f1(double t, double x, double y, double v), double f2(double t, double x, double y, double v)){
  // s is spatial variable, v is velocity associated with s
  double h = m_h;

  double K1s = f1(t,x,y,v);
  double K1v = f2(t,x,y,v);
  double K2s = f1(t + h/2,x + h*K1s/2,y,v + h*K1v/2);
  double K2v = f2(t + h/2,x + h*K1s/2,y,v + h*K1v/2);
  double K3s = f1(t + h/2,x + h*K2s/2,y,v + h*K2v/2);
  double K3v = f2(t + h/2,x + h*K2s/2,y,v + h*K2v/2);
  double K4s = f1(t + h/2,x + h*K3s/2,y,v + h*K3v/2);
  double K4v = f2(t + h/2,x + h*K3s/2,y,v + h*K3v/2);

  double s_new = x + h/6*(K1s + 2*K2s + 2*K3s + K4s);
  double v_new = v + h/6*(K1v + 2*K2v + 2*K3v + K4v);

  vec m_xupdate  = zeros<vec>(2);
  m_xupdate(0) = s_new ; m_xupdate(1) = v_new;
};

void ParticleSolver::RK4_yupdate(double t,double x,double y, double v, double f1(double t, double x, double y, double v), double f2(double t, double x, double y, double v)){
  // s is spatial variable, v is velocity associated with s
  double h = m_h;

  double K1s = f1(t,x,y,v);
  double K1v = f2(t,x,y,v);
  double K2s = f1(t + h/2,x,y + h*K1s/2,v + h*K1v/2);
  double K2v = f2(t + h/2,x,y + h*K1s/2,v + h*K1v/2);
  double K3s = f1(t + h/2,x,y + h*K2s/2,v + h*K2v/2);
  double K3v = f2(t + h/2,x,y + h*K2s/2,v + h*K2v/2);
  double K4s = f1(t + h/2,x,y + h*K3s/2,v + h*K3v/2);
  double K4v = f2(t + h/2,x,y + h*K3s/2,v + h*K3v/2);

  double s_new = x + h/6*(K1s + 2*K2s + 2*K3s + K4s);
  double v_new = v + h/6*(K1v + 2*K2v + 2*K3v + K4v);

  vec m_yupdate  = zeros<vec>(2);
  m_yupdate(0) = s_new ; m_yupdate(1) = v_new;
};

void ParticleSolver::RK4(double f1(double t, double x, double y, double v),double force(double t, double x, double y, double v)){
  //Write Runge Kutta here
  for (int j = 0; j < m_k; ++j){ // for time
    for (int i = 0; i < m_N; ++i){ //for planets
      //RK4_y_consts()
      //RK4_x_consts()
    };
  };
};

/* void ParticleSolver::EulerChromer(){
};

double get_force() {
}
*/
