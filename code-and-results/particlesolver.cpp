#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(double beta, int N, int k, int T){

  m_N = N;
  m_k = k;
  m_T = T;
  m_T0 = 0.0;
  m_h = (m_T-m_T0)/(m_k - 1); // time step
  m_beta = beta;

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

  //for (int i = 0; i < ) fill initial conditions for position, velcocity and acceleration
};

double ParticleSolver::force_a(vec pos, int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double a = 0;
  double diffx,diffy,diffz,diffr,r;
  for (int i = 0; i < m_N; i++){ //for planets
    if (i != l){                 //l is the index of the specific planet we are looking at
      diffx = m_X(l*m_k+j) - m_X(i*m_k+j);   //j is the given timestep
      diffy = m_Y(l*m_k+j) - m_Y(i*m_k+j);
      diffz = m_Z(l*m_k+j) - m_Z(i*m_k+j);
      diffr = diffx*diffx + diffy*diffy + diffz*diffz;
      r = pow(diffr,(m_beta+1)/2);
      a += ((pos(l*m_k+j)-pos(i*m_k+j))*G*m_masses(i))/r;
    }
  }
  return -a;
}

void ParticleSolver::verlet(){
  double h = m_h;
  for (int j = 0; j < m_k-1; j++){ // for time
    for (int i = 0; i < m_N; i++){ //for planets
      m_X(i*m_k+j+1) = m_X(i*m_k+j) + h*m_Vx(i*m_k+j) + (1./2.)*h*h*m_ax(i*m_k+j);
      m_Y(i*m_k+j+1) = m_Y(i*m_k+j) + h*m_Vy(i*m_k+j) + (1./2.)*h*h*m_ay(i*m_k+j);
      m_Z(i*m_k+j+1) = m_Z(i*m_k+j) + h*m_Vz(i*m_k+j) + (1./2.)*h*h*m_az(i*m_k+j);
      //cout << m_X(i*m_k+j+1) << "\n";
      //cout << m_Y(i*m_k+j+1) << "\n";
      //cout << m_Z(i*m_k+j+1) << "\n";
    }
    for (int i = 0; i < m_N; i++){
      m_ax(i*m_k+j+1) = force_a(m_X,i,j+1);
      m_Vx(i*m_k+j+1) = m_Vx(i*m_k+j) + (1./2)*h*(m_ax(i*m_k+j+1) + m_ax(i*m_k+j));

      m_ay(i*m_k+j+1) = force_a(m_Y,i,j+1);
      m_Vy(i*m_k+j+1) = m_Vy(i*m_k+j) + (1./2)*h*(m_ay(i*m_k+j+1) + m_ay(i*m_k+j));

      m_az(i*m_k+j+1) = force_a(m_Z,i,j+1);
      m_Vz(i*m_k+j+1) = m_Vz(i*m_k+j) + (1./2)*h*(m_az(i*m_k+j+1) + m_az(i*m_k+j));
    }
  }
}


void ParticleSolver::eulerchromer(){
  double h = m_h;
  for (int j = 0; j < m_k-1; j++){ // for time
    for (int i = 0; i < m_N; i++){ //for planets
    m_ax(i*m_k+j) = force_a(m_X,i,j);
    m_Vx(i*m_k+j+1) = m_Vx(i*m_k+j) + h*m_ax(i*m_k+j);
    m_X(i*m_k+j+1) = m_X(i*m_k+j) +   h*m_Vx(i*m_k+j+1);

    m_ay(i*m_k+j) = force_a(m_Y,i,j);
    m_Vy(i*m_k+j+1) = m_Vy(i*m_k+j) + h*m_ay(i*m_k+j);
    m_Y(i*m_k+j+1) = m_Y(i*m_k+j) +   h*m_Vy(i*m_k+j+1);

    m_az(i*m_k+j) = force_a(m_Y,i,j);
    m_Vz(i*m_k+j+1) = m_Vz(i*m_k+j) + h*m_az(i*m_k+j);
    m_Z(i*m_k+j+1) = m_Z(i*m_k+j) +   h*m_Vz(i*m_k+j+1);
    };
  };
};


/*void ParticleSolver::RK4_xupdate(double t, double x,double y, double v, double f1(double t, double x, double y, double v), double f2(double t, double x, double y, double v)){
  // s is spatial variable, v is velocity associated with s
  double h = m_h;
  double K1s = f1(v);
  double K1v = f2(x,y);
  double K2s = f1(v + h*K1v/2);
  double K2v = f2(x + h*K1s/2,y);
  double K3s = f1(v + h*K2v/2);
  double K3v = f2(x + h*K2s/2,y);
  double K4s = f1(v + h*K3v/2);
  double K4v = f2(x + h*K3x/2,y);

  double s_new = x + h/6*(K1s + 2*K2s + 2*K3s + K4s);
  double v_new = v + h/6*(K1v + 2*K2v + 2*K3v + K4v);

  vec m_xupdate  = zeros<vec>(2);
  m_xupdate(0) = s_new ; m_xupdate(1) = v_new;
};

void ParticleSolver::RK4_yupdate(double t,double x,double y, double v, double f1(double t, double x, double y, double v), double f2(double t, double x, double y, double v)){
  // s is spatial variable, v is velocity associated with s
  double h = m_h;

  double h = m_h;
  double K1s = f1(v);
  double K1v = f2(x,y);
  double K2s = f1(v + h*K1v/2);
  double K2v = f2(x,y + h*K1s/2);
  double K3s = f1(v + h*K2v/2);
  double K3v = f2(x,y + h*K2s/2);
  double K4s = f1(v + h*K3v/2);
  double K4v = f2(x,y + h*K3x/2);

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
*/
