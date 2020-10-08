#include "particlesolver.hpp"
#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;


//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

void ParticleSolver::initialize(double beta, int N, int k, int m_T){

  m_N = N;
  m_k = k;
  double step = (m_T-m_T0)/(m_k - 1); // time step
  m_h = step;
  m_beta = beta;

  //initialize vectors
  vec m_X = zeros<vec>(N*m_k);  // N number of planets, m_k number of time step
  vec m_Y = zeros<vec>(N*m_k);
  vec m_Z = zeros<vec>(N*m_k);

  vec m_Vx = zeros<vec>(N*m_k);
  vec m_Vy = zeros<vec>(N*m_k);
  vec m_Vz = zeros<vec>(N*m_k);

  vec m_ax = zeros<vec>(N*m_k);
  vec m_ay = zeros<vec>(N*m_k);
  vec m_az = zeros<vec>(N*m_k);

};

double ParticleSolver::force_a(vec pos, double x, double y, double z, int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double a = 0;
  double diffx,diffy,diffz,diffr,r;
  for (int i = 0; i < m_N; ++i){ //for planets
    if (i != l){                 //l is the index of the specific planet we are looking at
      diffx = x - m_X(i*m_k+j);   //j is the given timestep
      diffy = y - m_Y(i*m_k+j);
      diffz = z - m_Z(i*m_k+j);
      diffr = diffx*diffx + diffy*diffy + diffz*diffz;
      r = pow(diffr,(m_beta+1)/2);
      a += ((pos(l*m_k+j+1)-pos(i*m_k+j))*G*m_masses(i))/r;
    }
  }
  return a;
}

void ParticleSolver::verlet(double force(double s, double x, double y)){
  double h = m_h;
  for (int j = 0; j < m_k; ++j){ // for time
    for (int i = 0; i < m_N; ++i){ //for planets
      m_X(i*m_k+j+1) = m_X(i*m_k+j) + h*m_Vx(i*m_k+j) + (1./2)*h*h*m_ax(i*m_k+j);
      m_Y(i*m_k+j+1) = m_Y(i*m_k+j) + h*m_Vy(i*m_k+j) + (1./2)*h*h*m_ay(i*m_k+j);
      m_Z(i*m_k+j+1) = m_Z(i*m_k+j) + h*m_Vz(i*m_k+j) + (1./2)*h*h*m_az(i*m_k+j);

      double pos_x = m_X(i*m_k+j+1);
      double pos_y = m_Y(i*m_k+j+1);
      double pos_z = m_Z(i*m_k+j+1);

      m_ax(i*m_k+j+1) = force_a(m_X,pos_x,pos_y,pos_z,i,j);
      m_Vx(i*m_k+j+1) = m_Vx(i*m_k+j) + (1./2)*h*(m_ax(i*m_k+j+1) + m_ax(i*m_k+j));

      m_ay(i*m_k+j+1) = force_a(m_Y,pos_x,pos_y,pos_z,i,j);
      m_Vy(i*m_k+j+1) = m_Vy(i*m_k+j) + (1./2)*h*(m_ay(i*m_k+j+1) + m_ay(i*m_k+j));

      m_az(i*m_k+j+1) = force_a(m_Z,pos_x,pos_y,pos_z,i,j);
      m_Vz(i*m_k+j+1) = m_Vz(i*m_k+j) + (1./2)*h*(m_az(i*m_k+j+1) + m_az(i*m_k+j));
    };
  };
};



/*void ParticleSolver::RK4_xupdate(double t, double x,double y, double v, double f1(double t, double x, double y, double v), double f2(double t, double x, double y, double v)){
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

 void ParticleSolver::EulerChromer(){
};

double get_force() {
}
*/
