#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>

void MercurySunSolver::init(vector<string> names, double beta, int N, int k, double T){
  initialize(beta, N, k, T);

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);
  m_masses = zeros<vec>(m_N);

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
  }
  m_X(0) = m_Y(0) = m_Z(0) = 0.0;
  m_Vx(0) = m_Vy(0) = m_Vz(0) = 0.0;
  m_ax(0) = m_ay(0) = m_az(0) = 0.0;

  m_X(m_k) = 0.3075; m_Y(m_k) = m_Z(m_k) = 0.0;
  m_Vx(m_k) = 0.0; m_Vy(m_k) = 12.44; m_Vz(m_k) = 0.0;
  m_ax(m_k) = m_ay(m_k) = m_az(m_k) = 0.0;

};



void MercurySunSolver::solve_mercury_sun(){

};

void MercurySunSolver::write_pos_to_file(){
  ofstream x;
  ofstream y;
  ofstream z;

  string filename_1("./results/position_x" + to_string(m_k) + ".txt");
  string filename_2("./results/position_y" + to_string(m_k) + ".txt");
  string filename_3("./results/position_z" + to_string(m_k) + ".txt");
  x.open(filename_1);
  y.open(filename_2);
  z.open(filename_3);
  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      x << m_X(i*m_k+j) << " ";
      y << m_Y(i*m_k+j) << " ";
      z << m_Z(i*m_k+j) << " ";
    }
    x << "\n";
    y << "\n";
    z << "\n";
  }
  x.close();
  y.close();
  z.close();
}
