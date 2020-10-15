#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>

void MercurySunSolver::init(vector<string> m_names, double beta, int N, int k, double T){
  initialize(beta, N, k, T);

  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);
  m_masses = zeros<vec>(m_N);

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
  };
  m_X(0) = m_Y(0) = m_Z(0) = 0.0;
  m_Vx(0) = m_Vy(0) = m_Vz(0) = 0.0;
  m_ax(0) = m_ay(0) = m_az(0) = 0.0;

  m_X(m_k) = 0.3075; m_Y(m_k) = m_Z(m_k) = 0.0;
  m_Vx(m_k) = 12.44; m_Vy(m_k) = m_Vz(m_k) = 0.0;
  m_ax(m_k) = m_ay(m_k) = m_az(m_k) = 0.0;

};

void MercurySunSolver::solve_mercury_sun(){

};

void MercurySunSolver::write_pos_to_file(){
};
