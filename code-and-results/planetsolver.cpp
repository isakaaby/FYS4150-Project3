#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>

void PlanetSolver::init(double beta, int N, int k, int T){
  initialize(beta, N, k, T);

  vector<string> all_names;
  all_names.push_back("Sun");
  all_names.push_back("Earth");
  all_names.push_back("Jupiter");
  all_names.push_back("Mars");
  all_names.push_back("Venus");
  all_names.push_back("Saturn");
  all_names.push_back("Mercury");
  all_names.push_back("Uranus");
  all_names.push_back("Neptune");
  all_names.push_back("Pluto");

  vector<string> m_names;

  for (int i = 0; i < m_N; i++){
    m_names.push_back(all_names[i]);
  }
  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);

  m_masses = zeros<vec>(m_N);

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
  }

  /*double M = 0;
  double posMx = 0;
  double posMy = 0;
  double posMz = 0;
  double velx_sun = 0;
  double vely_sun = 0;
  double velz_sun = 0;

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
    m_X(i*m_k) = params(1); m_Y(i*m_k) = params(2); m_Z(i*m_k) = params(3);
    m_Vx(i*m_k) = 365*params(4); m_Vy(i*m_k) = 365*params(5); m_Vz(i*m_k) = 365*params(6);
    m_ax(i*m_k) = m_ay(i) = m_az(i*m_k) = 0.0;
    M += m_masses[i];
    posMx += m_masses[i]*m_X(i*m_k);
    posMy += m_masses[i]*m_Y(i*m_k);
    posMz += m_masses[i]*m_Z(i*m_k);
    velx_sun -= m_masses[i]*m_Vx(i*m_k);
    vely_sun -= m_masses[i]*m_Vy(i*m_k);
    velz_sun -= m_masses[i]*m_Vz(i*m_k);
  }
  //m_Vx(0) = velx_sun;
  //m_Vy(0) = vely_sun;
  //m_Vz(0) = velz_sun;

  for (int i = 0; i < m_N; i++){
    m_X(i*m_k) -= posMx/M;
    m_Y(i*m_k) -= posMy/M;
    m_Z(i*m_k) -= posMz/M;
  }*/

  m_X(0*m_k) = 0.; m_Y(0*m_k) = 0.; m_Z(0*m_k) = 0.;
  m_X(1*m_k) = 1.; m_Y(1*m_k) = 0.; m_Z(1*m_k) = 0.;
  m_Vx(0*m_k) = 0.; m_Vy(0*m_k) = 0.; m_Vz(0*m_k) = 0.;
  m_Vx(1*m_k) = 2*M_PI; m_Vy(1*m_k) = 0.; m_Vz(1*m_k) = 0.;


}

void PlanetSolver::solvesystem(){
  verlet();

};

void PlanetSolver::write_pos_to_file(){
  ofstream x;
  ofstream y;
  ofstream z;

  string filename_1("./results/position_x.txt");
  string filename_2("./results/position_y.txt");
  string filename_3("./results/position_z.txt");
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