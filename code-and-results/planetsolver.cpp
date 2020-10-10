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

  m_masses = vec(m_N);

  for (int i = 0; i < m_N; ++i){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
    m_X(i) = params(1); m_Y(i) = params(2); m_Z(i) = params(3);
    m_Vx(i) = params(4); m_Vy(i) = params(5); m_Vz(i) = params(6);
    m_ax(i) = m_ay(i) = m_az(i) = 0;
  }


};

void PlanetSolver::solvesystem(){
  verlet();

};

void PlanetSolver::write_pos_to_file(){
  ofstream myfile;
  string filename("./results/Earth_sun.txt");
  myfile.open(filename);
  myfile << m_X << "\n";
  myfile << m_Y << "\n";
  myfile.close();

}
