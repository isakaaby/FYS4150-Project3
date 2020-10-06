#include "particlesolver.hpp"
#include "planets.hpp"

void PlanetSolver::init(vector<string> names){
  vector<string> m_names = names; // std type vector

  //to be used in all derived classes
  vec masses = vec(m_N);
  vec distances = vec(m_N);
  vec params = vec(2);
  Planets Planet;
  for (int i = 0; i < m_N; ++i){
    params = Planet.initialize(m_names[i]);
    masses(i) = params(0);
    distances(i) = params(1);
  }
};

void PlanetSolver::solvesystem(){

};
