#include "jacobimethodsolver.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;

//Setting up the superclass for Jacobi's method with rotational algorithm to be used in all derived classes

mat SolarSystem::initialize(int N){
  m_N = N;
  //initialize variables to set up Jacobis algorithm
  //to be used in all derived classes

}


void SolarSystem::grav_force(){


}

void SolarSystem::Verlet(){


}

void SolarSystem::EulerChromer(){


}

void SolarSystem::RungeKutta4(){


}
