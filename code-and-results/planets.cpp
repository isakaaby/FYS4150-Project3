#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;

void Planets::intialize(double mass1, double mass2, double r){
  m_mass1 = mass1;
  m_mass2 = mass2;
  m_r = r;
}

void Planets::g_force(){
  force = (mass1*mass2*G)/(m_r*m_r);


  //defining distances from sun for the different planets in AU


}
