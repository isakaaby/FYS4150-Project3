#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
#include<stdlib.h>
#include<armadillo>
using namespace std;
using namespace chrono;
using namespace arma;


double sun = 1.;

vec Planets::initialize(string planet){
  //defining masses
  if(planet == "Earth"){
    mass = m_earth;
    cout << mass << endl;
    distance_sun = r_earth;
  } else if(planet == "Jupiter"){
    mass = m_jupiter;
    distance_sun = r_jupiter;
  } else if(planet == "Mars") {
    mass = m_mars;
    distance_sun = r_mars;
  } else if(planet == "Venus") {
    mass = m_venus;
    distance_sun = r_venus;
  } else if(planet == "Saturn") {
    mass = m_saturn;
    distance_sun = r_saturn;
  } else if(planet == "Mercury") {
    mass = m_mercury;
    distance_sun = r_mercury;
  } else if(planet == "Uranus") {
    mass = m_uranus;
    distance_sun = r_uranus;
  } else if(planet == "Neptune") {
    mass = m_neptune;
    distance_sun = r_neptune;
  } else if(planet == "Pluto") {
    mass = m_pluto;
    distance_sun = r_pluto;
  }
  vec l = zeros<vec>(2);
  l(0) = m_earth;
  l(1) = distance_sun;
  cout << (6E+24)/(2E+30) << endl;
  //cout << mass << endl;
  return l;
}

//void Planets::initialize(double m, double r){
//  mass = m;
//  distance_sun = r;
//}


/*double Planets::get_Mass(string name) {
  double names = "m_" + name;
  //cout << names << "\n";
  double value = atof(names.cstr());
  //cout << names << "\n";
  return value;
}

double Planets::get_Distance(string name) {
  string names = "r_" + name;

  double value = atof(names.c_str());

  return value;
}*/
