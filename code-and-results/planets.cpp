#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
#include<stdlib.h>
#include<armadillo>
using namespace std;
using namespace chrono;
using namespace arma;

<<<<<<< HEAD

double m_sun = 2E+30;
double m_earth = (6E+24)/m_sun;
double m_jupiter = (1.9E+27)/m_sun;
double m_mars = (6.6E23)/m_sun;
double m_venus = (4.9E24)/m_sun;
double m_saturn = (5.5E26)/m_sun;
double m_mercury = (3.3E23)/m_sun;
double m_uranus = (8.8E25)/m_sun;
double m_neptune = (1.03E26)/m_sun;
double m_pluto = (1.31E22)/m_sun;

//Defining distances from sun for the different planets in AU
double r_earth = 1;
double r_jupiter = 5.20;
double r_mars = 1.52;
double r_venus = 0.72;
double r_saturn = 9.54;
double r_mercury = 0.39;
double r_uranus = 19.19;
double r_neptune = 30.065;
double r_pluto = 39.53;

double mass;
double distance_sun;
double sun = 1;
=======
>>>>>>> 88476de5d03e874f28b3a9bb1f5cbeca349ade23

double sun = 1.;

vec Planets::initialize(string planet){
  //defining masses
  if(planet == "Earth"){
    mass = m_earth;
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
  l(0) = mass;
  l(1) = distance_sun;
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
<<<<<<< HEAD
>>>>>>> ff47c7ac2e14fce3b1ba7362944ca93824ac2612

=======
>>>>>>> 88476de5d03e874f28b3a9bb1f5cbeca349ade23
}
*/

/*double Planets::get_Distance(string name) {
  string names = "r_" + name;

  double value = atof(names.c_str());

  return value;
}*/
