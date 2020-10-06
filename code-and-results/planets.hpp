#ifndef Planets_HPP
#define Planets_HPP

#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

//Solar system planets class
// Work in units: time [years], distance [AU], mass [solar masses]

class Planets {

private:


protected:

  double m_sun = 2E+30;
  double m_earth = (6E+24)/m_sun;
  //double earth = (6E+24)/m_sun;
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

  double sun;
  double mass;
  double distance_sun;

  vector<string> names;
  vector<double> system_mass;
  vector<double> system_r;


public:
  //double get_Mass(string name);      // calculate force between stellar objects
  //double get_Distance(string name);
  vec initialize(string planet);
  //void initialize(double m, double r);
};



#endif
