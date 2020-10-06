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
  double sun;
  double m_earth, r_earth;
  double m_jupiter, r_jupiter;
  double m_mars, r_mars;
  double m_venus, r_venus;
  double m_saturn, r_saturn;
  double m_mercury, r_mercury;
  double m_uranus, r_uranus;
  double m_neptune, r_neptune;
  double m_pluto, r_pluto;

  double mass;
  double distance_sun;

  vector<string> names;
  vector<double> system_mass;
  vector<double> system_r;


public:
  //double get_Mass(string name);      // calculate force between stellar objects
  //double get_Distance(string name);
  vec initialize(string planet);
  void initialize(double m, double r);
};



#endif
