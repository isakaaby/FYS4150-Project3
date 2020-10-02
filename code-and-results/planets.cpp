#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
using namespace std;
using namespace chrono;

void Planets::mass(){
  //defining masses
  double m_sun = 2*10**30;
  double m_earth = (6*10**(24))/m_sun;
  double m_jupiter = (1.9*10**(27))/m_sun;
  double m_mars = (6.6*10**(23))/m_sun;
  double m_venus = (4.9*10**(24))/m_sun;
  double m_saturn = (5.5*10**(26))/m_sun;
  double m_mercury = (3.3*10**(23))/m_sun;
  double m_uranus = (8.8*10**(25))/m_sun;
  double m_neptun = (1.03*10**(26))/m_sun;
  double m_pluto = (1.31*10**(22))/m_sun;

  double m_sun = 1;
}

void Planets::distance(){
  //defining distances from sun for the different planets in AU
  double r_earth = 1;
  double r_jupiter = 5.20;
  double r_mars = 1.52;
  double r_venus = 0.72;
  double r_saturn = 9.54;
  double r_mercury = 0.39;
  double r_uranus = 19.19;
  double r_neptun = 30.065;
  double r_pluto = 39.53;

}
