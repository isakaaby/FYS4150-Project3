#include "planets.hpp"
#include<iostream>
#include<string>
#include<chrono>
#include<stdlib.h>
#include<armadillo>
#include <fstream>
#include<cstdio>
using namespace std;
using namespace chrono;
using namespace arma;

void Planets::read_pos_vel(){
  int Nparticles = 10;
  x0 = new double[Nparticles];
  y0 = new double[Nparticles];
  z0 = new double[Nparticles];
  vx0 = new double[Nparticles];
  vy0 = new double[Nparticles];
  vz0 = new double[Nparticles];
  char const *filename_pos_vel = "./data/NASA_pos_vel.txt";

  //Open files
  FILE *fp_init = fopen(filename_pos_vel, "r"); //Open file to read, specified by "r".

  //Loop over each particle and extract its mass and initial conditions:
  for (int i = 0; i < Nparticles; i++){
    fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x0[i], &y0[i], &z0[i], &vx0[i], &vy0[i], &vz0[i]);
  }
  fclose(fp_init); //Close file with initial conditions

}

vec Planets::initialize(string planet){
  //Read off txt file from NASA here or in method below (return as vector)


  //defining masses
  double x, y, z, vx, vy, vz;
  if (planet == "Sun"){
    mass = sun;
    x = x0[0]; y = y0[0]; z = z0[0]; vx = vx0[0]; vy = vy0[0]; vz = vz0[0];
  }
  if(planet == "Earth"){
    mass = m_earth;
    x = x0[1]; y = y0[1]; z = z0[1]; vx = vx0[1]; vy = vy0[1]; vz = vz0[1];
  }
  if(planet == "Jupiter"){
    mass = m_jupiter;
    x = x0[2]; y = y0[2]; z = z0[2]; vx = vx0[2]; vy = vy0[2]; vz = vz0[2];
  }
  if(planet == "Mars") {
    mass = m_mars;
    x = x0[3]; y = y0[3]; z = z0[3]; vx = vx0[3]; vy = vy0[3]; vz = vz0[3];
  }
  if(planet == "Venus") {
    mass = m_venus;
    x = x0[4]; y = y0[4]; z = z0[4]; vx = vx0[4]; vy = vy0[4]; vz = vz0[4];
  }
  if(planet == "Saturn") {
    mass = m_saturn;
    x = x0[5]; y = y0[5]; z = z0[5]; vx = vx0[5]; vy = vy0[5]; vz = vz0[5];
  }
  if(planet == "Mercury") {
    mass = m_mercury;
    x = x0[6]; y = y0[6]; z = z0[6]; vx = vx0[6]; vy = vy0[6]; vz = vz0[6];
  }
  if(planet == "Uranus") {
    mass = m_uranus;
    x = x0[7]; y = y0[7]; z = z0[7]; vx = vx0[7]; vy = vy0[7]; vz = vz0[7];
  }
  if(planet == "Neptune") {
    mass = m_neptune;
    x = x0[8]; y = y0[8]; z = z0[8]; vx = vx0[8]; vy = vy0[8]; vz = vz0[8];
  }
  if(planet == "Pluto") {
    mass = m_pluto;
    x = x0[9]; y = y0[9]; z = z0[9]; vx = vx0[9]; vy = vy0[9]; vz = vz0[9];
  }
  vec l = zeros<vec>(7);
  l(0) = mass; l(1) = x; l(2) = y; l(3) = z; l(4) = vx; l(5) = vy; l(6) = vz;
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
}*/


/*double Planets::get_Distance(string name) {
  string names = "r_" + name;
  double value = atof(names.c_str());
  return value;
}*/
