// declare force here PS. Remember unit time [years], length in [AU] and mass [solar mass]
// tried to make a general class system which can be used on many body problems
// this uses class planets
#include<cmath>
#include<armadillo>
#include "particlesolver.hpp"
#include <vector>
#include <string>

using namespace arma;
using namespace std;

//double grav_force();

int main(int argc, char const *argv[]){
  int task;
  cout << "Press 1 to run for Earth-Sun system \n";
  cout << "Press 2 to run for Earth-Jupiter-Sun system \n";
  cout << "Press 3 to run for all planets \n";           //asking which task you want to run
  cout << "Enter number:" << " ";
  cin >> task;

  //int N;
  int k = 1000;
  double beta = 2;
  double T = 1;

  if (task==1){
    int N = 2;
    PlanetSolver solver;
    solver.init(beta,N,k,T);
    solver.solvesystem();
    solver.write_pos_to_file();

  }
  if (task==2){
    int N = 3;
    PlanetSolver solver;
    solver.init(beta,N,k,T);
    solver.solvesystem();
  }

  if (task==3){
    int N = 10;
    PlanetSolver solver;
    solver.init(beta,N,k,T);
    solver.solvesystem();
    solver.write_pos_to_file();
  }


}

/*
double grav_force(double Ma, double Mb,double x, double y, double z):
  double m_G = 39.478 //AU^(3)*yr^(-2)*M^(-1)
  double g_force= m_G*(Ma*Mb)/cmath::sqrt(x*x + y*y + z*z)
*/
