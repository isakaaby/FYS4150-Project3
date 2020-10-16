// declare force here PS. Remember unit time [years], length in [AU] and mass [solar mass]
// tried to make a general class system which can be used on many body problems
// this uses class planets
#include <cmath>
#include <armadillo>
#include "particlesolver.hpp"
#include <vector>
#include <string>

using namespace arma;
using namespace std;

//double grav_force();
void menu();

int main(int argc, char const *argv[]){
  menu();

}

void menu() {
  int task;
  int NASA;
  cout << "Press 0 to run the program with initialized values from NASA and 1 to input position, velocity and mass yourself \n";
  cin >> NASA;
  cout << "Press 1 to run for Earth-Sun system \n";
  cout << "Press 2 to run for Earth-Jupiter-Sun system \n";
  cout << "Press 3 to run for all planets \n";           //asking which task you want to run
  cout << "Press 4 to run for a system of preferred objects \n";
  cout << "Enter number:" << " ";
  cin >> task;

  int N;
  int k = 3;
  double beta = 2;
  double T = 250;
  vector<string> planets;
  vector<double> x, y, z, vx, vy, vz, masses;
  double x_, y_, z_, vx_, vy_, vz_, mass;


  if (task==1){
    int N = 2;
    if(NASA == 0) {
      planets.push_back("Earth");
      planets.push_back("Sun");
      PlanetSolver solver;
      solver.init(beta,N,k,T, planets);
      solver.solvesystem();
      solver.write_pos_to_file();
    }



  }
  if (task==2){
    int N = 3;
    planets.push_back("Earth");
    planets.push_back("Sun");
    planets.push_back("Jupiter");
    PlanetSolver solver;
    solver.init(beta,N,k,T, planets);
    solver.solvesystem();
  }

  if (task==3){
    int N = 10;
    planets.push_back("Earth");
    planets.push_back("Sun");
    planets.push_back("Jupiter");
    planets.push_back("Mars");
    planets.push_back("Venus");
    planets.push_back("Saturn");
    planets.push_back("Mercury");
    planets.push_back("Uranus");
    planets.push_back("Neptune");
    planets.push_back("Pluto");
    PlanetSolver solver;
    solver.init(beta,N,k,T, planets);
    solver.solvesystem();
    solver.write_pos_to_file();
  }

  if(task == 4) {
    string param;
    cout << "Number of objects ";
    cin >> N;
    cout << "Enter planets \n";
    if (NASA == 1) {
      for(int i = 0; i < N; i++) {
        cin >> param >> x_ >> y_ >> z_ >> vx_ >> vy_ >> vz_ >> mass;
        planets.push_back(param);
        x.push_back(x_);
        y.push_back(y_);
        z.push_back(z_);
        vx.push_back(vx_);
        vy.push_back(vy_);
        vz.push_back(vz_);
        masses.push_back(mass);
      }
    }else {
      for(int i = 0; i < N; i++) {
        cin >> param;
        planets.push_back(param);
      }
    }

  PlanetSolver solver;
    solver.init(beta, N, k, T, planets, x, y, z, vx, vy, vz, masses);
    solver.solvesystem();
    solver.write_pos_to_file();

    //int N = input("")
  }
}

/*
double grav_force(double Ma, double Mb,double x, double y, double z):
  double m_G = 39.478 //AU^(3)*yr^(-2)*M^(-1)
  double g_force= m_G*(Ma*Mb)/cmath::sqrt(x*x + y*y + z*z)
*/
