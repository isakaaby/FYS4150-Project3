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
  cout << "Press 1 to run for Earth-Sun system \n";
  cout << "Press 2 to run for Earth-Jupiter-Sun system \n";
  cout << "Press 3 to run for all planets \n";           //asking which task you want to run
  cout << "Press 4 to run for Mercury-Sun system with relativistic correction \n";
  cout << "Enter number:" << " ";
  cin >> task;

  vector<string> all_names;
  all_names.push_back("Sun");
  all_names.push_back("Earth");
  all_names.push_back("Jupiter");
  all_names.push_back("Mars");
  all_names.push_back("Venus");
  all_names.push_back("Saturn");
  all_names.push_back("Mercury");
  all_names.push_back("Uranus");
  all_names.push_back("Neptune");
  all_names.push_back("Pluto");

  //int N;
  double beta = 2;
  double T;
  int N;
  int k;
  vector<string> planets;


  if (task==1){
    N = 2;
    T = 1;
    k = 10000;
    planets.push_back("Sun");
    planets.push_back("Earth");
    PlanetSolver solver;
    solver.init_sun_center(planets,beta,N,k,T);
    bool sun_center = true;
    solver.solvesystem_verlet(sun_center);
    solver.write_pos_to_file();
  }

  if (task==2){
    N = 3;
    T = 23;                  //orbit time for Jupiter
    k = 10000;
    planets.push_back("Sun");
    planets.push_back("Earth");
    planets.push_back("Jupiter");
    PlanetSolver solver;
    solver.init_sun_center(planets,beta,N,k,T);
    bool sun_center = true;
    solver.solvesystem_verlet(sun_center);
    solver.write_pos_to_file();
  }

  if (task==3){
    N = 10;
    T = 250;
    k = 10000;
    planets.push_back("Sun");
    planets.push_back("Earth");
    planets.push_back("Jupiter");
    planets.push_back("Mars");
    planets.push_back("Venus");
    planets.push_back("Saturn");
    planets.push_back("Mercury");
    planets.push_back("Uranus");
    planets.push_back("Neptune");
    planets.push_back("Pluto");
    PlanetSolver solver;
    solver.init(planets,beta,N,k,T);
    bool sun_center = false;
    solver.solvesystem_verlet(sun_center);
    solver.write_pos_to_file();
  }

  if (task==4){
    N = 2;
    T = 100;
    k = 100000;
    //T = 200.;          //orbit time for mercury
    //T = 24.1095;
    planets.push_back("Sun");
    planets.push_back("Mercury");
    MercurySunSolver solver;
    solver.init(planets,beta,N,k,T);
    solver.solve_mercury_sun_verlet();
    //solver.solve_mercury_sun_eulerchromer();
    solver.write_pos_to_file();
  }
}
