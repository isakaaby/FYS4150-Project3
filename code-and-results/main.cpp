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
  cout << "Press 4 to run for a system of preferred objects \n";
  cout << "Press 5 to run a convergence test (Earth-Sun system) \n";
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
  int k = 10000;
  double beta = 2;
  double T = 250;
  vector<string> planets;
  vector<string> m_names;



  if (task==1){
    int N = 2;
    T = 1;                    //orbit time for earth
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.init_sun_center(m_names,beta,N,k,T);
    bool sun_center = true;
    solver.solvesystem(sun_center);
    solver.write_pos_to_file();

    int do_;
    cout << "Press 1 to test conservation of energy\n";
    cout << "Press 2 to test conservation of angular momentum\n";
    cin >> do_;
    if (do_ == 1){
      solver.test_constant_energy();
    }

    if (do_ == 2){
      solver.test_constant_angular();
    }
  }

  if (task==2){
    int N = 3;
    T = 23;                  //orbit time for Jupiter
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.init_sun_center(m_names,beta,N,k,T);
    bool sun_center = true;
    solver.solvesystem(sun_center);
    solver.write_pos_to_file();

  }

  if (task==3){
    int N = 10;
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
    solver.init(m_names,beta,N,k,T);
    bool sun_center = false;
    solver.solvesystem(sun_center);
    solver.write_pos_to_file();
  }

  if (task==4){
    int N = 2;
    T = 100;
    //T = 200.;          //orbit time for mercury
    //T = 24.1095;
    m_names.push_back(all_names[0]);
    m_names.push_back(all_names[6]);
    MercurySunSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solve_mercury_sun_verlet();
    //solver.solve_mercury_sun_eulerchromer();
    solver.write_pos_to_file();

    //int N = input("")
  }

  if (task == 5){
    int N_experiments;
    cout << "Enter number of experiments:" << " ";
    cin >> N_experiments;
    cout << "Use 2pi as inital velocity (1axis) to get circular motion \n";

    int N = 2;
    T = 1;                    //orbit time for earth
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.test_convergence(m_names,beta,N,k,T,N_experiments);
  }
}

/*
double grav_force(double Ma, double Mb,double x, double y, double z):
  double m_G = 39.478 //AU^(3)*yr^(-2)*M^(-1)
  double g_force= m_G*(Ma*Mb)/cmath::sqrt(x*x + y*y + z*z)
*/
