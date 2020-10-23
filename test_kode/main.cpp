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
  cout << "Press 3 to run for all planets \n";
  cout << "Press 4 to run for Mercury-Sun system\n";
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
  int run_verlet = 1;
  int run_euler_forward = 2;
  int run_euler_chromer = 3;



  if (task==1){
    int N = 2;
    T = 30;                    //orbit time for earth
    k = 20000;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);
    solver.write_pos_to_file();

    int do_;
    double tol;
    cout << "Press 1 to test convervation laws \n";
    cout << "Enter number:" << " ";
    cin >> do_;
    if (do_ == 1){
      tol = 1e-07;
      solver.test_constant_energy(tol);

      tol = 1e-12;
      solver.test_constant_angular(tol);

      tol = 1e-03;
      solver.test_circular_orbit(tol);
    }
  }

  if (task==2){
    int N = 3;
    T = 23;                  //orbit time for Jupiter
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);
    solver.write_pos_to_file();

  }

  if (task==3){
    int N = 10;
    int k = 10000;
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
    solver.solvesystem(sun_center,run_verlet);
    solver.write_pos_to_file();
  }

  if (task==4){
    int N = 2;
    T = 100;
    k = 1000000;
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

    int N = 2;
    T = 100;                    //orbit time for earth
    int k = 10000;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.test_convergence(m_names,beta,N,k,T,N_experiments, run_euler_chromer);
  }
}
