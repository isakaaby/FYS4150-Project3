#include <cmath>
#include <armadillo>
#include "particlesolver.hpp"
#include <vector>
#include <string>

using namespace arma;
using namespace std;

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
  int k;
  int N;
  double T;
  double beta = 2;
  int do_;
  double tol;

  vector<string> planets;
  vector<string> m_names;
  int run_verlet = 1;
  int run_euler_forward = 2;
  int run_euler_chromer = 3;



  if (task==1){
    N = 2;
    T = 1;                    //orbit time for earth
    k = 10e6;
    double beta = 2.0;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);
    solver.write_pos_to_file();

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
    N = 3;
    T = 50;
    k = 1e5;     //orbit time for Jupiter
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);
    solver.write_pos_to_file();

    cout << "Press 1 to test convervation laws \n";
    cout << "Enter number:" << " ";
    cin >> do_;
    if (do_ == 1){
      tol = 1e-07;
      solver.test_constant_energy(tol);

      tol = 1e-12;
      solver.test_constant_angular(tol);

      tol = 1e-02;
      solver.test_circular_orbit(tol);
    }
  }

  if (task==3){
    N = 10;
    k = 1e4;
    T = 250;
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
    N = 2;
    T = 100;
    k = 1e6;
    m_names.push_back(all_names[0]);
    m_names.push_back(all_names[6]);
    MercurySunSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solve_mercury_sun_verlet();
    solver.write_pos_to_file();

  }

  if (task == 5){
    int N_experiments;
    cout << "Enter number of experiments:" << " ";
    cin >> N_experiments;

    int N = 2;
    T = 1;
    int k = 50000;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.test_convergence(m_names,beta,N,k,T,N_experiments, run_euler_forward);
  }
}
