// tried to make a general class system which can be used on many body problems
// this uses class planets
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
  //Asks for which task to run.
  cout << "Press 1 to run for Earth-Sun system \n";
  cout << "Press 2 to run for Earth-Jupiter-Sun system \n";
  cout << "Press 3 to run for all planets \n";
  cout << "Press 4 to run for Mercury-Sun system\n";
  cout << "Press 5 to run a convergence test (Earth-Sun system) \n";
  cout << "Enter number:" << " ";
  cin >> task;

  //Creates array of solar system.
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

  //Initializing and creating values to be used as paramters down the line.
  int k;
  int N;
  double T;
  double beta = 2;

  vector<string> planets;
  vector<string> m_names;
  int run_verlet = 1;
  int run_euler_forward = 2;
  int run_euler_chromer = 3;


  //Code for Earth-Sun system.
  if (task==1){
    N = 2;                    //Number of planetary objects
    T = 1;                    //orbit time for earth
    k = 1e6;                  //time step
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);  //Filling array with names of the planets.
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);   //Initializing the system.
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);      //Solving the system using the velocity Verlet method.
    solver.write_pos_to_file();                     //Writing positions to files.
    solver.write_energy_to_file();                  //Writing potential, kinetic and total energy to files.

    int do_;
    double tol;
    cout << "Press 1 to test convervation laws \n";
    cout << "Enter number:" << " ";
    cin >> do_;
    if (do_ == 1){
      tol = 1e-07;
      solver.test_constant_energy(tol);     //Testing if energy is constant.

      tol = 1e-12;
      solver.test_constant_angular(tol);    //Testing if angular momentum is constant.

      tol = 1e-02;
      solver.test_circular_orbit(tol);      //Testing if the orbit is cicular.
    }
  }

  //Code for Earth-Sun-Jupiter system.
  if (task==2){
    N = 3;
    T = 23;         //orbit time for Jupiter
    k = 100000;     //Time step
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);  //Creating array for the planets.
    }
    PlanetSolver solver;
    bool no_error_test = false;
    solver.init_sun_center(m_names,beta,N,k,T,no_error_test);   //Initializing the system.
    bool sun_center = true;
    solver.solvesystem(sun_center,run_verlet);    //Solving the system using the velocity Verlet method.
    solver.write_pos_to_file();                   //Writing positions to files.
    solver.write_energy_to_file();                //Writing potential, kinetic and total energy to files.
  }

  //Code for all planets.
  if (task==3){
    N = 10;             //Number of planetary objects
    k = 1e4;            //Time step
    T = 250;            //Years (using earth years)

    //Filling array with all the planets.
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
    solver.init(planets,beta,N,k,T);        //Initializing the system.
    bool sun_center = false;
    solver.solvesystem(sun_center, run_euler_forward);    //Solving the system using the Euler-Cromer method.
    solver.write_pos_to_file();                           //Writing positions to files.
    solver.write_energy_to_file();                        //Writing potential, kinetic and total energy to files.
  }

  //Code for Mercury-Sun system.
  if (task==4){
    N = 2;     //Number of objects
    T = 100;  //orbit time for mercury
    k = 100;  //Time step

    m_names.push_back(all_names[0]);
    m_names.push_back(all_names[6]);
    MercurySunSolver solver;
    solver.init(m_names,beta,N,k,T);    //Initializing the system.
    solver.solve_mercury_sun_verlet();  //Solving the system using the velocity Verlet method.
    //solver.solve_mercury_sun_eulerchromer();
    solver.write_pos_to_file();         //Writing positions to files.
    solver.write_energy_to_file();      //Writing potential, kinetic and total energy to files.
  }

    //Code for convergence test for Earth-Sun system.
  if (task == 5){
    int N_experiments;
    cout << "Enter number of experiments:" << " ";
    cin >> N_experiments;

    int N = 2;      //Number of objects.
    T = 1;          //orbit time for earth
    int k = 12000;  //Time step
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);    //Filling array with planet names.
    }
    PlanetSolver solver;
    solver.test_convergence(m_names,beta,N,k,T,N_experiments, run_euler_chromer);   //Solver that tests the convergence of the Euler-Cromer method.
  }
}
