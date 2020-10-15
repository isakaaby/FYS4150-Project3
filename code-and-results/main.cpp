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
  cout << "Press 4 to run for Mercury and the Sun only \n";
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
  double T;

  if (task==1){
    int N = 2;
    T = 1;                    //orbit time for earth
    vector<string> m_names;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solvesystem();
    solver.write_pos_to_file();



  }
  if (task==2){
    int N = 3;
    T = 12;                  //orbit time for Jupiter
    vector<string> m_names;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solvesystem();
  }

  if (task==3){
    int N = 10;
    T = 250;               //orbit time for Pluto
    vector<string> m_names;
    for (int i = 0; i < N; i++){
      m_names.push_back(all_names[i]);
    }
    PlanetSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solvesystem();
    solver.write_pos_to_file();
  }

  if (task==4){
    int N = 2;
    T = 1;           //orbit time for mercury
    vector<string> m_names;
    m_names.push_back(all_names[0]);
    m_names.push_back(all_names[6]);
    PlanetSolver solver;
    solver.init(m_names,beta,N,k,T);
    solver.solvesystem();
    solver.write_pos_to_file();
  }

}

/*
double grav_force(double Ma, double Mb,double x, double y, double z):
  double m_G = 39.478 //AU^(3)*yr^(-2)*M^(-1)
  double g_force= m_G*(Ma*Mb)/cmath::sqrt(x*x + y*y + z*z)
*/
