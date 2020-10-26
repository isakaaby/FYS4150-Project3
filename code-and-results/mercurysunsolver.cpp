#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>
#include <iomanip>

//Method that initializes the solver, returns nothing.
void MercurySunSolver::init(vector<string> names, double beta, int N, int k, double T){
  initialize(beta, N, k, T);    //Calls on particlesolver.cpp's initializer.

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();      //Reads NASA positions and velocities from .txt file.
  vec params = vec(7);
  m_masses = zeros<vec>(m_N);

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);   //Initializing Mercury and the Sun.
    m_masses(i) = params(0);
  }

  m_X(0) = m_Y(0) = m_Z(0) = 0.0;
  m_Vx(0) = m_Vy(0) = m_Vz(0) = 0.0;
  m_ax(0) = m_ay(0) = m_az(0) = 0.0;

  m_X(m_k) = 0.3075; m_Y(m_k) = m_Z(m_k) = 0.0;
  m_Vx(m_k) = 0.0; m_Vy(m_k) = 12.44; m_Vz(m_k) = 0.0;
  m_ax(m_k) = m_ay(m_k) = m_az(m_k) = 0.0;


  /*double posMx,posMy,posMz,velMx,velMy,velMz;
  double M = 0;
  M = m_masses[0] + m_masses[1];
  posMx = m_masses[1]*m_X(m_k);
  posMy = m_masses[1]*m_Y(m_k);
  posMz = m_masses[1]*m_Z(m_k);
  velMx = m_masses[1]*m_Vx(m_k);
  velMy = m_masses[1]*m_Vy(m_k);
  velMz = m_masses[1]*m_Vz(m_k);
  m_X(m_k) -= posMx/M;
  m_Y(m_k) -= posMy/M;
  m_Z(m_k) -= posMz/M;
  m_Vx(m_k) -= velMx/M;
  m_Vy(m_k) -= velMy/M;
  m_Vz(m_k) -= velMz/M;*/
};

//Method that computes the gravitational force, returns nothing.
void MercurySunSolver::force_mercury_rel(int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double mass_sun = m_masses(0);
  double V = 0;
  double diffl = m_Lx(l*m_k)*m_Lx(l*m_k) + m_Ly(l*m_k)*m_Ly(l*m_k) + m_Lz(l*m_k)*m_Lz(l*m_k);
  double c = 173*365;      //AU yr^(-1)
  double cc = c*c;
  double diffr,r,r_term,rel_term,a;
  diffr = m_X(l*m_k+j)*m_X(l*m_k+j) + m_Y(l*m_k+j)*m_Y(l*m_k+j) + m_Z(l*m_k+j)*m_Z(l*m_k+j);
  r = pow(diffr,0.5);
  r_term = pow(r,(m_beta+1));
  rel_term = 3*diffl/(diffr*cc);

  //calculate gravitational acceleration
  double ax = (m_X(l*m_k+j)*G*mass_sun)/r_term + ((m_X(l*m_k+j)*G*mass_sun)/r_term)*rel_term;
  double ay = (m_Y(l*m_k+j)*G*mass_sun)/r_term + ((m_Y(l*m_k+j)*G*mass_sun)/r_term)*rel_term;
  double az = (m_Z(l*m_k+j)*G*mass_sun)/r_term + ((m_Z(l*m_k+j)*G*mass_sun)/r_term)*rel_term;
  m_ax(l*m_k+j) = -ax;
  m_ay(l*m_k+j) = -ay;
  m_az(l*m_k+j) = -az;
}


//Method that solves the problem using the velocity Verlet method, returns nothing.
void MercurySunSolver::solve_mercury_sun_verlet(){
  double l = 1;
  get_angular_momentum(0);    //Calculating the angular momentum.
  for (int j = 0; j < m_k-1; j++){ // for time
    verlet_pos(l,j);                //Updating the positions.
    force_mercury_rel(l,j+1);       //Calculating the gravitational force.
    verlet_vel(l,j);                //Updating the velocities.
  }


  //double min_distance = min(r_vec);
  //uvec indices = find(r_vec == min(r_vec));
  //cout << indices << "\n";
  //int I = 40236;
  //cout << m_X(l*m_k+I) << "\n";
  //cout << m_Y(l*m_k+I) << "\n";

  //cout << min_distance << "\n";
  //double theta_rad = atan(1);
  //cout << theta_rad << "\n";
  //cout << theta_rad*3600 << "\n";
};

//Mehod that calculates the angular momentum, returns angular momentum.
double MercurySunSolver::angular_momentum(double pos1, double v1, double pos2, double v2){
  return pos1*v2 - pos2*v1;
}

//Method that updates the angular momentum, returns nothing.
void MercurySunSolver::get_angular_momentum(int j){
  for (int i = 0; i < m_N; i++){ //for planets
    m_Lx(i*m_k+j) = angular_momentum(m_Y(i*m_k+j), m_Vy(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Ly(i*m_k+j) = angular_momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Z(i*m_k+j), m_Vz(i*m_k+j));
    m_Lz(i*m_k+j) = angular_momentum(m_X(i*m_k+j), m_Vx(i*m_k+j), m_Y(i*m_k+j), m_Vy(i*m_k+j));
  }
}
/*void MercurySunSolver::solve_mercury_sun_eulerchromer(){
  double l = 1;
  vec r_vec = vec(m_k);
  r_vec(0) = m_X(m_k);
  get_angular_momentum(0);
  for (int j = 0; j < m_k-1; j++){ // for time
    vec ux = force_mercury_rel(m_X,l,j);
    r_vec(j+1) = ux(0);
    m_ax(l*m_k+j) = ux(1);
    vec uy = force_mercury_rel(m_Y,l,j);
    m_ay(l*m_k+j) = uy(1);
    vec uz = force_mercury_rel(m_Z,l,j);
    m_az(l*m_k+j) = uz(1);
    eulerchromer(l,j);
  }
  //double min_distance = min(r_vec);
  //uvec indices = find(r_vec == min(r_vec));
  //cout << indices << "\n";
  //int I = 98858;
  //cout << m_X(l*m_k+I) << "\n";
  //cout << m_Y(l*m_k+I) << "\n";
  //cout << r_vec << "\n";
  //cout << min_distance << "\n";
  //double theta_rad = atan(1);
  //cout << theta_rad << "\n";
  //cout << theta_rad*3600 << "\n";
};*/


//Method that writes the positions to a file, returns nothing.
void MercurySunSolver::write_pos_to_file(){
  ofstream x;
  ofstream y;
  ofstream z;
  ofstream planet_names;

  string filename_1("./results/position_x.txt");
  string filename_2("./results/position_y.txt");
  string filename_3("./results/position_z.txt");
  string filename_4("./results/planet_names.txt");
  x.open(filename_1);
  y.open(filename_2);
  z.open(filename_3);
  planet_names.open(filename_4);
  for (int j = 0; j < m_k; j++){
    for (int i = 0; i < m_N; i++){
      x << setprecision(15) << m_X(i*m_k+j) << " ";
      y << setprecision(15) << m_Y(i*m_k+j) << " ";
      z << setprecision(15) << m_Z(i*m_k+j) << " ";
    }
    x << "\n";
    y << "\n";
    z << "\n";
  }
  for (int i = 0; i < m_N; i++){
    planet_names << m_names[i] << "\n";
  }
  x.close();
  y.close();
  z.close();
  planet_names.close();
}

//Method that writes the angular momentum to file, returns nothing. 
void MercurySunSolver::write_angular_momentum_to_file(vec x, vec y, vec z) {
  ofstream x_;
  ofstream y_;
  ofstream z_;

  string filename_1("./results/ang_x.txt");
  string filename_2("./results/ang_y.txt");
  string filename_3("./results/ang_z.txt");

  x_.open(filename_1);
  y_.open(filename_2);
  z_.open(filename_3);

  for (int j = 0; j < m_N; j++) {
    for (int i = 0; i < m_N; i++){
      x_ << x(i*m_k + j) << " ";
      y_ << y(i*m_k + j) << " ";
      z_ << z(i*m_k + j) << " ";
    }
    x_ << "\n";
    y_ << "\n";
    z_ << "\n";

  }
  x_.close();
  y_.close();
  z_.close();
}
