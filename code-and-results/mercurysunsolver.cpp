#include "particlesolver.hpp"
#include "planets.hpp"
#include <vector>
#include <string>

void MercurySunSolver::init(vector<string> names, double beta, int N, int k, double T){
  initialize(beta, N, k, T);

  m_names = names;

  Planets Planet;
  Planet.read_pos_vel();
  vec params = vec(7);
  m_masses = zeros<vec>(m_N);

  for (int i = 0; i < m_N; i++){
    params = Planet.initialize(m_names[i]);
    m_masses(i) = params(0);
  }
  m_X(0) = m_Y(0) = m_Z(0) = 0.0;
  m_Vx(0) = m_Vy(0) = m_Vz(0) = 0.0;
  m_ax(0) = m_ay(0) = m_az(0) = 0.0;

  m_X(m_k) = 0.3075; m_Y(m_k) = m_Z(m_k) = 0.0;
  m_Vx(m_k) = 0.0; m_Vy(m_k) = 12.44; m_Vz(m_k) = 0.0;
  m_ax(m_k) = m_ay(m_k) = m_az(m_k) = 0.0;


  //double posMx,posMy,posMz,velMx,velMy,velMz;
  //double M = 0;
  //M = m_masses[0] + m_masses[1];
  //posMx = m_masses[1]*m_X(m_k);
  //posMy = m_masses[1]*m_Y(m_k);
  //posMz = m_masses[1]*m_Z(m_k);
  //velMx = m_masses[1]*m_Vx(m_k);
  //velMy = m_masses[1]*m_Vy(m_k);
  //velMz = m_masses[1]*m_Vz(m_k);

  //m_X(m_k) -= posMx/M;
  //m_Y(m_k) -= posMy/M;
  //m_Z(m_k) -= posMz/M;
  //m_Vx(m_k) -= velMx/M;
  //m_Vy(m_k) -= velMy/M;
  //m_Vz(m_k) -= velMz/M;
};

double MercurySunSolver::force_mercury_rel(vec pos, int l, int j){
  double G = 4*M_PI*M_PI; //AU^(3)*yr^(-2)*M(sol)^(-1);
  double mass_sun = m_masses(0);
  double c = 173*365;      //AU yr^(-1)
  double diffr,r,r_term,diffl,rel_term,a;
  diffr = m_X(l*m_k+j)*m_X(l*m_k+j) + m_Y(l*m_k+j)*m_Y(l*m_k+j) + m_Z(l*m_k+j)*m_Z(l*m_k+j);
  r = pow(diffr,0.5);
  r_term = pow(r,(m_beta+1));
  diffl = m_Lx(l*m_k+j)*m_Lx(l*m_k+j) + m_Ly(l*m_k+j)*m_Ly(l*m_k+j) + m_Lz(l*m_k+j)*m_Lz(l*m_k+j);
  rel_term = 3*diffl/(diffr*c*c);
  //calculate gravitational acceleration
  a = (pos(l*m_k+j)*G*mass_sun)/r_term + ((pos(l*m_k+j)*G*mass_sun)/r_term)*rel_term;

  return -a;
}



void MercurySunSolver::solve_mercury_sun_verlet(){
  double l = 1;
  for (int j = 0; j < m_k-1; j++){ // for time
    verlet_pos(l,j);
    m_ax(l*m_k+j+1) = force_a(m_X,l,j+1);
    m_ay(l*m_k+j+1) = force_a(m_Y,l,j+1);
    m_az(l*m_k+j+1) = force_a(m_Z,l,j+1);
    verlet_vel(l,j);
    get_angular_momentum(j+1);
    m_ax(l*m_k+j+1) = force_mercury_rel(m_X,l,j+1);
    m_ay(l*m_k+j+1) = force_mercury_rel(m_Y,l,j+1);
    m_az(l*m_k+j+1) = force_mercury_rel(m_Z,l,j+1);
    verlet_vel(l,j);
  }
};

void MercurySunSolver::solve_mercury_sun_eulerchromer(){
  double l = 1;
  for (int j = 0; j < m_k-1; j++){ // for time
    get_angular_momentum(j);
    m_ax(l*m_k+j) = force_mercury_rel(m_X,l,j);
    m_ay(l*m_k+j) = force_mercury_rel(m_Y,l,j);
    m_az(l*m_k+j) = force_mercury_rel(m_Z,l,j);
    eulerchromer(l,j);
  }
  cout << atan(m_Y(l*m_k+m_k-1)/m_X(l*m_k+m_k-1)) << "\n";
};



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
      x << m_X(i*m_k+j) << " ";
      y << m_Y(i*m_k+j) << " ";
      z << m_Z(i*m_k+j) << " ";
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
