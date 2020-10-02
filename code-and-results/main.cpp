// declare force here PS. Remember unit time [years], length in [AU] and mass [solar mass]
// tried to make a general class system which can be used on many body problems
// this uses class planets
#include<cmath>

double grav_force();

int main(int argc, char const *argv[]){
  //use classes here

}

double grav_force(double Ma, double Mb,double x, double y, double z):
  double m_G = 39.478 //AU^(3)*yr^(-2)*M^(-1)
  double g_force= m_G*(Ma*Mb)/cmath::sqrt(x*x + y*y + z*z)
