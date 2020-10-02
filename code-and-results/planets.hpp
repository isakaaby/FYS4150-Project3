#ifndef Planets_HPP
#define Planets_HPP

#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

//Solar system planets class
// Work in units: time [years], distance [AU], mass [solar masses]

class Planets {

private:

protected:



public:
  void mass();      // calculate force between stellar objects
  void distance();
};



#endif
