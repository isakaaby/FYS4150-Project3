#include "catch.hpp"
#include "particlesolver.hpp"
#include "planetsolver.hpp"
#include "planets.hpp"  // Remove if unnecesseray
#include<iostream>

//TEST_CASE("Potential Energy")

//TEST_CASE("Kinetic Energy")

TEST_CASE("Conservation of total energy"){

}

TEST_CASE("Conservation of angular momentum"){
  REQUIRE((val - exact) < tol);
}
