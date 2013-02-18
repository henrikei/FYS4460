#include <iostream>
#include <armadillo>
#include <sstream>
#include <cmath>
#include "atom.h"
#include "system.h"

using namespace std;
using namespace arma;

int main()
{
    // Input in SI units
    string atomType = "Ar";
    int nAtomsPerDim = 7;
    double mass = 39.948*1.66E-27;
    double fccLength = 5.260E-10;
    double temperature = 2000;
    double endTime = 1.0E-12;
    double timeStep = 1.0E-14;

    System test("Ar", nAtomsPerDim, mass, fccLength, temperature, endTime, timeStep);
    test.generate();



//    test.integrate();
//    test.writeVelHist();
    return 0;
}
