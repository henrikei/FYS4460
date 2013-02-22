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
    int nAtomsPerDim = 8;
    double mass = 39.948*1.66E-27;
    double fccLength = 5.090E-10;
    double temperature = 120;
    double endTime = 2.0E-12;
    double timeStep = 2.0E-14;

    System test("Ar", nAtomsPerDim, mass, fccLength, temperature, endTime, timeStep);
    test.generate();
    test.integrate();
//    test.writeVelHist();
    return 0;
}
