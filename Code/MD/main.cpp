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
    System test("Ar", 8, 39.948, 1.7246, 2);
    test.generate();
    test.integrate(0.05,0.0001);
    test.writeVelHist();
    return 0;
}
