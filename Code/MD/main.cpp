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
    System test("Ar", 7, 39.948, 1.7246, 2);
    test.generate();
    test.integrate(0.1,0.001);
    test.writeVelHist();
    return 0;
}
