#include <iostream>
#include "Pores/pores.h"
#include "Pores/cylindricalpore.h"
#include "Pores/circularpores.h"
#include "Pores/reducedensity.h"

using namespace std;

int main()
{
    ReduceDensity reducedensity;
    reducedensity.readFile("circularPores0.xyz");
    reducedensity.make(0.5);
    reducedensity.writeFile("cirularPoresHalfDensity.xyz");
//    double rMin = 20/3.405;
//    double rMax = 30/3.405;
//    int nBalls = 20;
//    double boxLength = 20*5.72/3.405;
//    CircularPores pore(boxLength);
//    pore.readFile("state2000.xyz");
//    pore.initializeFree();
//    pore.makeSeveral(nBalls, rMin, rMax);
//    pore.writeFile("circularPores0.xyz");
//    return 0;
}
