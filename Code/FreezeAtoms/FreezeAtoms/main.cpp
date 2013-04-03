#include <iostream>
#include "Pores/pores.h"
#include "Pores/cylindricalpore.h"
#include "Pores/circularpores.h"
#include "Pores/reducedensity.h"

using namespace std;

int main()
{
    vec2 center;
    center << 6*5.72/3.405 << 6*5.72/3.405;
    double radius = 14/3.405;
    CylindricalPore pore;
    pore.readFile("state1000_12x12x12.xyz");
    pore.make(center, radius);
    pore.writeFile("cylindricalPoreR14.xyz");

    ReduceDensity red;
    double reductionRatio = 0.5;
    red.readFile("cylindricalPoreR14.xyz");
    red.make(reductionRatio);
    red.writeFile("cylindricalPoreR14HalfDensity.xyz");

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
