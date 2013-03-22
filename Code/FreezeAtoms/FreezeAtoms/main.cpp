#include <iostream>
#include "Pores/pores.h"
#include "Pores/cylindricalpore.h"
#include "Pores/circularpores.h"

using namespace std;

int main()
{
    CircularPores pore;
    double rMin = 20/3.405;
    double rMax = 30/3.405;
    int nBalls = 20;
    double boxLength = 20*5.72/3.405;
    pore.readFile("state400.xyz");
    pore.initializeFree();
    pore.makeSeveral(nBalls, rMin, rMax, boxLength);
    pore.writeFile("circularPore0.xyz");
    return 0;
}
