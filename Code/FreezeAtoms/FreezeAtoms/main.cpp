#include <iostream>
#include "Pores/pores.h"

using namespace std;

int main()
{
    CylindricalPore pore;
    vec2 center;
    center << 10*5.72/3.405 << 10*5.72/3.405;
    double radius = 20/3.405;
    pore.readFile("state400.xyz");
    pore.make(center, radius);
    pore.writeFile("cylindricalPore0.xyz");
    return 0;
}
