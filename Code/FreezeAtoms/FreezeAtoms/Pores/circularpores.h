#ifndef CIRCULARPORES_H
#define CIRCULARPORES_H

#include "Pores/pores.h"


class CircularPores : public Pores
{
public:
    CircularPores();
    void makeSingle(vec3 center, double radius);
    void makeSeveral(int numOf, double rMin, double rMax, double boxLength);
};

#endif // CIRCULARPORES_H
