#ifndef CIRCULARPORES_H
#define CIRCULARPORES_H

#include "Pores/pores.h"


class CircularPores : public Pores
{
public:
    CircularPores(double boxLen);
    void makeSingle(vec3 center, double radius);
    void makeSeveral(int numOf, double rMin, double rMax);
private:
    double boxLength;
    mat minimalImage;
};

#endif // CIRCULARPORES_H
