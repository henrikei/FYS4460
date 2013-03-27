#ifndef REDUCEDENSITY_H
#define REDUCEDENSITY_H

#include "Pores/pores.h"


class ReduceDensity : public Pores
{
public:
    ReduceDensity();
    void make(double reductionRatio);
};

#endif // REDUCEDENSITY_H
