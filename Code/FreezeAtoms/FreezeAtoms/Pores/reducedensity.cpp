#include "reducedensity.h"

ReduceDensity::ReduceDensity()
{
}

void ReduceDensity::make(double reductionRatio){
    int reductionNumber = 0;
    for (int i = 0; i < nAtoms; i++){
        if (randu() < reductionRatio && free.at(i-reductionNumber)){
            positions.erase(positions.begin() + i - reductionNumber);
            velocities.erase(velocities.begin() + i - reductionNumber);
            free.erase(free.begin() + i - reductionNumber);
            reductionNumber += 1;
        }
    }
    nAtoms = nAtoms - reductionNumber;
}
