#include "Pores/cylindricalpore.h"

CylindricalPore::CylindricalPore()
{
}

void CylindricalPore::make(vec2 center, double radius){
    vec2 distance;
    for (int i = 0; i < nAtoms; i++){
        distance << (positions.at(i)(0) - center(0)) << (positions.at(i)(1) - center(1));
        if (dot(distance, distance) > radius*radius){
            free(i) = 0;
        }
    }
}
