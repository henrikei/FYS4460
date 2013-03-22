#include "circularpores.h"

CircularPores::CircularPores()
{
}

void CircularPores::makeSingle(vec3 center, double radius){
    vec3 distance;
    for (int i = 0; i < nAtoms; i++){
        distance << (positions.at(i)(0) - center(0)) << (positions.at(i)(1) - center(1))
                 << (positions.at(i)(2) - center(2));
        if (dot(distance, distance) < radius*radius){
            free(i) = 0;
        }
    }
}

void CircularPores::makeSeveral(int numOf, double rMin, double rMax, double boxLength){
    vec3 center;
    double radius;
    for (int i = 0; i < numOf; i++){
        radius = rMin + randu()*(rMax - rMin);
        center = randu(3)*boxLength;
        makeSingle(center, radius);
        cout << "Hei" << endl;
    }
}
