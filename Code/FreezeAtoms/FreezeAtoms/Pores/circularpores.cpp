#include "circularpores.h"

CircularPores::CircularPores(double boxLen)
{
    boxLength = boxLen;
    minimalImage = zeros(3,27);
    int counter = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j ++){
            for (int k = 0; k < 3; k++){
                minimalImage(0, counter) = (i-1)*boxLength;
                minimalImage(1, counter) = (j-1)*boxLength;
                minimalImage(2, counter) = (k-1)*boxLength;
                counter += 1;
            }
        }
    }
}

void CircularPores::makeSingle(vec3 center, double radius){
    vec3 distanceVecTest;
    double distanceTest;
    for (int i = 0; i < nAtoms; i++){
        double distance = 10*boxLength*boxLength;
        for (int j = 0; j < 27; j++){
            distanceVecTest << (positions.at(i)(0) - center(0) - minimalImage(0,j)) << (positions.at(i)(1) - center(1) - minimalImage(1,j))
                            << (positions.at(i)(2) - center(2) - minimalImage(2,j));
            distanceTest = dot(distanceVecTest, distanceVecTest);
            if (distanceTest < distance){
                distance = distanceTest;
            }
        }
        if (distance < radius*radius){
            free.at(i) = 0;
        }
    }
}

void CircularPores::makeSeveral(int numOf, double rMin, double rMax){
    vec3 center;
    double radius;
    for (int i = 0; i < numOf; i++){
        radius = rMin + randu()*(rMax - rMin);
        center = randu(3)*boxLength;
        makeSingle(center, radius);
    }
}
