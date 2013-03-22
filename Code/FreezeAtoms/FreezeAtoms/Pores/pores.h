#ifndef PORES_H
#define PORES_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


class Pores
{
public:
    Pores();
    void readFile(string inFileName);
    void writeFile(string outFileName);
    void initializeFree();
protected:
    int nAtoms;
    string text;
    string atomType;
    vector<vec3> positions;
    vector<vec3> velocities;
    ivec free;
};

#endif // PORES_H
