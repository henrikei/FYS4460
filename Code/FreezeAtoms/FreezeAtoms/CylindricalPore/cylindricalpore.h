#ifndef CYLINDRICALPORE_H
#define CYLINDRICALPORE_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


class CylindricalPore
{
public:
    CylindricalPore();
    void readFile(string inFileName);
    void make(vec2 center, double radius);
    void writeFile(string outFileName);
private:
    int nAtoms;
    string text;
    string atomType;
    vector<vec3> positions;
    vector<vec3> velocities;
    vector<bool> free;
};

#endif // CYLINDRICALPORE_H
