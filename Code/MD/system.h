#ifndef SYSTEM_H
#define SYSTEM_H

#include<iostream>
#include<armadillo>
#include"atom.h"

using namespace std;
using namespace arma;

class System
{
public:
    System(string, int, double, double, double, double, double);
    void generate();
    void writeState(string);
    void integrate();
    void calculateForce(int);
    void writeVelHist();
    mat minimalImageConv;
private:
    string atomType;
    int nAtomsPerDim;
    double atomMass;
    double dist;
    double temperature;
    double timeEnd;
    double timeStep;
    vector<Atom*> atoms;
};

#endif // SYSTEM_H
