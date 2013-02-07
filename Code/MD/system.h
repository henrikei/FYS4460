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
    System(string, int, double, double, double);
    void generate();
    void writeState(string);
    void integrate(double, double);
    void calculateForce(Atom **, int);
    void writeVelHist();
    mat minimalImageConv;
private:
    string atomType;
    int nAtomsPerDim;
    double atomMass;
    double dist;
    double temperature;
    Atom **atoms;
};

#endif // SYSTEM_H
