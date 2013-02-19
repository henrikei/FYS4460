#ifndef SYSTEM_H
#define SYSTEM_H

#include<iostream>
#include<armadillo>
#include"atom.h"
#include"cell.h"

using namespace std;
using namespace arma;

class System
{
public:
    System(string, int, double, double, double, double, double);
    void generate();
    void writeState(string);
    void integrate();
    void calculateForce();
    void populateCells();
    double getKineticEnergy();
    double getPotentialEnergy();
    void writeVelHist();
    void writeEnergy(ofstream &, double);
private:
    string atomType;
    int nAtomsPerDim;
    double atomMass;
    double dist;
    double temperature;
    double timeEnd;
    double timeStep;
    double cellSize;
    int nCells;
    vector<Atom*> atoms;
    vector<Cell*> cells;
};

#endif // SYSTEM_H
