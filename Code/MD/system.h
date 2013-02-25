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
    double getMeanSquareDisplacement();
    void writeVelHist();
    void writeObservables(ofstream &, double);
private:
    string atomType;
    int nAtomsPerDim;
    int nAtoms;
    double atomMass;
    double fccLength;
    double temperature;
    double pressure;
    double timeEnd;
    double timeStep;
    double cellSize;
    int nCells;
    vector<Atom*> atoms;
    vector<Cell*> cells;

    // Conversion factors
    double m0;                      // mass
    double sigma;                   // length
    double T0;                      // temperature
    double t0;                      // time
    double epsilon;                 // energy
};

#endif // SYSTEM_H
