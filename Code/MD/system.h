#ifndef SYSTEM_H
#define SYSTEM_H

#include<iostream>
#include<armadillo>
#include"atom.h"
#include"cell.h"
#include"Modifier/modifier.h"

using namespace std;
using namespace arma;

class System
{
public:
    System(string, int, double, double, double, double, double);
    void generate();
    void writeState(string);
    void readState(string);
    void integrate();
    void calculateForce();
    void populateCells();
    void addModifier(Modifier *);
    void addExternalForce(vec3);
    double getTimeStep();
    int getNumberOfAtoms();
    vector<Atom*> getAtoms();
    double getKineticEnergy();
    double getPotentialEnergy();
    double getTemperature();
    double getMeanSquareDisplacement();
    void writeVelHist();
    void writeObservables(ofstream &, double);
private:
    string atomType;
    int nAtomsPerDim;
    int nAtoms;
    int nAtomsFree;
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
    vector<Modifier*> modifiers;
    vec3 extForce;

    // Conversion factors
    double m0;                      // mass
    double sigma;                   // length
    double T0;                      // temperature
    double t0;                      // time
    double epsilon;                 // energy
};

#endif // SYSTEM_H
