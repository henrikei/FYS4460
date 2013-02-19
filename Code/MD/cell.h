#ifndef CELL_H
#define CELL_H

#include "atom.h"
#include<armadillo>

using namespace arma;

class Cell
{
public:
    Cell(const ivec3 &, const double &);
    void addAtom(Atom*);
    void addNeighbour(Cell*);
    ivec3 getPositionIndices();
    vector<Atom*> getAtoms();
    vector<Cell*> getNeighbours();
    void setDistanceCorrection(const double &, const int &, const int &);
    vec3 getDistanceCorrection(const int &);
    void deleteAtoms();
    vector<Cell*> neighbours;
private:
    ivec3 positionIndices;
    double sideLength;
    vector<Atom*> atoms;
    mat distanceCorrection;
};

#endif // CELL_H
