#include "cell.h"

Cell::Cell(const ivec3 &posInd, const double &sideL)
{
    positionIndices = posInd;
    sideLength = sideL;
    distanceCorrection = zeros<mat>(3,26);
}

void Cell::addAtom(Atom *a){
    atoms.push_back(a);
}

void Cell::addNeighbour(Cell *n){
    neighbours.push_back(n);
}

ivec3 Cell::getPositionIndices(){
    return positionIndices;
}

vector<Atom*> Cell::getAtoms(){
    return atoms;
}

vector<Cell*> Cell::getNeighbours(){
    return neighbours;
}

void Cell::setDistanceCorrection(const double &value, const int &i, const int &j){
    distanceCorrection(i,j) = value;
}

vec3 Cell::getDistanceCorrection(const int &n){
    return distanceCorrection.col(n);
}

void Cell::deleteAtoms(){
    atoms.clear();
}
