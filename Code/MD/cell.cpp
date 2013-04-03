#include "cell.h"

Cell::Cell(const ivec3 &posInd, const double &sideL)
{
    positionIndices = posInd;
    sideLength = sideL;
    for (int i = 0; i < 26; i++){
        distanceCorrection.push_back(zeros(3));
    }
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

vector<Atom*>& Cell::getAtoms(){
    return atoms;
}

vector<Cell*>& Cell::getNeighbours(){
    return neighbours;
}

void Cell::setDistanceCorrection(const double &value, const int &i, const int &j){
    distanceCorrection.at(j)(i) = value;
}

const vec3& Cell::getDistanceCorrection(const int &n){
    return distanceCorrection.at(n);
}

void Cell::deleteAtoms(){
    atoms.clear();
}
