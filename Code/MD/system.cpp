#include <iostream>
#include <armadillo>
#include <sstream>
#include "system.h"
#include "atom.h"
#include "cell.h"
#include "lib.h"

using namespace std;
using namespace arma;


System::System(string a, int N, double m, double b, double T, double tEnd, double tStep)
{
    double m0 = 39.948*1.66E-27;
    double sigma = 3.405E-10;
    double T0 = 119.8;
    double t0 = 2.1569E-12;
    double v0 = sigma/t0;

    atomType = a;
    nAtomsPerDim = N;
    atomMass = m/m0;
    dist = b/sigma;
    temperature = T/T0;
    timeEnd = tEnd/t0;
    timeStep = tStep/t0;
    cellSize = 2.5;
}


void System::generate(){
    // Initialize atoms with initial positions and velocities
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    double velocityStdDev = sqrt(temperature/atomMass);
    vec3 r = zeros(3);
    vec3 r_mod = zeros(3);
    vec3 v = zeros(3);
    vec3 dr = zeros(3);
    for (int i = 0; i < nAtomsPerDim; i++){
        r(0) = i*dist;
        for (int j = 0; j < nAtomsPerDim; j++){
            r(1) = j*dist;
            for (int k = 0; k < nAtomsPerDim; k++){
                r(2) = k*dist;
                v = velocityStdDev*randn(3);
                atoms.push_back(new Atom(atomType, atomMass, r, v));
                dr << dist/2 << dist/2 << 0;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms.push_back(new Atom(atomType, atomMass, r_mod, v));
                dr << 0 << dist/2 << dist/2;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms.push_back(new Atom(atomType, atomMass, r_mod, v));
                dr << dist/2 << 0 << dist/2;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms.push_back(new Atom(atomType, atomMass, r_mod, v));
            }
        }
    }
    // Remove translational drift
    vec3 velocitySum = zeros(3);
    for (int i = 0; i < nAtoms; i++){
        velocitySum =+ atoms.at(i)->getVelocity();
    }
    velocitySum = -velocitySum/nAtoms;
    for (int i = 0; i < nAtoms; i++){
        atoms.at(i)->addVelocity(velocitySum);
    }


    // Initialize cells

    // cell positions (position = cell indices (integers))
    int nCellsPerDim = (int) (nAtomsPerDim*dist/cellSize);
    // redefine cellSize
    cellSize = nAtomsPerDim*dist/nCellsPerDim;
    nCells = nCellsPerDim*nCellsPerDim*nCellsPerDim;
    ivec3 positionIndices;
    for (int i = 0; i < nCellsPerDim; i++){
        for (int j = 0; j < nCellsPerDim; j++){
            for (int k = 0; k < nCellsPerDim; k++){
                positionIndices << i << j << k;
                cells.push_back(new Cell(positionIndices, cellSize));
            }
        }
    }

    // add neighbours
    imat directionVecs = zeros<imat>(3,26);
    int counter = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                if (!(i == 1 && j == 1 && k == 1))
                {
                    directionVecs(0, counter) = i-1;
                    directionVecs(1, counter) = j-1;
                    directionVecs(2, counter) = k-1;
                    counter += 1;}
            }
        }
    }
    for (int i = 0; i < nCells; i++){
        for (int j = 0; j < 26; j++){
            ivec3 neighbourPointer = cells.at(i)->getPositionIndices() + directionVecs.col(j);
            for (int k = 0; k < 3; k++){
                if (neighbourPointer(k) < 0){
                    neighbourPointer(k) += nCellsPerDim;
                    cells.at(i)->setDistanceCorrection(-nCellsPerDim*cellSize, k, j);
                } else if (neighbourPointer(k) >= nCellsPerDim){
                    neighbourPointer(k) -= nCellsPerDim;
                    cells.at(i)->setDistanceCorrection(nCellsPerDim*cellSize, k, j);
                }
            }
            for (int k = 0; k < nCells; k++){
                if (neighbourPointer(0) == cells.at(k)->getPositionIndices()(0) &&
                    neighbourPointer(1) == cells.at(k)->getPositionIndices()(1) &&
                    neighbourPointer(2) == cells.at(k)->getPositionIndices()(2)){
                    cells.at(i)->addNeighbour(cells.at(k));
                }
            }
        }
    }
    populateCells();
}

void System::writeState(string fn){
    string filename = fn;
    ofstream ofile;
    ofile.open(filename);
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    ofile << nAtoms << endl;
    ofile << "Lattice of Argon atoms" << endl;
    for (int i = 0; i < nAtoms; i++){
        ofile << atoms.at(i)->getName() << " ";
        for (int j = 0; j < 3; j++){
            ofile << atoms.at(i)->getPosition()(j) << " ";
        }
        for (int j = 0; j < 3; j++){
            ofile << atoms.at(i)->getVelocity()(j) << " ";
        }
        ofile << endl;
    }
    ofile.close();
}

void System::integrate(){
    double boxLength = nAtomsPerDim*dist;
    double time = 0;
    vec3 velMidpoint;
    vec3 newPos;
    vec3 newVel;
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    int nSteps = (int)(timeEnd/timeStep);
    ofstream energyOut;
    energyOut.open("out/energy.dat");
    writeEnergy(energyOut, time);
    writeState("out/state0.xyz");
    calculateForce();
    for (int i = 0; i < nSteps; i++){
        time += timeStep;
        for (int j = 0; j < nAtoms; j++){
            velMidpoint = atoms.at(j)->getVelocity() + atoms.at(j)->getForce()*timeStep/2;
            atoms.at(j)->setVelocity(velMidpoint);
            newPos = atoms.at(j)->getPosition() + velMidpoint*timeStep;
            newPos(0) = fmod(newPos(0) + 5*boxLength, boxLength);
            newPos(1) = fmod(newPos(1) + 5*boxLength, boxLength);
            newPos(2) = fmod(newPos(2) + 5*boxLength, boxLength);
            atoms.at(j)->setPosition(newPos);
        }
        populateCells();
        calculateForce();
        for (int j = 0; j < nAtoms; j++){
            newVel = atoms.at(j)->getVelocity() + atoms.at(j)->getForce()*timeStep/2;
            atoms.at(j)->setVelocity(newVel);
        }
        writeEnergy(energyOut, time);
        stringstream outName;
        outName << "out/state" << (i+1) << ".xyz";
        writeState(outName.str());
    }
    energyOut.close();
}


void System::calculateForce(){
    for (int i = 0; i < nCells; i++){
        vector<Atom*> residents = cells.at(i)->getAtoms();
        vector<Cell*> neighbourCells = cells.at(i)->getNeighbours();
        int nResidents = residents.size();
        int nNeighbourCells = neighbourCells.size();
        vec3 Force = zeros(3);
        vec3 radialVec = zeros(3);
        double radialDist = 0;
        for (int j = 0; j < nResidents; j++){
            residents.at(j)->setForce(Force);
        }
        for (int j = 0; j < nResidents; j++){
            for (int k = j + 1; k < nResidents; k++){
                radialVec = residents.at(j)->getPosition() - residents.at(k)->getPosition();
                radialDist = sqrt(dot(radialVec, radialVec));
                Force = 24*(2/(pow(radialDist,14)) - 1/pow(radialDist,8))*radialVec;
                residents.at(j)->addForce(Force);
                Force = -Force;
                residents.at(k)->addForce(Force);
            }
            for (int k = 0; k < nNeighbourCells; k++){
                vector<Atom*> neighbourAtoms = neighbourCells.at(k)->getAtoms();
                int nNeighbourAtoms = neighbourAtoms.size();
                for (int l = 0; l < nNeighbourAtoms; l++){
                    radialVec = residents.at(j)->getPosition() - (neighbourAtoms.at(l)->getPosition() + cells.at(i)->getDistanceCorrection(k));
                    radialDist = sqrt(dot(radialVec, radialVec));
                    Force = 24*(2/(pow(radialDist,14)) - 1/pow(radialDist,8))*radialVec;
                    residents.at(j)->addForce(Force);
                }
            }
        }
    }
}

void System::populateCells(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    for (int i = 0; i < nCells; i++){
        cells.at(i)->deleteAtoms();
        for (int j = 0; j < nAtoms; j++){
            ivec3 cellIndices = cells.at(i)->getPositionIndices();
            ivec3 atomIndices;
            for (int k = 0; k < 3; k++){
                atomIndices(k) = (int) atoms.at(j)->getPosition()(k)/cellSize;
            }
            if (cellIndices(0) <= atomIndices(0) && cellIndices(0) + 1 > atomIndices(0) &&
                cellIndices(1) <= atomIndices(1) && cellIndices(1) + 1 > atomIndices(1) &&
                cellIndices(2) <= atomIndices(2) && cellIndices(2) + 1 > atomIndices(2)){
                cells.at(i)->addAtom(atoms.at(j));
            }
        }
    }
}

double System::getKineticEnergy(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    double kineticEnergy = 0;
    for (int i = 0; i < nAtoms; i++){
        vec3 velocity = atoms.at(i)->getVelocity();
        kineticEnergy += 0.5*atomMass*dot(velocity, velocity);
    }
    return kineticEnergy;
}

double System::getPotentialEnergy(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    vec3 radialVec = zeros(3);
    double radialDist = 0;
    double potentialEnergy = 0;
    for (int i = 0; i < nAtoms; i++){
        for (int j = i + 1; j < nAtoms; j++){
            radialVec = atoms.at(i)->getPosition() - atoms.at(j)->getPosition();
            radialDist = sqrt(dot(radialVec, radialVec));
            potentialEnergy += 4*(1/pow(radialDist,12) - 1/pow(radialDist,6));
        }
    }
    return potentialEnergy;
}

void System::writeVelHist(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    double velocityStdDev = sqrt(temperature/atomMass);
    double maxVelocity = 4*velocityStdDev;
    int nBins = 100;
    double binSize1 = 2*maxVelocity/nBins;
    double binSize2 = binSize1/2;
    double binLower1 = -maxVelocity;
    double binLower2 = 0;
    vec velx = zeros(nBins); vec vely = zeros(nBins); vec velz = zeros(nBins); vec velMagnitude = zeros(nBins);
    int nx = 0; int ny = 0; int nz = 0; int nmag = 0;
    double velMag = 0;
    ofstream ofile;
    ofile.open("velHist.dat");
    for(int i = 0; i < nBins; i++){
        for(int j = 0; j< nAtoms; j++){
            if (binLower1 < atoms.at(j)->getVelocity()(0) && atoms.at(j)->getVelocity()(0) <= (binLower1 + binSize1))
                velx(i) += 1; nx += 1;
            if (binLower1 < atoms.at(j)->getVelocity()(1) && atoms.at(j)->getVelocity()(1) <= (binLower1 + binSize1))
                vely(i) += 1; ny +=1;
            if (binLower1 < atoms.at(j)->getVelocity()(2) && atoms.at(j)->getVelocity()(2) <= (binLower1 + binSize1))
                velz(i) += 1; nz += 1;
            velMag = sqrt(atoms.at(j)->getVelocity()(0)*atoms.at(j)->getVelocity()(0) + atoms.at(j)->getVelocity()(1)*atoms.at(j)->getVelocity()(1) + atoms.at(j)->getVelocity()(2)*atoms.at(j)->getVelocity()(2));
            if (binLower2 < velMag && velMag <= binLower2 + binSize2)
                velMagnitude(i) +=1; nmag += 1;
        }
        binLower1 += binSize1;
        binLower2 += binSize2;
    }
    binLower1 = -maxVelocity;
    binLower2 = 0;
    for (int i = 0; i < nBins; i++)
        ofile << (binLower1 + (i+1)*binSize1)/2 << " " << velx(i)/(nx*binSize1) << " " << vely(i)/(ny*binSize1)
              << " " << velz(i)/(nz*binSize1) << " " << (binLower2 + (i+1)*binSize2)/2 << " " << velMagnitude(i)/(nmag*binSize2) << endl;
    ofile.close();
}

void System::writeEnergy(ofstream &ofile, double time){
    ofile << time << "  "  << getKineticEnergy() << "  " << getPotentialEnergy() << endl;
}
