#include <iostream>
#include <armadillo>
#include <sstream>
#include "system.h"
#include "atom.h"
#include "lib.h"

using namespace std;
using namespace arma;


System::System(string a, int N, double m, double b, double T)
{
    atomType = a;
    nAtomsPerDim = N;
    atomMass = m; // In atomic units
    dist = b;
    temperature = T;
    minimalImageConv = zeros(3,27);
    double boxLength = nAtomsPerDim*dist;
    int counter = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                minimalImageConv(0,counter) = (i-1)*boxLength;
                minimalImageConv(1,counter) = (j-1)*boxLength;
                minimalImageConv(2,counter) = (k-1)*boxLength;
                counter += 1;
            }
        }
    }
}


void System::generate(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    double velocityStdDev = 52.6632*sqrt(temperature/atomMass);
    atoms = new Atom*[nAtoms];
    vec3 r = zeros(3);
    vec3 r_mod = zeros(3);
    vec3 v = zeros(3);
    vec3 dr = zeros(3);
    int counter = 0;
    for (int i = 0; i < nAtomsPerDim; i++){
        r(0) = i*dist;
        for (int j = 0; j < nAtomsPerDim; j++){
            r(1) = j*dist;
            for (int k = 0; k < nAtomsPerDim; k++){
                r(2) = k*dist;
                v = velocityStdDev*randn(3);
                atoms[counter] = new Atom(atomType, atomMass, r, v);
                counter += 1;
                dr << dist/2 << dist/2 << 0;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms[counter] = new Atom(atomType, atomMass, r_mod, v);
                counter += 1;
                dr << 0 << dist/2 << dist/2;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms[counter] = new Atom(atomType, atomMass, r_mod, v);
                counter += 1;
                dr << dist/2 << 0 << dist/2;
                r_mod = r + dr;
                v = velocityStdDev*randn(3);
                atoms[counter] = new Atom(atomType, atomMass, r_mod, v);
                counter += 1;
            }
        }
    }
    // Remove translational drift
    vec3 velocitySum = zeros(3);
    for (int i = 0; i < nAtoms; i++){
        velocitySum =+ atoms[i]->getVelocity();
    }
    velocitySum = -velocitySum/nAtoms;
    for (int i = 0; i < nAtoms; i++){
        atoms[i]->addVelocity(velocitySum);
    }
}

void System::writeState(string fn){
    string filename = fn;
    ofstream ofile;
    ofile.open(filename);
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    ofile << nAtoms << endl;
    ofile << "Lattice of Argon atoms" << endl;
    for (int i = 0; i < nAtoms; i++){
        ofile << atoms[i]->getName() << " ";
        for (int j = 0; j < 3; j++){
            ofile << atoms[i]->getPosition()(j) << " ";
        }
        for (int j = 0; j < 3; j++){
            ofile << atoms[i]->getVelocity()(j) << " ";
        }
        ofile << endl;
    }
    ofile.close();
}

void System::integrate(double endtime, double timestep){
    double t = endtime;
    double dt = timestep;
    double boxLength = nAtomsPerDim*dist;
    vec3 velMidpoint;
    vec3 newPos;
    vec3 newVel;
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    int nSteps = (int)(t/dt);
    writeState("out/state0.xyz");
    calculateForce(atoms, nAtoms);
    for (int i = 0; i < nSteps; i++){
        for (int j = 0; j < nAtoms; j++){
            velMidpoint = atoms[j]->getVelocity() + atoms[j]->getForce()*dt/2;
            atoms[j]->setVelocity(velMidpoint);
            newPos = atoms[j]->getPosition() + velMidpoint*dt;
            newPos(0) = fmod(newPos(0) + 5*boxLength, boxLength);
            newPos(1) = fmod(newPos(1) + 5*boxLength, boxLength);
            newPos(2) = fmod(newPos(2) + 5*boxLength, boxLength);
            atoms[j]->setPosition(newPos);
        }
        calculateForce(atoms, nAtoms);
        for (int j = 0; j < nAtoms; j++){
            newVel = atoms[j]->getVelocity() + atoms[j]->getForce()*dt/2;
            atoms[j]->setVelocity(newVel);
        }
        stringstream outName;
        outName << "out/state" << (i+1) << ".xyz";
        writeState(outName.str());
    }
}

void System::calculateForce(Atom **a, int n){
    Atom **atoms = a;
    int nAtoms = n;
    vec3 initForce = zeros(3);
    vec3 radialVec;
    vec3 radialVecTest;
    double radialDist;
    double radialDistTest;
    for (int i = 0; i < nAtoms; i++){
        atoms[i]->setForce(initForce);
    }
    for (int i = 0; i < nAtoms; i++){
        for (int j = i + 1; j < nAtoms; j++){
            radialDist = nAtomsPerDim*dist;
            for (int k = 0; k < 27; k++){
                radialVecTest = atoms[i]->getPosition() - atoms[j]->getPosition() + minimalImageConv.col(k);
                radialDistTest = sqrt(dot(radialVecTest,radialVecTest));
                if (radialDistTest <= radialDist){
                    radialVec = radialVecTest;
                    radialDist = radialDistTest;
                }
            }
            initForce = 24*(2/(pow(radialDist,14)) - 1/pow(radialDist,8))*radialVec;
            atoms[i]->addForce(initForce);
            initForce = -initForce;
            atoms[j]->addForce(initForce);
        }
    }
}

void System::writeVelHist(){
    int nAtoms = 4*nAtomsPerDim*nAtomsPerDim*nAtomsPerDim;
    double velocityStdDev = 52.6632*sqrt(temperature/atomMass);
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
            if (binLower1 < atoms[j]->getVelocity()(0) && atoms[j]->getVelocity()(0) <= (binLower1 + binSize1))
                velx(i) += 1; nx += 1;
            if (binLower1 < atoms[j]->getVelocity()(1) && atoms[j]->getVelocity()(1) <= (binLower1 + binSize1))
                vely(i) += 1; ny +=1;
            if (binLower1 < atoms[j]->getVelocity()(2) && atoms[j]->getVelocity()(2) <= (binLower1 + binSize1))
                velz(i) += 1; nz += 1;
            velMag = sqrt(atoms[j]->getVelocity()(0)*atoms[j]->getVelocity()(0) + atoms[j]->getVelocity()(1)*atoms[j]->getVelocity()(1) + atoms[j]->getVelocity()(2)*atoms[j]->getVelocity()(2));
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
