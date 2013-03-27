#include "pores.h"
#include <fstream>

Pores::Pores()
{
}

void Pores::readFile(string inFileName){
    vec3 position = zeros(3);
    vec3 velocity = zeros(3);
    ifstream inFile;
    int freee;
    inFile.open(inFileName.c_str());
    if (!inFile){ // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    inFile >> nAtoms;
    getline(inFile, text); // For unkown reason, two commands are needed...
    getline(inFile, text);
    int counter = 0;
    while(counter < nAtoms){
        inFile >> atomType >> position[0] >> position[1] >> position[2]
               >> velocity[0] >> velocity[1] >> velocity[2] >> freee;
        positions.push_back(position);
        velocities.push_back(velocity);
        free.push_back(freee);
        counter += 1;
    }
}

void Pores::writeFile(string outFileName){
    ofstream outFile;
    outFile.open(outFileName.c_str());
    outFile << nAtoms << endl;
    outFile << text << endl;
    for (int i = 0; i < nAtoms; i++){
        outFile << atomType <<"  "<< positions.at(i)(0) <<"  "<< positions.at(i)(1)<<"  "<< positions.at(i)(2) <<"  "
                << velocities.at(i)(0) <<"  "<< velocities.at(i)(1) <<"  "<< velocities.at(i)(2) <<"  "<< free.at(i) << endl;
    }
}
