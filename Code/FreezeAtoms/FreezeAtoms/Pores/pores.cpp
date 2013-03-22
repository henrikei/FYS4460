#include "cylindricalpore.h"
//#include <string>
#include <fstream>

CylindricalPore::CylindricalPore()
{
}

void CylindricalPore::readFile(string inFileName){
    vec3 position = zeros(3);
    vec3 velocity = zeros(3);
    ifstream inFile;
    inFile.open(inFileName.c_str());
    if (!inFile){ // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    inFile >> nAtoms;
    getline(inFile, text); // For unkown reason, two commands are needed...
    getline(inFile, text);
    while(!inFile.eof()){
        inFile >> atomType >> position[0] >> position[1] >> position[2]
               >> velocity[0] >> velocity[1] >> velocity[2];
        positions.push_back(position);
        velocities.push_back(velocity);
    }
}

void CylindricalPore::make(vec2 center, double radius){
    vec2 distance;
    for (int i = 0; i < nAtoms; i++){
        distance << (positions.at(i)(0) - center(0)) << (positions.at(i)(1) - center(1));
        if (dot(distance, distance) > radius*radius){
            free.push_back(false);
        } else {
            free.push_back(true);
        }
    }
}

void CylindricalPore::writeFile(string outFileName){
    ofstream outFile;
    outFile.open(outFileName.c_str());
    outFile << nAtoms << endl;
    outFile << text << endl;
    for (int i = 0; i < nAtoms; i++){
        outFile << atomType <<"  "<< positions.at(i)(0) <<"  "<< positions.at(i)(1)<<"  "<< positions.at(i)(2) <<"  "
                << velocities.at(i)(0) <<"  "<< velocities.at(i)(1) <<"  "<< velocities.at(i)(2) <<"  "<< free.at(i) << endl;
    }
}
