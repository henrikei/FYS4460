#include <iostream>
#include <fstream>
#include <armadillo>


using namespace arma;
using namespace std;

int main()
{
    // Load positions from .xyz file
    ifstream inFile;
    inFile.open("state300_50K.xyz");
    if(!inFile) { // file couldn't be opened
          cerr << "Error: file could not be opened" << endl;
          exit(1);
       }
    int nAtoms;
    string dummy;
    vec3 position;
    vec3 velocity;
    vector<vec3> positions;
    vector<vec3> velocities;
    inFile >> nAtoms >> dummy >> dummy >> dummy >> dummy;
    while (!inFile.eof()){
        inFile >> dummy >> position[0] >> position[1] >> position[2] >> velocity[0] >> velocity[1] >> velocity[2];
        positions.push_back(position);
        velocities.push_back(velocity);
    }

    // Minimal image convension
    mat minimalImageConv = zeros(3,27);
    double boxLength = 8*5.26/3.405;
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

    // Find distances
    int nDistances = nAtoms*nAtoms/2;
    vec distances = ones(nDistances)*50*boxLength;
    vec3 distanceVecTest = zeros(3);
    double distanceTest = 0;

    counter = 0;
    for (uint i = 0; i < positions.size()-1; i++){
        for (uint j = i + 1; j < positions.size()-1; j++){
            for (int k = 0; k < 27; k++){
                distanceVecTest = positions.at(i) - positions.at(j) + minimalImageConv.col(k);
                distanceTest = sqrt(distanceVecTest(0)*distanceVecTest(0) + distanceVecTest(1)*distanceVecTest(1) + distanceVecTest(2)*distanceVecTest(2));
                if (distanceTest < distances(counter)){
                    distances(counter) = distanceTest;
                }
            }
            counter += 1;
        }
    }

    // write distances to file
    ofstream outFile;
    outFile.open("distances.dat");
    for (int i = 0; i < nDistances; i++){
        outFile << distances(i) << endl;
    }
    cout << nDistances << "  " << positions.size();
    return 0;
}

