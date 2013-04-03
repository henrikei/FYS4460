#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Atom
{
public:
    Atom();
    Atom(string &, double &, vec3 &, vec3 &);
    void setPosition(vec3 &);
    void setVelocity(vec3 &);
    void addVelocity(vec3 &);
    void addDisplacement(vec3 &);
    void setForce(vec3 &);
    void addForce(vec3 &);
    void setFree(int &);
    const vec3 &getPosition();
    const vec3 &getVelocity();
    const vec3 &getDisplacement();
    const vec3 &getForce();
    const int &getFree();
    const double &getMass();
    const string &getName();
private:
    string name;
    double mass;
    vec3 position;
    vec3 velocity;
    vec3 displacement;
    vec3 force;
    int free;
};

#endif // ATOM_H
