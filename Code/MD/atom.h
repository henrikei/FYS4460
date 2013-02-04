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
    void setForce(vec3 &);
    void addForce(vec3 &);
    vec3 getPosition();
    vec3 getVelocity();
    vec3 getForce();
    double getMass();
    string getName();
private:
    string name;
    double mass;
    vec3 position;
    vec3 velocity;
    vec3 force;
};

#endif // ATOM_H
