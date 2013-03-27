#include "atom.h"

Atom::Atom()
{
}

Atom::Atom(string &a, double &m, vec3 &r, vec3 &v)
{
    name = a;
    mass = m;
    position = r;
    velocity = v;
    force = zeros(3);
    displacement = zeros(3);
    free = 1;
}

void Atom::setPosition(vec3 &r){
    position = r;
}

void Atom::setVelocity(vec3 &v){
    velocity = v;
}

void Atom::addVelocity(vec3 &v){
    velocity += v;
}

void Atom::addDisplacement(vec3 &d){
    displacement += d;
}

void Atom::setForce(vec3 &f){
    force = f;
}

void Atom::addForce(vec3 &f){
    force(0) += f(0);
    force(1) += f(1);
    force(2) += f(2);
}

void Atom::setFree(int &f){
    free = f;
}

vec3& Atom::getPosition(){
    return position;
}

vec3& Atom::getVelocity(){
    return velocity;
}

vec3& Atom::getDisplacement(){
    return displacement;
}

vec3& Atom::getForce(){
    return force;
}

int Atom::getFree(){
    return free;
}

double& Atom::getMass(){
    return mass;
}

string& Atom::getName()
{
    return name;
}
