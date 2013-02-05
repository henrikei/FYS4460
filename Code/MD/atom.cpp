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

void Atom::setForce(vec3 &f){
    force = f;
}

void Atom::addForce(vec3 &f){
    force += f;
}

vec3 Atom::getPosition(){
    return position;
}

vec3 Atom::getVelocity(){
    return velocity;
}

vec3 Atom::getForce(){
    return force;
}

double Atom::getMass(){
    return mass;
}

string Atom::getName()
{
    return name;
}
