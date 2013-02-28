#ifndef MODIFIER_H
#define MODIFIER_H

#include "atom.h"
class System;

class Modifier
{
public:
    Modifier(System *);
    virtual void apply()=0;
protected:
    System *system;
};

#endif // MODIFIER_H
