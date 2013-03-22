#ifndef CYLINDRICALPORE_H
#define CYLINDRICALPORE_H

#include "Pores/pores.h"



class CylindricalPore : public Pores
{
public:
    CylindricalPore();
    void make(vec2 center, double radius);
};

#endif // CYLINDRICALPORE_H
