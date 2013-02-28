#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H

#include "system.h"
#include "atom.h"

class BerendsenThermostat: public Modifier
{
public:
    BerendsenThermostat(System *);
    void setTargetTemperature(double);
    void setRelaxationTime(double);
    void setTurnOffTime(double);
    void apply();
private:
    vector<Atom*> atoms;
    double gamma;
    double targetTemperature;
    double tau;
    double tOff;
    double time;
    int nAtoms;
};

#endif // BERENDSENTHERMOSTAT_H
