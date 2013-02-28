#ifndef ANDERSENTHERMOSTAT_H
#define ANDERSENTHERMOSTAT_H

#include "system.h"
#include "atom.h"

class AndersenThermostat: public Modifier
{
public:
    AndersenThermostat(System *);
    void setTargetTemperature(double);
    void setRelaxationTime(double);
    void setTurnOffTime(double);
    void apply();
private:
    vector<Atom*> atoms;
    double targetTemperature;
    double tau;
    double tOff;
    double time;
    double newVelocity;
    int nAtoms;
};

#endif // ANDERSENTHERMOSTAT_H
