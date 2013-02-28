#include "berendsenthermostat.h"
#include "system.h"
#include "atom.h"

BerendsenThermostat::BerendsenThermostat(System *system) :
    Modifier(system)
{
    tau = 10;
    time = 0;
}

void BerendsenThermostat::setTargetTemperature(double T){
    targetTemperature = T/119.74;
}

void BerendsenThermostat::setRelaxationTime(double relaxationTime){
    tau = relaxationTime;
}

void BerendsenThermostat::setTurnOffTime(double t){
    tOff = t/(2.1569E-12);
}

void BerendsenThermostat::apply(){
    if (time < tOff){
        atoms = system->getAtoms();
        nAtoms = system->getNumberOfAtoms();
        gamma = sqrt(1 + (targetTemperature/system->getTemperature() - 1)/tau);
        for (int i = 0; i < nAtoms; i++){
            vec3 newVelocity = gamma*atoms.at(i)->getVelocity();
            atoms.at(i)->setVelocity(newVelocity);
        }
        time += system->getTimeStep();
        cout << "Thermo on" << endl;
    }
}
