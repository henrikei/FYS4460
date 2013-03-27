#include "andersenthermostat.h"
#include "system.h"
#include "atom.h"
#include "armadillo"

using namespace arma;

AndersenThermostat::AndersenThermostat(System *system) :
    Modifier(system)
{
    tau = 20;
}

void AndersenThermostat::setTargetTemperature(double T){
    targetTemperature = T/119.74;
}

void AndersenThermostat::setRelaxationTime(double relaxationTime){
    tau = relaxationTime;
}

void AndersenThermostat::setTurnOffTime(double t){
    tOff = t/(2.1569E-12);
}

void AndersenThermostat::apply(){
    if (time < tOff){
        atoms = system->getAtoms();
        nAtoms = system->getNumberOfAtoms();
        for (int i = 0; i < nAtoms; i++){
            if(atoms.at(i)->getFree()){
                vec3 newVelocity = randn(3)*sqrt(targetTemperature/atoms.at(i)->getMass());
                double randomCheck = randu();
                if (randomCheck < 1/tau){
                    atoms.at(i)->setVelocity(newVelocity);
                }
            }
        }
        time += system->getTimeStep();
        cout << "Thermo on" << endl;
    }
}
