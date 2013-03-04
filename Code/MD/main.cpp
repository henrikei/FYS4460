#include <iostream>
#include <armadillo>
#include <sstream>
#include <cmath>
#include "atom.h"
#include "system.h"
#include "Modifier/modifier.h"
#include "Modifier/berendsenthermostat.h"
#include "Modifier/andersenthermostat.h"

using namespace std;
using namespace arma;

int main()
{
    // Input in SI units
    string atomType = "Ar";
    int nAtomsPerDim = 8;
    double mass = 39.948*1.66E-27;
    double fccLength = 5.260E-10;
    double temperature = 50;
    double endTime = 4.0E-12;
    double ThermoTurnOffTime = 2.0E-12;
    double timeStep = 1.0E-14;// 4.862E-14;
    string thermostat = "berendsenThermostat";
    double relaxationTime = 10; // In units of timeSteps

    System test("Ar", nAtomsPerDim, mass, fccLength, temperature, endTime, timeStep);

    if (thermostat == "berendsenThermostat"){
        BerendsenThermostat *thermos = new BerendsenThermostat(& test);
        thermos->setTargetTemperature(temperature);
        thermos->setRelaxationTime(relaxationTime);
        thermos->setTurnOffTime(ThermoTurnOffTime);
        test.addModifier(thermos);
    }

    if (thermostat == "andersenThermostat"){
        AndersenThermostat *thermos = new AndersenThermostat(& test);
        thermos->setTargetTemperature(temperature);
        thermos->setRelaxationTime(relaxationTime);
        thermos->setTurnOffTime(ThermoTurnOffTime);
        test.addModifier(thermos);
    }

    test.generate();
    test.integrate();
    return 0;
}
