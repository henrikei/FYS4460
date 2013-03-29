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
    int nAtomsPerDim = 12;
    double mass = 39.948*1.66E-27;
    double fccLength = 5.720E-10;
    double temperature = 0.851*119.74;
    double endTime = 2.5E-11;
    double ThermoTurnOffTime = 1.0E-11;
    double timeStep = 2.5E-14;// 4.862E-14;
    string thermostat = "berendsenThermostat";
    double relaxationTime = 10; // In units of timeSteps
    //vec3 externalForce;
    //externalForce << 0 << 0 << 0.1*1.6531E-21/3.405E-10;

    System test("Ar", nAtomsPerDim, mass, fccLength, temperature, endTime, timeStep);
    //test.addExternalForce(externalForce);

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
    //test.readState("circularPores0.xyz");
    test.integrate();
    return 0;
}
