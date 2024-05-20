//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "Mercury3D.h"
#include "cmath"

/**
 * Set parameters for the contact model paper
 */
class Initialization : public Mercury3D {
protected:
    double radius = 50e-6;
    double density = 1050;
    double elasticModulus = 23.6e8;
    double poissonRatio = 0.4;
    double latentHeat = 56.4e3;// /(1-pow(1-relMoltenLayerThicknessMax,3));
    double solidHeatCapacity = 1200;
    double liquidHeatCapacity = 1.2*solidHeatCapacity;
    double temperatureInterval = 20;
    double meltingTemperature = 451.15;
    double thermalConductivity = 0.12;
    double thermalConvection = 150;
    double emmisivity = 0.9;
    double ambientTemperature = meltingTemperature-2.0*temperatureInterval;
    double surfaceTension = 34.4e-3;
    double activationEnergy = 32.5e3;
    double mass = density*4./3.*constants::pi*pow(radius,3);
    double dissipation = 0.9;
    double referenceViscosity = 1e-2;
    double gravity = -9.81;
    double timeStep = 19e-9;

    void setScaleFactor (double scaleFactor) {
        this->scaleFactor = scaleFactor;
        density *= scaleFactor;
        surfaceTension *= sqrt(scaleFactor); //0.035 * 1e-3 * sqrt(scaleMass);// * 2e-3;
        referenceViscosity *= sqrt(scaleFactor); //2e-4 * sqrt(scaleMass);
        gravity /= scaleFactor;
        solidHeatCapacity /= scaleFactor;
        liquidHeatCapacity /= scaleFactor;
        timeStep *= sqrt(scaleFactor);
    }

    double getScaleFactor() {
        return scaleFactor;
    }

private:
    double scaleFactor = 1.0;
};

