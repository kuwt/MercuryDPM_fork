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
    double radius = 115e-6;
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
    double ambientTemperature = 428.15;
    double surfaceTension = 34.4e-3;
    double arrheniusCoefficient = 0.01;
    double appActivationEnergy = 32.5;
    double mass = density*4./3.*constants::pi*pow(radius,3);
    // other parameters
//    Mdouble maxSimTime = 4.0;
//    Mdouble particleInitialTemp = 155 + 273.15;
//    Mdouble damping = 0.9;
//    Mdouble refVis = 1e-4;
//    Mdouble aTemp = particleInitialTemp;
//    Mdouble mTemp = 178 + 273.15; // -> 451.15 (K)
//    Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
//    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
//    int thermalExpansionOnOff = 0;
//    Mdouble gravity = -9.81;
//    Mdouble cp_SL = (0.5*(solidHeatCapacity+liquidHeatCapacity)) + (latentHeat/deltaT);
//    Mdouble Ep = (mTemp-aTemp+(deltaT/2.0))*mass*cp_SL;//*(heatingTime/maxSimTime);
//    Mdouble pulseDuration = 1e-3;//1e-3;
//    Mdouble heatInputTUE = (384e-6)/pulseDuration; //(384e-6)/pulseDuration;
//    Mdouble heatSource = heatInputTUE;//(100e-6/pulseDuration) * scaleMass;
//    Mdouble heatingTime = 0.0;
//    Mdouble coolingTime = 2*pulseDuration; //1.0;//2*pulseDuration + heatingTime;
//    int saveCount = 1e4; //1e6; //2;//4;//1e8;
};

