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

class InitializationPA12 : public Mercury3D {
public:
    //
    Mdouble maxSimTime = 6.0;//6.0;
    // material properties and simulation parameters:
    Mdouble pRadius1 = 115e-6;//115e-6;//30e-6;//115e-6;//140e-6;//115e-6; //50e-6;
    Mdouble pRadius2 = 115e-6;//15e-6;//30e-6;//115e-6;//20e-6;//50e-6;//115e-6;//115e-6;//115e-6;//140e-6;//115e-6; //50e-6;
    //
    Mdouble minRadius = pRadius2; //0.4*pRadius2; //pRadius2;
    Mdouble radiusForHeatInput = pRadius2;//115e-6;//pRadius1;
    //
    Mdouble particleInitialTemp = 155 + 273.15; //20 + 273.15;//428.15
    //
    Mdouble MatDensity = 1050;//1400.0;
    Mdouble elasticModulus = 0.8*2.95e9;//1e9;//8e7; --- ///1e7; // J: 2.95e9;
    Mdouble possionsRation = 0.4;//0.2;
    Mdouble damping = 0.95;//0.9;//5e-6*10.0;
    //
    Mdouble latentHeat = 56.4e3;//245e3;
    Mdouble heatCapacity = 1200;//2270;
    //
    Mdouble surfaceTension = 34.3e-3;//0.035;
    Mdouble refVis = 2e-3;//2e-3;
    // cond , conv, rad
    Mdouble thermalcond = 0.12;
    Mdouble thermalConvCoff = 75.0;//50.0;//150;//300.0;//15;//25.0;// 10.0; 75.0;//150.0;
    Mdouble emmisivity = 0.9;//0.65;//0.9;//0.1;//0.7;//0.9;
    //
    Mdouble aTemp = 155 + 273.15;//428.15
    Mdouble mTemp = 178 + 273.15;//451.15;
    Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
    Mdouble liquidHeatCapacity = 1.2 * heatCapacity;//2.0*heatCapacity;
    Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//4.0*heatCapacity;
    Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
    Mdouble vaporizationTemp = 3.5 * mTemp;
    //
    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
    int thermalExpansionOnOff = 0;
    //------------------------------------------------------------------------------------------------
    Mdouble heatingTime = 0.0;//0.05;//0.0;//2*1e-3;//1e-3;//3.0;//3;//1e-3;//800e-3;
    //------------------------------------
    Mdouble Lm = (0.5 * deltaT * (heatCapacity + liquidHeatCapacity)) + latentHeat;
    Mdouble mass = MatDensity*((4.0/3.0)*constants::pi*mathsFunc::cubic(radiusForHeatInput));
    // heat source in watt
    // Heat input based on melt layer def:
    //Mdouble heatSource = (1.0/4.0)*mass*Lm;//(1.0/heatingTime)*
    ////------------------------------------
    //Mdouble meltLimit = 0.99;
    //Mdouble solidR = radiusForHeatInput*(1-meltLimit);
    //Mdouble heatSource = (1.0/maxSimTime)*MatDensity*(4.0*constants::pi*meltLimit*radiusForHeatInput*mathsFunc::square(solidR))*Lm;
    ////------------------------------------
    Mdouble cp_SL = (0.5*(heatCapacity+liquidHeatCapacity)) + (latentHeat/deltaT);
    Mdouble Ep = (mTemp-aTemp+(deltaT/2.0))*mass*cp_SL;//*(heatingTime/maxSimTime);
    //Mdouble heatSource = Ep;// /heatingTime;//(mTemp-aTemp+(deltaT/2.0))*mass*cp_SL*(1/maxSimTime);//(heatingTime/maxSimTime);
    ////------------------------------------
    Mdouble pulseDuration = 1e-3;//1e-3;
    Mdouble heatSource = (384e-6)/pulseDuration;//384e-6;// = 192e-6*2.0;///pulseDuration;//(384e-6)/1e-3;//100e-6/800e-3;//1.5;//1.5;//2.5;//2.5e-4;// -> used case in presentation;//
    //(50e-6)/pulseDuration;//
    ////------------------------------------
    Mdouble coolingTime = 2*pulseDuration + heatingTime;//2*1e-3;// 1*1e-3;//2*heatingTime;//0.0;//0.001*2*5;//0.075*2.0;
    //------------------------------------------------------------------------------------------------
    // LASER:
    //------------------------------------------------------------------------------------------------
    int saveCount = 1e6;//10000;//100000
    //int restartSaveCount = 20000;
    //
    Mdouble gravityAcc = -9.81;
    //
    Mdouble SlidingFrictionCoeff = 0.0;//0.5;
    Mdouble RollingFrictionCoeff = 0.0;//0.5;//0.05;
    //
};

