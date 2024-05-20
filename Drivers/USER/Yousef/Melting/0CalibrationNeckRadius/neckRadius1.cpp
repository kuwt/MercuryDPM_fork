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
#include "Species/MeltableSpecies.h"
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include <Walls/InfiniteWall.h>


//TU/e PA12 Experiment: paper: -- ADD HERE --

class neckRadius1: public Mercury3D
{
private:
    //
    Mdouble timeStep = 1e-7;
    Mdouble maxSimTime = 6.0;
    // material properties and simulation parameters:
    Mdouble pRadius1 = 115e-6; //20e-6;
    Mdouble pRadius2 = 115e-6; //20e-6;
    //
    Mdouble minRadius = pRadius2;
    Mdouble radiusForHeatInput = pRadius2;
    // --------------------------------------------------------------------
    Mdouble particleInitialTemp = 155 + 273.15; // -> 428.15 (K)
    //
    Mdouble scaleMass = 1; //1e4;//1e6; //1e8;
    //
    Mdouble MatDensity = 1050 * scaleMass;
    Mdouble elasticModulus = 0.8*2.95e9;
    Mdouble possionsRation = 0.4;
    Mdouble damping = 0.9; //0.9; //0.9; //0.95;
    //
    Mdouble latentHeat = 56.4e3;
    Mdouble heatCapacity = 1200;
    //
    Mdouble scaleTension1 = 1e-3;
    Mdouble surfaceTension = 34.3e-3;// * scaleTension;// * 1e-3; // * 2e-2 *  sqrt(scaleMass); //2e-2 //2e-3  //3e-2;
    //Mdouble refVis = 1e-4; //2 * 2e-3 * sqrt(scaleMass);; //3e-1; //3e-1; // 2 * 2e-3 * sqrt(scaleMass); //4e-3 //2e-2 // 2e-3
    // cond , conv, rad
    Mdouble thermalcond = 0.12 * scaleMass;
    Mdouble thermalConvCoff = 0.0;//40.0 * scaleMass; //100
    Mdouble emmisivity = 0.0;//0.9 * scaleMass;
    //
    Mdouble aTemp = particleInitialTemp;
    Mdouble mTemp = 178 + 273.15; // -> 451.15 (K)
    Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
    Mdouble liquidHeatCapacity = 1.2 * heatCapacity;//<- assumed, 2.0*heatCapacity;
    Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//<- assumed, 4.0*heatCapacity;
    Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
    Mdouble vaporizationTemp = 3.5 * mTemp;
    //
    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
    int thermalExpansionOnOff = 0;
    //
    Mdouble gravityAcc = -9.81/sqrt(scaleMass); //scaleMass
    //------------------------------------------------------------------------------------------------
    Mdouble Lm = (0.5 * deltaT * (heatCapacity + liquidHeatCapacity)) + latentHeat;
    Mdouble mass = MatDensity*((4.0/3.0)*constants::pi*mathsFunc::cubic(radiusForHeatInput));
    // -1-
    // Heat input based on melt layer def:
    //Mdouble heatSource = mass*Lm; //(1.0/1e-3)*mass*Lm; //(1.0/heatingTime)*
    // -2-
    //Mdouble meltLimit = 0.99;
    //Mdouble solidR = radiusForHeatInput*(1-meltLimit);
    //Mdouble heatSource = (1.0/maxSimTime)*MatDensity*(4.0*constants::pi*meltLimit*radiusForHeatInput*mathsFunc::square(solidR))*Lm;
    // -3-
    Mdouble cp_SL = (0.5*(heatCapacity+liquidHeatCapacity)) + (latentHeat/deltaT);
    Mdouble Ep = (mTemp-aTemp+(deltaT/2.0))*mass*cp_SL;//*(heatingTime/maxSimTime);
    //Mdouble heatSource = 0.0;//Ep;
    ////------------------------------------
    Mdouble pulseDuration = 1e-3;//1e-3;
    Mdouble heatInputTUE = (384e-6)/pulseDuration; //(384e-6)/pulseDuration;
    Mdouble heatSource = heatInputTUE * scaleMass;//(100e-6/pulseDuration) * scaleMass;
    ////------------------------------------
    Mdouble heatingTime = 0.0;
    Mdouble coolingTime = 2*pulseDuration; //1.0;//2*pulseDuration + heatingTime;
    //------------------------------------------------------------------------------------------------
    int saveCount = 1500000;//40/1e-7;

public:

    //void setupInitialConditions() override {}
    //
    explicit neckRadius1(Mdouble key1, Mdouble refVis, Mdouble scaleTension)
    {

        std::string key = helpers::toString(key1);
        std::string refViscosity = helpers::toString(refVis);
        std::string scaleTensionFactor = helpers::toString(scaleTension);

        setName("GLNeck_"+key+"_"+refViscosity+"_"+ scaleTensionFactor);


        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, gravityAcc));
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(0.0002);
        setYMax(0.0002);
        setZMax(0.0002);
        setTimeMax(maxSimTime);
        setHGridMaxLevels(3);
        //
        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        //
        species->setDensity(MatDensity);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(possionsRation);
        species->setDissipation(damping);
        //
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(heatCapacity);
        species->setThermalConductivityCoefficient(thermalcond);
        species->setThermalConvectionCoefficient(thermalConvCoff);
        species->setAmbientTemperature(aTemp);
        species->setMeltingTemperature(mTemp);
        species->setMaterialEmissivity(emmisivity);
        //
        species->setHeatingTime(heatingTime);
        species->setHeatInput(heatSource);
        species->setCoolingTime(coolingTime);
        //
        species->setSurfaceTension(surfaceTension * std::abs(scaleTension));
        species->setRefViscosity(std::abs(refVis));
        //

        species->setDeltaT(deltaT);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        species->setVaporizationHeatCapacity(vaporizationHeatCapacity);
        species->setVaporizationLatentHeat(vaporizationLatentHeat);
        species->setVaporizationTemperature(vaporizationTemp);
        //
        species->setThermalExpansionCoefficient(thermalExpansionCoeff);
        species->setThermalExpansionOnOff(thermalExpansionOnOff);


        setTimeStep(timeStep);

        //
        setSaveCount(saveCount);

        dataFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::NO_FILE);

        //setParticlesWriteVTK(true);
        ////setWallsWriteVTK(FileType::MULTIPLE_FILES);

        // Sim setup:
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set({0.0, 0.0, -1.0}, {0.0, 0.0, 0.0});
        wallHandler.copyAndAddObject(w0);

        //
        ThermalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setTemperature(particleInitialTemp);
        p0.setRadius(pRadius1); //particle-1 radius 50 um
        p0.setPosition(Vec3D(pRadius1, 0.0, pRadius1));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0); // 1st particle created
        //
        ThermalParticle p1;
        p1.setSpecies(speciesHandler.getObject(0));
        p1.setTemperature(particleInitialTemp);
        p1.setRadius(pRadius2); // particle-2 radius
        // as averaged exp:
        p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.01*pRadius1),pRadius2));// 2-> (0.01*pRadius1) // 1-> (0.05*pRadius1)
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p1); // 2nd particle created
    }

    //
    void writeFstatHeader(std::ostream& os) const override
    {
        //write interaction properties at every time step
        auto i = interactionHandler.getObject(0); //getLastObject(); // getObject(0);
        logger.assert_debug(i,"No interaction exists");
        auto mi = dynamic_cast<const MeltableInteraction *>(i);
        logger.assert_debug(mi,"Interaction is not of type MeltableInteraction");

        // header
        if (getTime() == 0){
            os << '#' <<" " << "time" << " " << "neck" <<'\n';
        }

        os      <<  getTime()
                << " " << mi->getContactRadiusMeltable()/pRadius1
                << '\n';

    }


    void printTime() const override
    {
        auto i = interactionHandler.getObject(0); //getLastObject(); // getObject(0);
        logger.assert_debug(i,"No interaction exists");
        auto mi = dynamic_cast<const MeltableInteraction *>(i);
        logger.assert_debug(mi,"Interaction is not of type MeltableInteraction");

        logger(INFO,"Time % TimeMax % neck %", getTime(), getTimeMax(), mi->getContactRadiusMeltable()/pRadius1);
    }
};


int main(int argc , char* argv[] )
{
    Mdouble key = helpers::readFromCommandLine(argc,argv,"-key",0);
    Mdouble calibrationRefVis =  helpers::readFromCommandLine(argc,argv,"-calibrationRefVis",1e-4);
    Mdouble calibrationTension =  helpers::readFromCommandLine(argc,argv,"-calibrationTension",1e-3);

    neckRadius1 problem(key,calibrationRefVis,calibrationTension);

    problem.solve();

    return 0;
}
