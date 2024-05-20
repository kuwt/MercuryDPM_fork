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
#include "Species/MeltableFrictionSpecies.h"
//#include "Species/MeltableSpecies.h"
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Particles/BaseParticle.h"
//#include <Walls/InfiniteWall.h>
//#include "Interactions/NormalForceInteractions/MeltableInteraction.h"

class TestForceLaw : public Mercury3D
{
public:
    //
    //bool bonded = true;
    Mdouble simTime = 1.0; //1.5; //8.0;//0.5; //0.006; //0.005; //0.05; //0.5
    unsigned int fstatSetSaveCount =  10; //50000; //10000; //5000; //500 //1; //10;
    unsigned int saveCount = 1e5; //1e3; //1e7;//1e5;//1e4;
    bool vtkOut = true;
    //
    Mdouble pRadius1 = 20e-6; //20e-6; //50e-6; //11.5e-6; //100e-6;//11.5e-6; // 50-e-6
    Mdouble pRadius2 = 20e-6; //20e-6; //50e-6; //11.5e-6; //100e-6;//11.5e-6;
    // in case pR = 11.5e-6 -> timeStep = 0.3*getTimeStep
    //
    // OLD values? 3.45e-6 // R = 10e-6 -> 0.3e-6/ R = 50e-6 -> 1.5e-6
    // NEW: //0.85*pRadius2; //0.7*pRadius2; //0.03*pRadius2; //0.005*pRadius2;
    //- if overlap = 2R then particles centers are at the same point i.e. |rij| = 0
    Mdouble initialOverlap = 1.3*pRadius2;//1.6*pRadius2;//+(0.001*pRadius2); //0.03*pRadius2; //1.5*pRadius2;//1.8*pRadius2+(0.01*pRadius2);
    // 1.15e-6 //0.5e-6 -> solidOverlap=1.5e-6-1.0e-6=0.5e-6
    //
    // 0.9*pRadius1; //0.8*pRadius2; //0.01*pRadius2; //0.1*pRadius2;
    Mdouble moltenLaterThickness1 = 0.7*pRadius1;
    Mdouble moltenLaterThickness2 = 0.7*pRadius2;
    //delta_s = initialOverlap-(2.0*moltenLaterThickness); //0.0;  //2.3e-6; //1.0e-6;
    Mdouble initialBondingOverlap = 0.0;//1.4*pRadius2; //1.6*pRadius1;//0.0;
    //
    Mdouble aTemp = 155 + 273.15;//428.15 -- 443.15;//413.15;
    Mdouble mTemp = 178 + 273.15;//451.15; -- 451.15;
    Mdouble initialParticleTemperature = mTemp;
    //
    Mdouble MatDensity = 1050.0; // 1400
    Mdouble elasticModulus = 0.8*2.95e9; //8e4; //0.8*2.95e9; //8e7;
    Mdouble possionsRation = 0.4;//0.2;
    //
    Mdouble elasticTerm = (1.0-(mathsFunc::square(possionsRation)))/elasticModulus;
    Mdouble Eeff = 1.0/(elasticTerm+elasticTerm);
    //
    Mdouble damping = 0.95;//0.95; //0.1;//0.95;//0.5;//0.85;
    //
    //
    Mdouble latentHeat = 56.4e3;//245e3;
    Mdouble heatCapacity = 1200;//2270;
    Mdouble thermalcond = 0.0;//0.12;//0.36;
    Mdouble thermalConvCoff = 0.0;//100.0;//15;
    Mdouble emmisivity = 0.0;//0.9;
    //
    //Mdouble heatingTime = 0.0;//0.075;//0.0075;
    //Mdouble heatSource = 0.0;//2.5e-4;
    //Mdouble coolingTime = 0.0;//0.075*2.0;
    //
    //Mdouble SlidingFrictionCoeff = 0;//0.25;
    //Mdouble RollingFrictionCoeff = 0;//0.05;
    //
    Mdouble surfaceTension = 0.035;
    Mdouble refVis = 2e-3; //2e-3; //0.0;
    //Mdouble refTemp = 224;//122 + 272.15;
    //Mdouble tempCoeff = 0.053;
    //
    Mdouble deltaT = 20.0;
    //Mdouble liquidHeatCapacity = 2.0*heatCapacity;
    //Mdouble vaporizationHeatCapacity = 4.0*heatCapacity;
    //Mdouble vaporizationLatentHeat = 20.0*latentHeat;
    //Mdouble vaporizationTemp = 2.0*mTemp;
    //
    //
    Mdouble gravityAcc = 0.0; //0.0;//-9.81;
    //
    //-----------------------------------------------------------
    // TIME STEP criterion:
    //
/*    Mdouble Reff = (2.0*pRadius1*pRadius2)/(pRadius1+pRadius2);
    Mdouble volume = (4.0/3.0)*constants::pi*mathsFunc::cubic(pRadius2);
    Mdouble mass = MatDensity*volume;
    //
    Mdouble vMax = 1;
    //
    Mdouble C = 0.3;
    Mdouble timeStepYade =   C * pRadius2 *  sqrt(MatDensity/elasticModulus);
    Mdouble timeStepYade2 =  C * pRadius2 *  sqrt(MatDensity/Eeff);
    //
    Mdouble timeStepHertz = 2.87*pow(mathsFunc::square(mass)/(Reff*mathsFunc::square(Eeff)*vMax),0.2);
    //
    Mdouble timeStep = timeStepYade;*/
    //-----------------------------------------------------------
    //
    void setupInitialConditions() override {
        //
        setName("TestForceLaw");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0,gravityAcc));
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(0.0002);
        setYMax(0.0002);
        setZMax(0.0002);
        setTimeMax(simTime);
        //
        auto species = speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
        /*//
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        //
        p.setTemperature(initialParticleTemperature);
        p.setRadius(pRadius1); //particle-1 radius 50 um
        p.setPosition(Vec3D(pRadius1, 0.0, pRadius1));
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //p.setMoltenLayerThicknessRateParticle(-3e-6);
        p.addMoltenLayerThickness(moltenLaterThickness);
        //p.setMoltenLayerThickness(3e-5);
        auto p0 = particleHandler.copyAndAddObject(p); // 1st particle created
        //
        p.setTemperature(initialParticleTemperature);//);
        p.setRadius(pRadius2); // particle-2 radius
        p.setPosition(Vec3D(pRadius1, (pRadius1+pRadius2)-initialOverlap, pRadius1));//
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //p.setMoltenLayerThicknessRateParticle(-3e-6);
        p.addMoltenLayerThickness(moltenLaterThickness);
        //p.setMoltenLayerThickness(3e-5);
        auto p1 = particleHandler.copyAndAddObject(p); // 4nd particle created
        //
        // Calculate forces and set initial bonding overlap:
        auto i = interactionHandler.getInteraction(p0,p1,1);
        auto mi = dynamic_cast<MeltableInteraction*>(i);
        logger.assert_always(mi!= nullptr,"mi not set");
        mi->computeNormalForce();
        mi->setBondingOverlap(initialBondingOverlap);
        //*/
        /*std::cout << " Overlap " << mi->getOverlap() << " "
        << " Force " << mi->getAbsoluteNormalForce() << " "
        << " bondingOverlap " << mi->getBondingOverlap() << " "
        << " solidOverlap " << mi->getSolidOverlap() << std::endl;*/
        //
        species->setDensity(MatDensity);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(possionsRation);
        species->setDissipation(damping);
        species->setAmbientTemperature(aTemp);
        species->setMeltingTemperature(mTemp);
        //
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(heatCapacity);
        species->setThermalConductivityCoefficient(thermalcond);
        species->setThermalConvectionCoefficient(thermalConvCoff);
        species->setMaterialEmissivity(emmisivity);
        //species->setHeatingTime(heatingTime);
        //species->setHeatInput(heatSource);
        //species->setCoolingTime(coolingTime);
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(refVis);
        //species->setRefTemperature(refTemp);
        //species->setTemperatureSensitivityCoeff(tempCoeff);
        //
        species->setDeltaT(deltaT);
/*        species->setLiquidHeatCapacity(liquidHeatCapacity);
        species->setVaporizationHeatCapacity(vaporizationHeatCapacity);
        species->setVaporizationLatentHeat(vaporizationLatentHeat);
        species->setVaporizationTemperature(vaporizationTemp);*/
        //
        species->setThermalExpansionOnOff(0);
        //
        species->setMaxTimeStep(pRadius2,MatDensity,Eeff);//elasticModulus);
        setTimeStep(0.3*species->getMaxTimeStep()); //6e-10);//1e-9); //getAdaptiveTimeStep());//timeStep); // 4e-9 //1e-8;
        // species->getMaxTimeStep());// 0.3*species->getMaxTimeStep());
        //
/*  species->setSlidingStiffness(2./7.*species->getEffectiveElasticModulus());
    species->setSlidingDissipation(2./7.*species->getViscosityCoefficient());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getElasticModulus());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);*/
        //
        setSaveCount(saveCount);//5000);
        //
        //restartFile.setSaveCount(10000);
        //restartFile.setFileType(FileType::MULTIPLE_FILES);
        //
        fStatFile.setSaveCount(fstatSetSaveCount);//1);
        //eneFile.setSaveCount(1);
        //fStatFile.setFileType(FileType::MULTIPLE_FILES);
        dataFile.setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::ONE_FILE);//NO_FILE);
        //
        setParticlesWriteVTK(vtkOut);
        //setWallsWriteVTK(FileType::MULTIPLE_FILES);
        //
        //
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        //
        p.setTemperature(initialParticleTemperature);
        p.setRadius(pRadius1); //particle-1 radius 50 um
        p.setPosition(Vec3D(pRadius1, 0.0, pRadius1));
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //p.setMoltenLayerThicknessRateParticle(-3e-6);
        p.addMoltenLayerThickness(moltenLaterThickness1);
        //p.setMoltenLayerThickness(3e-5);
        auto p0 = particleHandler.copyAndAddObject(p); // 1st particle created
        //
/*        p.setTemperature(initialParticleTemperature);//);
        p.setRadius(pRadius2); // particle-2 radius
        p.setPosition(Vec3D(pRadius1, (pRadius1+pRadius2)-initialOverlap, pRadius1));//
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //p.setMoltenLayerThicknessRateParticle(-3e-6);
        p.addMoltenLayerThickness(moltenLaterThickness);
        //p.setMoltenLayerThickness(3e-5);
        auto p1 = particleHandler.copyAndAddObject(p); // 2nd particle created*/
        //--------
        ThermalParticle I;
        I.setSpecies(speciesHandler.getObject(0));
        I.setTemperature(initialParticleTemperature);//);
        I.setRadius(pRadius2); // particle-2 radius
        I.setPosition(Vec3D(pRadius1, (pRadius1+pRadius2)-initialOverlap, pRadius1));//
        I.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //p.setMoltenLayerThicknessRateParticle(-3e-6);
        I.addMoltenLayerThickness(moltenLaterThickness2);
        //p.setMoltenLayerThickness(3e-5);
        auto p1 = particleHandler.copyAndAddObject(I); // 2nd particle created
        //------------
        // Calculate forces and set initial bonding overlap:
        auto i = interactionHandler.getInteraction(p0,p1,1);
        auto mi = dynamic_cast<MeltableInteraction*>(i);
        logger.assert_always(mi!= nullptr,"mi not set");
        mi->computeNormalForce();
        //mi->setBondingOverlap(initialBondingOverlap);
        mi->addBondingOverlap(initialBondingOverlap);
        //
    }
/*    void writeFstatHeader(std::ostream& os) const override
    {
        for (const auto &q: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(q);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            double force = i0->getForce().y()/i0->getNormal().y();
            os << force << " " << i0->getOverlap() << " " << i0->getSolidOverlap() << " " << i0->getBondingOverlap() << std::endl;
        }
    }*/

    void printTime() const override {
        for (const auto& q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            logger(INFO, "time % temperature % moltenLayer % meltRate %",
                   getTime(), p0->getTemperature(), p0->getMoltenLayerThickness(), p0->getMoltenLayerThicknessRateParticle());
        }
        //
        for (const auto &q: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(q);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            Mdouble forceX = i0->getForce().x()*i0->getNormal().x();
            Mdouble forceY = i0->getForce().y()*i0->getNormal().y();
            Mdouble forceZ = i0->getForce().z()*i0->getNormal().z();
            logger(INFO, "time % forceX % forceY % forceZ % overlap % solidOverlap % bondingOverlap % meltRateP % meltRateI % "
                         "relativeVelocity % contactRadius % neck % "
                         "solidContactRadius % bondingContactRadius % "
                         "tensionForce % visForce % bondingForce % dampingForce % elasticForce % distance %",
                         getTime(), forceX, forceY, forceZ, i0->getOverlap(), i0->getSolidOverlap(), i0->getBondingOverlap(),
                         i0->getMoltenLayerThicknessRateP(), i0->getMoltenLayerThicknessRateI(),
                         i0->getNormalRelativeVelocity(),
                    i0->getContactRadiusMeltable(), i0->getContactRadiusMeltable()/pRadius1, i0->getSolidContactRadius(), i0->getBondingContactRadius(),
                    i0->getTensionForce(), i0->getViscousForce(), i0->getBondingForce(), i0->getDampingForce(), i0->getElasticForce(),
                    i0->getDistance());
        }
    }
    //
/*    Mdouble getAdaptiveTimeStep()
    {
        for (auto i : interactionHandler) {
            Mdouble vn = i->getNormalRelativeVelocity();
            adaptiveTimeStep_ = 2.87*pow(mathsFunc::square(mass)/(Reff*mathsFunc::square(Eeff)*vn),0.2);
        }
        if (adaptiveTimeStep_ <= 0)
            return 2.87*pow(mathsFunc::square(mass)/(Reff*mathsFunc::square(Eeff)*1.0),0.2);
        else
            return adaptiveTimeStep_;
    }
private:
    Mdouble adaptiveTimeStep_{};*/
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    TestForceLaw problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();

    return 0;
}
