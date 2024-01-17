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
#include <Walls/InfiniteWall.h>
//#include "Interactions/NormalForceInteractions/MeltableInteraction.h"
#include "../Initialization.h"

// Questions:
// what happen when heating to e.g. T=190 then stop the heat input, do particles continue to coalescence?

class Test2ParticlesBonding: public Mercury3D
{
public:
    //
    Mdouble timeStep = 1e-7;
    //
    Mdouble maxSimTime = 9.0; //9.0;//6.0;
    // material properties and simulation parameters:
    Mdouble pRadius1 = 37e-6; //R4,6: 31e-6; //33e-6;
    Mdouble pRadius2 = 37e-6; //R4,6: 31e-6; //33e-6;
    //
    Mdouble minRadius = pRadius2;//45e-6;//pRadius2;
    Mdouble radiusForHeatInput = pRadius1;//pRadius2;//115e-6;//pRadius1;
    //
    //
    Mdouble particleInitialTemp = 170 + 273.15; //160 + 273.15; //175 + 273.15;
    //
    Mdouble MatDensity = 1050;//1400.0;
    Mdouble elasticModulus =  0.8*2.95e9;
    Mdouble possionsRation = 0.4;//0.2;
    Mdouble damping = 0.9;//0.95;
    //
    Mdouble latentHeat = 102e3;//56.4e3;//245e3;
    Mdouble heatCapacity = 2900;//1200;//2270;
    //
    Mdouble surfaceTension = 34.3e-3 * 1e-5; // 8e-3; //0.035;
    Mdouble refVis = 2e-5; //2e-4;//2e-3;
    // cond , conv, rad
    Mdouble thermalcond = 0.12;
    Mdouble thermalConvCoff = 0.0;//150;//300.0;//15;//25.0;// 10.0; 75.0;//150.0;
    Mdouble emmisivity = 0.0;//0.9;//0.1;//0.7;//0.9;
    //
    Mdouble aTemp = particleInitialTemp;//155 + 273.15;//428.15
    Mdouble mTemp = 185 + 273.15;//178 + 273.15;//451.15;
    Mdouble deltaT = 20; //10; //20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
    Mdouble liquidHeatCapacity = 3400;//1.2 * heatCapacity;//2.0*heatCapacity;
    Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//4.0*heatCapacity;
    Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
    Mdouble vaporizationTemp = 3.5 * mTemp;
    //
    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
    int thermalExpansionOnOff = 0;
    //------------------------------------------------------------------------------------------------
    Mdouble Lm = (0.5 * deltaT * (heatCapacity + liquidHeatCapacity)) + latentHeat;
    Mdouble mass = MatDensity*((4.0/3.0)*constants::pi*mathsFunc::cubic(radiusForHeatInput));
    //
    Mdouble heatingTime = 0.0;
    // heat source in watt
    //Mdouble heatSource = mass*Lm;//(1.0/heatingTime)*
    //
    //Mdouble meltLimit = 0.99;
    //Mdouble solidR = radiusForHeatInput*(1-meltLimit);
    //Mdouble heatSource = (1.0/maxSimTime)*MatDensity*(4.0*constants::pi*meltLimit*radiusForHeatInput*mathsFunc::square(solidR))*Lm;
    //
    Mdouble cp_SL = (0.5*(heatCapacity+liquidHeatCapacity)) + (latentHeat/deltaT);
    Mdouble Ep = (mTemp-aTemp+(deltaT/2.0))*mass*cp_SL;//*(heatingTime/maxSimTime);
    //Mdouble heatSource = Ep/6.0; //6.0;// /heatingTime;//(mTemp-aTemp+(deltaT/2.0))*mass*cp_SL*(1/maxSimTime);//(heatingTime/maxSimTime);
    // averaged axp:
    Mdouble heatSource = Ep/1.0; //6.0;// /heatingTime;//(mTemp-aTemp+(deltaT/2.0))*mass*cp_SL*(1/maxSimTime);//(heatingTime/maxSimTime);
    //
    //Mdouble pulseDuration = 1e-3;//1e-3;
    Mdouble heatingRate = 2+273.15;
    //Mdouble heatSource = 0.0;//heatingRate*mass*heatCapacity;//getHeatInput();//(384e-6/2)/pulseDuration;
    //
    Mdouble coolingTime = 2.0; // //1,2-> 2.0;         //6.0;
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
    //
    void setupInitialConditions() override {
        //
        //
        setName("TwoParticlesHSM");
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
        auto species = speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
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
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(refVis);
        //
        species->setDeltaT(deltaT);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        species->setVaporizationHeatCapacity(vaporizationHeatCapacity);
        species->setVaporizationLatentHeat(vaporizationLatentHeat);
        species->setVaporizationTemperature(vaporizationTemp);
        //
        species->setThermalExpansionCoefficient(thermalExpansionCoeff);
        species->setThermalExpansionOnOff(thermalExpansionOnOff);
        //
        species->setMaxTimeStep(minRadius,MatDensity,elasticModulus);
        //
/*        species->setSlidingStiffness(2./7.*species->getElasticModulus());
        species->setSlidingDissipation(2./7.*species->getDissipation());
        species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
        //
        species->setRollingStiffness(2./5.*species->getElasticModulus());
        species->setRollingDissipation(2./5.*species->getDissipation());
        species->setRollingFrictionCoefficient(RollingFrictionCoeff);*/
        //
        setTimeStep(timeStep); //0.3*
        //
        setSaveCount(saveCount);
        //restartFile.setSaveCount(restartSaveCount);
        //fStatFile.setSaveCount(10000);
        dataFile.setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        //
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        //
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
        //-------------------------------------------
/*        // top of each other:
        //p1.setPosition(Vec3D(pRadius1, 0.0, (pRadius1*2.0 + pRadius2)));
        // same plan i.e. next to eaxh other: pRadius1 = pRadius2
        //p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0),pRadius2));
        // same plan with initial overlap
       //p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.02*pRadius1),pRadius2));//(0.03*pRadius1)
        // R 33e-6:
        p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.005*pRadius1),pRadius2));//(0.03*pRadius1)
        // smaller particle - same plan with initial overlap
        // 1
        //p1.setPosition(Vec3D(pRadius1,(pRadius1/2.0)+(2*pRadius2)-(0.05*pRadius1),pRadius2));
        // 2
        //p1.setPosition(Vec3D(pRadius1,(pRadius1/2.0)+(2*pRadius2)-(0.02*pRadius1),pRadius2));*/
        //-------------------------------------------
        // as averaged exp:
        p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.01*pRadius1),pRadius2));// 2-> (0.01*pRadius1) // 1-> (0.05*pRadius1)
        //-------------------------------------------
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p1); // 2nd particle created

/*        p.setTemperature(443.15);
        p.setRadius(0.000025); // particle-3 radius
        p.setPosition(Vec3D(0.00005, 0.0, (0.00005*4.0)+p.getRadius()));
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p); // 3nd particle created*/

/*        p.setTemperature(451.15);//);
        p.setRadius(0.00005); // particle-2 radius
        //p.setPosition(Vec3D(p.getRadius()*3.0, 0.0, p.getRadius()));
        p.setPosition(Vec3D(p.getRadius(), p.getRadius()*2.0, p.getRadius()));//
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        auto p1 = particleHandler.copyAndAddObject(p); // 4nd particle created*/
        //
/*        auto i = interactionHandler.getInteraction(p0,p1,getTime());
        auto mi = dynamic_cast<MeltableDPMInteraction*>(i);
        logger.assert_always(mi!= nullptr,"mi not set");
        mi->computeNormalForce();
        mi->setBondingOverlap(6e-6);
        //
        std::cout << " bondingOverlap " << mi->getBondingOverlap() << std::endl;
        std::cout << " solidOverlap " << mi->getSolidOverlap() << std::endl;*/
    }
    //
    void actionsAfterTimeStep() override {
        Mdouble heatingRate = 2;// C/s
        //
/*        for (const auto &q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
           //
           p0->addTemperature(heatingRate*getTimeStep());
        }*/
        //
        /*for (const auto &q: particleHandler)
        {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            //
            if (p0->getTemperature() < mTemp - (deltaT / 2.0))
            {
                p0->setLocalHeatInput(heatingRate*p0->getMass()*heatCapacity);
            }
            else if (mTemp - (deltaT / 2.0) <= p0->getTemperature()
                     & p0->getTemperature() <= mTemp + (deltaT / 2.0))
            {
                Mdouble heatCapacityL = 0.5*(heatCapacity+liquidHeatCapacity)+(latentHeat/deltaT);
                p0->setLocalHeatInput(heatingRate*p0->getMass()*heatCapacityL);
            }
            else if (mTemp + (deltaT / 2.0) < p0->getTemperature()
                     & p0->getTemperature() < vaporizationTemp - (deltaT / 2.0))
            {
                p0->setLocalHeatInput(heatingRate*p0->getMass()*liquidHeatCapacity);
            }
        }*/
    }
    //
    void writeFstatHeader(std::ostream& os) const override
    {
        for (const auto &q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            os << p0->getTemperature() << " " << p0->getMoltenLayerThickness()
               << " " << 00 << " " << 00 << " " << 00 <<  " " << 00 << " " << 00
               << " " << 00 << " " << 00 << " " << 00 << " " << 00 << " " << 00
               << " " << 00 << " " << 00 << std::endl;
        }
        //
        for (const auto &iH: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(iH);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            //double force = i0->getForce().z()/i0->getNormal().z();
            //Mdouble force = i0->getAbsoluteNormalForce();
            os << i0->getForce() << " " << i0->getOverlap()
               << " " << i0->getSolidOverlap() << " " << i0->getBondingOverlap() << " " << i0->getContactRadiusMeltable()
               << " " << i0->getElasticForce() << " " << i0->getDampingForce()
               << " " << i0->getBondingForce() << " " << i0->getViscousForce() << " " << i0->getTensionForce()
               << " " << i0->getSolidContactRadius() << " " << i0->getBondingContactRadius() << std::endl;
        }

    }
    //
/*    void printTime() const override
    {
        for (const auto& q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            logger(INFO, "time % temperature % moltenLayer %", getTime(), p0->getTemperature(), p0->getMoltenLayerThickness());
        }
        //
        for (const auto &iH: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(iH);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            //double force = i0->getForce().y()/i0->getNormal().y();
            logger(INFO, "time % force % overlap % solidOverlap % bondingOverlap % contactRadius % "
                         "solidContactRadius % bondingContactRadius % normalRelativeVelocity %", getTime(), i0->getForce(), i0->getOverlap(),
                   i0->getSolidOverlap(), i0->getBondingOverlap(), i0->getContactRadiusMeltable(),
                   i0->getSolidContactRadius(), i0->getBondingContactRadius(),i0->getNormalRelativeVelocity());
        }
    }*/
    //
/*    Mdouble getHeatInput()
    {
        Mdouble heatingRate = 2;// C/s
        for (const auto &q: particleHandler)
        {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            //
            if (p0->getTemperature() < mTemp - (deltaT / 2.0))
            {
                heatInput_ =  heatingRate*p0->getMass()*heatCapacity;
            }
            else if (mTemp - (deltaT / 2.0) <= p0->getTemperature()
                     & p0->getTemperature() <= mTemp + (deltaT / 2.0))
            {
                Mdouble heatCapacityL = 0.5*(heatCapacity+liquidHeatCapacity)+(latentHeat/deltaT);
                heatInput_ = heatingRate*p0->getMass()*heatCapacityL;
            }
            else if (mTemp + (deltaT / 2.0) < p0->getTemperature()
                     & p0->getTemperature() < vaporizationTemp - (deltaT / 2.0))
            {
                heatInput_ = heatingRate*p0->getMass()*liquidHeatCapacity;
            }
        }
        //
        return heatInput_;
    }
private:
    Mdouble heatInput_;*/
};


int main(int argc UNUSED, char* argv[] UNUSED)
{
    Test2ParticlesBonding problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();
    return 0;
}
