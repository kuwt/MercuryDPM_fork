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


//#include "Mercury3D.h"
//
//#include "Species/MeltableFrictionSpecies.h"
#include "Species/MeltableSpecies.h"
//
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Particles/BaseParticle.h"
#include <Walls/InfiniteWall.h>
//#include "Interactions/NormalForceInteractions/MeltableInteraction.h"
#include "InitializationPA12.h"
//#include "InitializationPS.h"


class TwoParticlesTUE: public InitializationPA12
{
public:
    //
    void setupInitialConditions() override {
        //
        //
        setName("TwoParticlesTUE");
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
        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());//MeltableFrictionSpecies());
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
        setTimeStep(0.3*species->getMaxTimeStep()); //species->getMaxTimeStep()
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
        // top of each other:
        //p1.setPosition(Vec3D(pRadius1, 0.0, (pRadius1*2.0 + pRadius2)));
        // same plan i.e. next to eaxh other: pRadius1 = pRadius2
        //p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0),pRadius2));
        // same plan with initial overlap
        // large
        p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.02*pRadius1),pRadius2));// (0.02*pRadius1) /(0.03*pRadius1)
        //small
        //p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(0.005*pRadius1),pRadius2));//(0.03*pRadius1)
        // --------- smaller particle with larger one - same plan with initial overlap
        // 1
        //p1.setPosition(Vec3D(pRadius1,(pRadius1/2.0)+(2*pRadius2)-(0.05*pRadius1),pRadius2));
        // 2
        //p1.setPosition(Vec3D(pRadius1,(pRadius1/2.0)+(2*pRadius2)-(0.02*pRadius1),pRadius2));
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
    void printTime() const override
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
    }
};


int main(int argc UNUSED, char* argv[] UNUSED)
{
    TwoParticlesTUE problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();
    return 0;
}
