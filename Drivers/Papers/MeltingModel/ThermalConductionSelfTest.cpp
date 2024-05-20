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

#include "Species/MeltableSpecies.h"
#include "Particles/MeltableParticle.h"
#include "Initialization.h"
using constants::pi;
using std::pow;

class ThermalConductionSelfTest : public Initialization
{
private:
    bool writeVTK = false;
    Mdouble temperature1 = meltingTemperature-temperatureInterval;
    Mdouble temperature2 = meltingTemperature-2.0*temperatureInterval;
    double timeStep = 1e-7;
    double timeMax = 2;

public:

    ThermalConductionSelfTest(int argc, char* argv[]) {
        // reduce time for selftest
        if (not helpers::readFromCommandLine(argc,argv,"-full")) {
            timeMax = 0.1;
        }
    }

    void setupInitialConditions() override {
        setName("ThermalConductionSelfTest");
        setDomain({-radius, -radius, -radius},
                  {radius, radius, radius});
        setTimeMax(timeMax);

        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species->setDensity(density);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(poissonRatio);
        species->setDissipation(0.5);
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(solidHeatCapacity);
        species->setThermalConductivityCoefficient(thermalConductivity);
        //species->setThermalConvectionCoefficient(thermalConvection);
        //species->setAmbientTemperature(ambientTemperature);
        species->setMeltingTemperature(meltingTemperature);
        //species->setMaterialEmissivity(emmisivity);
        //species->setHeatingTime(heatingTime);
        //species->setCoolingTime(coolingTime);
        //species->setHeatInput(heatSource);
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(referenceViscosity);
        species->setDeltaT(temperatureInterval);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        //species->setThermalExpansionCoefficient(thermalExpansionCoeff);

        setTimeStep(timeStep);
        setSaveCount(timeMax/timeStep/50);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);

        double overlap = 0.1*radius;

        MeltableParticle p0;
        p0.setSpecies(species);
        p0.setTemperature(temperature1);
        p0.setRadius(radius);
        p0.setPosition(Vec3D(0.0, 0.0,-radius+0.5*overlap));
        particleHandler.copyAndAddObject(p0);
        //
        MeltableParticle p1;
        p1.setSpecies(species);
        p1.setTemperature(temperature2);
        p1.setRadius(radius);
        p1.setPosition(Vec3D(0.0, 0.0, radius-0.5*overlap));
        particleHandler.copyAndAddObject(p1);

        auto i = dynamic_cast<MeltableInteraction*>(interactionHandler.getInteraction(
                particleHandler.getObject(0),particleHandler.getObject(1),0));
        i->setBondingOverlap(overlap);
        i->addNeckRadius2(radius*overlap);
    }

    //write out data to plot
    void actionsBeforeTimeStep() override {
        if (getNumberOfTimeSteps()==20) {
            write(std::cout);
            for (const auto& i: interactionHandler) {
                logger(INFO,"i % %",i->getOverlap(), i->getForce());
            }
        }
        // get particle
        auto p1 = dynamic_cast<MeltableParticle *>(particleHandler.getObject(0));
        auto p2 = dynamic_cast<MeltableParticle *>(particleHandler.getObject(1));
        logger.assert_debug(p1 != nullptr or p2 != nullptr, "MeltableParticle required");
        // compute energies
        //write to file
        static std::ofstream file(getName()+".out");
        if (getNumberOfTimeSteps()%dataFile.getSaveCount()==0)
            file << getTime()
                    << ' ' << p1->getTemperature()
                    << ' ' << p1->getMoltenLayerThickness()
                    << ' ' << p2->getTemperature()
                    << ' ' << p2->getMoltenLayerThickness()
                    << ' ' << p1->getPosition().Z
                    << ' ' << p1->getVelocity().Z
                    << ' ' << p1->getForce().Z
                    << std::endl;
    }

    void printTime() const override
    {
        auto p = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(0));
        auto q = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(1));
        auto i = dynamic_cast<const MeltableInteraction*>(interactionHandler.getObject(0));
        logger(INFO, "time % temperatures % % overlap % neck %",
               getTime(), p->getTemperature(), q->getTemperature(), i->getOverlap(), i->getNeckRadius());
    }
};

int main(int argc, char* argv[])
{
    ThermalConductionSelfTest problem(argc, argv);
    problem.solve();
    return 0;
}
