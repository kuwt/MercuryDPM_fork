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

class MeltableForceLawSelfTest : public Initialization
{
private:
    bool writeVTK = false;
    Mdouble temperature1 = meltingTemperature-temperatureInterval;
    Mdouble temperature2 = meltingTemperature-2.0*temperatureInterval;
    double timeMax = 1000*timeStep;
    unsigned saveCount = 1;
    std::ofstream file;

public:
    MeltableForceLawSelfTest(double dissipation) {
        this->dissipation = dissipation;
        solve();
    }

    void setupInitialConditions() override {
        setName("MeltableForceLawSelfTest"+helpers::toString(dissipation));
        file.open(getName()+".out");
        logger(INFO,"name %",getName());
        setDomain({-radius, -radius, -radius},
                  {radius, radius, radius});
        setTimeMax(timeMax);

        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species->setDensity(density);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(poissonRatio);
        species->setDissipation(dissipation);
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
        species->setActivationEnergy(activationEnergy);
        //species->setThermalExpansionCoefficient(thermalExpansionCoeff);

        setTimeStep(timeStep);
        setSaveCount(saveCount);
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
        i->setBondingOverlap(.55*overlap);
    }

    //write out data to plot
    void actionsBeforeTimeStep() override {
//        if (getNumberOfTimeSteps()==20) {
//            write(std::cout);
//            for (const auto& i: interactionHandler) {
//                logger(INFO,"i % %",i->getOverlap(), i->getForce());
//            }
//        }
        // get particle
        auto p1 = dynamic_cast<MeltableParticle *>(particleHandler.getObject(0));
        auto p2 = dynamic_cast<MeltableParticle *>(particleHandler.getObject(1));
        auto i = dynamic_cast<MeltableInteraction *>(interactionHandler.getObject(0));
        logger.assert_debug(p1 or p2, "MeltableParticle required");
        logger.assert_debug(i, "MeltableInteraction required");
        //write to file
        if (getNumberOfTimeSteps()%saveCount==0)
            file << getTime()
                << ' ' << p1->getVelocity().Z
                << ' ' << i->getOverlap()
                << ' ' << i->getForce().Z
                << std::endl;
    }

    void printTime() const override
    {
        auto p = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(0));
        auto q = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(1));
        auto i = dynamic_cast<const MeltableInteraction*>(interactionHandler.getObject(0));
        logger(INFO, "time % temperatures % % overlap %", getTime(), p->getTemperature(), q->getTemperature(), i->getOverlap());
    }
};

int main()
{
    MeltableForceLawSelfTest(0.0);
    MeltableForceLawSelfTest(0.1);
    return 0;
}
