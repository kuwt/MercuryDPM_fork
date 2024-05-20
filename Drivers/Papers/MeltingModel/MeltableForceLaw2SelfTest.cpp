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

class MeltableForceLaw2SelfTest : public Initialization
{
private:
    bool writeVTK = false;
    Mdouble temperature1 = meltingTemperature+0.5*temperatureInterval;
    Mdouble temperature2 = meltingTemperature+0.5*temperatureInterval;
    double timeMax = 0.1; //0.5 to see deltamax
    unsigned saveCount = timeMax/timeStep/1000;
    std::ofstream file;

public:
    MeltableForceLaw2SelfTest (double scaleFactor) {
        setScaleFactor(scaleFactor);
        saveCount = std::max(1.0,saveCount/sqrt(scaleFactor));
        solve();
    }

    void setupInitialConditions() override {
        setName("MeltableForceLaw2SelfTest_" + helpers::toString(getScaleFactor()));
        removeOldFiles();
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
        //species->setHeatInput(heatSource);
        //species->setThermalConvectionCoefficient(thermalConvection);
        //species->setMaterialEmissivity(emmisivity);
        species->setAmbientTemperature(ambientTemperature);
        species->setMeltingTemperature(meltingTemperature);
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(referenceViscosity);
        species->setDeltaT(temperatureInterval);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        species->setActivationEnergy(activationEnergy);
        //species->setThermalExpansionCoefficient(thermalExpansionCoeff);
        logger(INFO,"%",*species);
        species->analyseTimeScales(radius, species->getDensity(), temperature1);

        setTimeStep(timeStep);
        setSaveCount(saveCount);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);

        double overlap = 1e-10*radius;

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

//        auto i = dynamic_cast<MeltableInteraction*>(interactionHandler.getInteraction(
//                particleHandler.getObject(0),particleHandler.getObject(1),0));
//        i->setBondingOverlap(.55*overlap);
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
                << ' ' << p1->getMoltenLayerThickness()
                << ' ' << p2->getMoltenLayerThickness()
                << ' ' << i->getOverlap()
                << ' ' << i->getForce().Z
                << std::endl;
    }

    void printTime() const override
    {
        auto p = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(0));
        auto q = dynamic_cast<const MeltableParticle *>(particleHandler.getObject(1));
        auto i = dynamic_cast<const MeltableInteraction*>(interactionHandler.getObject(0));
        logger(INFO, "time % temperatures % % overlap % neck %", getTime(), p->getTemperature(), q->getTemperature(), i->getOverlap()/radius, i->getNeckRadius());
    }
};

int main(int argc, char* argv[])
{
    MeltableForceLaw2SelfTest(1e8);
    // run with argument "-full" to get output needed for paper
    if (helpers::readFromCommandLine(argc,argv,"-full")) {
        MeltableForceLaw2SelfTest(1e6);
        MeltableForceLaw2SelfTest(1e4);
        MeltableForceLaw2SelfTest(1e0);
    }
    return 0;
}
