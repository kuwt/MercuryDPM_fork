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
#include "Particles/MeltableParticle.h"
#include "Initialization.h"
using constants::pi;
using std::pow;

class OneParticleHeatingSelfTest : public Initialization
{
private:
    double heatSource = 0.33*2*2.3/(pi*pow(250e-6,2))*pi*pow(50e-6,2);
    double timeMax = (temperatureInterval*(solidHeatCapacity+liquidHeatCapacity) + latentHeat)*mass/heatSource;
    double timeStep = 1e-3*timeMax;
    double initialTemperature = meltingTemperature-temperatureInterval;
    unsigned saveCount = 15;

public:
    //
    void setupInitialConditions() override {
        //latentHeat *= (1-pow(1-relMoltenLayerThicknessMax,3));
        setName("OneParticleHeatingSelfTest");
        setDomain({-radius, -radius, -radius},
                  {radius, radius, radius});
        setTimeMax(timeMax);
        setTimeStep(timeStep);
        setSaveCount(saveCount);

        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species->setDensity(density);
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(solidHeatCapacity); // solid heat capacity
        species->setMeltingTemperature(meltingTemperature);
        species->setHeatInput(heatSource);
        species->setDeltaT(temperatureInterval);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        logger(INFO,"%",*species);

        MeltableParticle p;
        p.setSpecies(species);
        p.setTemperature(initialTemperature);
        p.setRadius(radius);
        particleHandler.copyAndAddObject(p);
    }

    //write out data to plot
    void actionsAfterTimeStep() override {
        // get particle
        auto p = dynamic_cast<MeltableParticle *>(particleHandler.getLastObject());
        logger.assert_debug(p != nullptr, "MeltableParticle required");
        // compute energies
        double eneMelt = latentHeat*(mass-4./3.*constants::pi*density*mathsFunc::cubic(p->getSolidRadius()));
        double T = p->getTemperature();
        double Tl = meltingTemperature+0.5*temperatureInterval;
        double Ts = meltingTemperature-0.5*temperatureInterval;
        double eneHeat = solidHeatCapacity*mass*(T-initialTemperature);
        if (T>Ts) {
            eneHeat += 0.5*(liquidHeatCapacity-solidHeatCapacity)*mass*(T-Ts);
        }
        if (T>Tl) {
            eneHeat += 0.5*(liquidHeatCapacity-solidHeatCapacity)*mass*(T-Tl);
        }
        //write to file
        static std::ofstream file(getName()+".out");
        if (getNumberOfTimeSteps()%saveCount==0)
            file << getTime() << ' ' << p->getTemperature() << ' ' << p->getMoltenLayerThickness() << ' ' << eneMelt << ' ' << eneHeat << std::endl;
    }

    void printTime() const override
    {
        for (const auto& q: particleHandler) {
            auto p0 = dynamic_cast<MeltableParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            logger(INFO, "time % temperature % moltenLayer % radius % solidRadius % invMass % mass % density % volume %",
                    getTime(), p0->getTemperature(), p0->getMoltenLayerThickness(),
                   p0->getRadius(), p0->getSolidRadius(), p0->getInvMass(), p0->getMass(), p0->getSpecies()->getDensity(), p0->getVolume());
        }
    }
};

class OneParticleCoolingSelfTest : public Initialization
{
    MeltableSpecies* s = nullptr;

public:
    void setupInitialConditions() override {
        setName("OneParticleHeatingSelfTest");
        logger.assert_always(readRestartFile(),"Could not restart from %", getName());
        setTime(0);
        setTimeMax(1.1);
        setSaveCount(getTimeMax()/getTimeStep()/50);
        //
        setName("OneParticleCoolingSelfTest");

        s = dynamic_cast<MeltableSpecies *>(speciesHandler.getObject(0));
        logger.assert_debug(s != nullptr, "MeltableSpecies required");
        s->setHeatInput(0.0);
        //s->setThermalConductivityCoefficient(thermalConductivity);
        s->setThermalConvectionCoefficient(thermalConvection);
        s->setMaterialEmissivity(0.9);
        s->setAmbientTemperature(ambientTemperature);
        logger(INFO,"%",*s);
    }

    //write out data to plot
    void actionsAfterTimeStep() override {
        // get particle
        auto p = dynamic_cast<MeltableParticle *>(particleHandler.getLastObject());
        logger.assert_debug(p != nullptr, "MeltableParticle required");
        // compute energies
        static double initialTemperature = p->getTemperature();
        double mass = p->getMass();
        double eneMelt = s->getLatentHeat()*(mass-4./3.*constants::pi*s->getDensity()*mathsFunc::cubic(p->getRadius()-p->getMoltenLayerThickness()));
        double T = p->getTemperature();
        double Tl = s->getMeltingTemperature()+0.5*s->getDeltaT();
        double Ts = s->getMeltingTemperature()-0.5*s->getDeltaT();
        double eneHeat = s->getSolidHeatCapacity()*mass*(T-initialTemperature);
        if (T>Ts) {
            eneHeat += 0.5*(s->getLiquidHeatCapacity()-s->getSolidHeatCapacity())*mass*(T-Ts);
        }
        if (T>Tl) {
            eneHeat += 0.5*(s->getLiquidHeatCapacity()-s->getSolidHeatCapacity())*mass*(T-Tl);
        }
        //write to file
        static std::ofstream file(getName()+".out");
        if (getNumberOfTimeSteps()%dataFile.getSaveCount()==0)
            file << getTime() << ' ' << p->getTemperature() << ' ' << p->getMoltenLayerThickness() << ' ' << eneMelt << ' ' << eneHeat << std::endl;
    }

    void printTime() const override
    {
        for (const auto& q: particleHandler) {
            auto p0 = dynamic_cast<MeltableParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            logger(INFO, "time % temperature % moltenLayer % radius % solidRadius % invMass % mass % density % volume %",
                   getTime(), p0->getTemperature(), p0->getMoltenLayerThickness(),
                   p0->getRadius(), p0->getSolidRadius(), p0->getInvMass(), p0->getMass(), p0->getSpecies()->getDensity(), p0->getVolume());
        }
    }
};

int main()
{
    {
        OneParticleHeatingSelfTest problem;
        problem.solve();
    }
    {
        OneParticleCoolingSelfTest problem;
        problem.solve();
    }
    return 0;
}
