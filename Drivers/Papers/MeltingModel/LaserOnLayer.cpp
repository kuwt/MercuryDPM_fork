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
#include "Initialization.h"
#include "Particles/MeltableParticle.h"

class LaserOnLayer: public Initialization
{
    double timeCool = 0.01;
    MeltableSpecies* species;

public:
    LaserOnLayer() {
        setScaleFactor(1e4);

        setName("InitialLayer");
        logger.assert_always(readRestartFile(),"Reading % failed", getName());
        setRestarted(false);
        setName("LaserOnLayer");
        removeOldFiles();

        species = dynamic_cast<MeltableSpecies *>(speciesHandler.getObject(0));
        species->setDensity(density);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(poissonRatio);
        species->setDissipation(dissipation);
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(solidHeatCapacity);
        species->setThermalConductivityCoefficient(thermalConductivity);
        species->setThermalConvectionCoefficient(thermalConvection);
        species->setAmbientTemperature(ambientTemperature);
        species->setMeltingTemperature(meltingTemperature);
        species->setMaterialEmissivity(emmisivity);
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(referenceViscosity);
        species->setDeltaT(temperatureInterval);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        //species->setWallTemperature(wallTemperature);
        //species->setThermalExpansionCoefficient();
        logger(INFO,"%", *species);

        using mathsFunc::square;
        const double laserPower = 5000; // 2.3 (W) -> Sintratec Kit laser power
        const double spotRadius = 250e-6; // 250
        const double porosity = 0.55;
        const double extinctionCoefficient = 3.0*(1-porosity)/(4.0*porosity*particleHandler.getMeanRadius());
        const double peakIntensity = 2.0*laserPower/constants::pi/square(spotRadius);
        const double surfaceHeight = getZMax();
        const double hatchWidth = 2.0*spotRadius;
        const double hatchLength = getXMax()-2.0*spotRadius;
        const double laserSpeed = 0.1;
        const double halfPeriod = (hatchWidth+hatchLength)/laserSpeed;
        const double corner = hatchLength/(hatchLength+hatchWidth)*halfPeriod; //when corner is reached
        const double laserTime = floor(hatchLength/hatchWidth)*halfPeriod+corner;
        logger(INFO,"peakIntensity % spotRadius % extinctionCoefficient % surfaceHeight % hatchWidth % hatchLength % laserSpeed % halfPeriod % corner % laserTime % meanRadius %", peakIntensity, spotRadius, extinctionCoefficient, surfaceHeight, hatchWidth, hatchLength, laserSpeed, halfPeriod, corner, laserTime, particleHandler.getMeanRadius());
        const double absorptivity = 0.33;

        setTimeMax(timeCool+laserTime);
        setSaveCount(getTimeMax()/getTimeStep()/500);

        // apply laser
        std::function<double(const BaseParticle*)> heatInput = [peakIntensity, spotRadius, extinctionCoefficient, surfaceHeight, hatchWidth, hatchLength, laserSpeed, halfPeriod, laserTime, corner, absorptivity] (const BaseParticle* p) {
            const double time = p->getHandler()->getDPMBase()->getTime();
            if (time>laserTime) return 0.0; // stop laser
            const int line = floor(time/halfPeriod) + 1;
            const double progress = std::fmod(time,halfPeriod);
            const double centerX = line % 2 == 0
                    ? spotRadius + (corner-std::min(progress,corner))*laserSpeed // even lines
                    : spotRadius + std::min(progress,corner)*laserSpeed; // odd lines
            const double centerY = spotRadius+(line-1)*hatchWidth + std::max(progress-corner,0.0) * laserSpeed;
            const double z = std::max(0.0,surfaceHeight-p->getPosition().Z);
            const double r2 = square(p->getPosition().X-centerX)+square(p->getPosition().Y-centerY);
            const double intensity = peakIntensity*exp(-extinctionCoefficient*z-2.0*r2/square(spotRadius));
            return absorptivity * constants::pi * square(p->getRadius()) * intensity;
        };
        species->setHeatInput(heatInput);

        for (const auto p : particleHandler) {
            auto q = dynamic_cast<MeltableParticle*>(p);
            q->setTemperature(meltingTemperature-0.5*temperatureInterval);
        }
        write(std::cout,false);

        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);
    };

    void printTime() const override {
        double maxTemperature = 0;
        double meanTemperature = 0;
        double volume = 0;
        for (auto p : particleHandler) {
            auto q = dynamic_cast<const MeltableParticle*>(p);
            //logger.assert_debug(q,"MeltableParticle needed");
            maxTemperature = std::max(maxTemperature,q->getTemperature());
            meanTemperature += q->getTemperature()*q->getVolume();
            volume += q->getVolume();
        }
        meanTemperature /= volume;
        logger(INFO, "time % timeMax % eneRatio % com % tempMean % tempMax %", getTime(), getTimeMax(), getKineticEnergy()/getGravitationalEnergy(),getCentreOfMass().Z,meanTemperature-meltingTemperature, maxTemperature-meltingTemperature);
    }
};

int main()
{
    LaserOnLayer problem;
    problem.solve();
    return 0;
}
