//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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
#include "Particles/HeatFluidCoupledParticle.h"
#include "Species/HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Walls/InfiniteWall.h"

/// Tests whether the evaporation model works
/// Based on verification test in Sahar's qualifier
class EvaporationAndHeatTest : public Mercury3D {
public:

    void setupInitialConditions () override {
        // Set name of output files
        setName("EvaporationAndHeatUnitTest");
        // set gravity
        setGravity(Vec3D(0, 0, -9.8));
        // Radius of the particle
        double radius = 3e-3;
        //Contact law and density
        auto species = speciesHandler.copyAndAddObject(HeatFluidCoupledLinearViscoelasticFrictionLiquidMigrationWilletSpecies());
        species->setDensity(1270); // Density of the particle (kg/m^3)
        double mass = species->getMassFromRadius(radius);
        species->setStiffness(mass * 9.8 / (0.01 * radius));
        species->setRestitutionCoefficient(0.4, mass);
        species->setSlidingFrictionCoefficient(0.5);
        species->setSlidingStiffness(2. / 7. * species->getStiffness());
        species->setSlidingDissipation(2. / 7. * species->getDissipation());
        species->setRollingFrictionCoefficient(0.5);
        species->setRollingStiffness(2. / 5. * species->getStiffness());
        species->setRollingDissipation(2. / 5. * species->getDissipation());
        //! I am removing liquid bridges for now because they are not evaporating yet
        species->setLiquidBridgeVolumeMax(1e-9);
        species->setSurfaceTension(72.8e-3);
        //species->setHeatCapacity(1e100); // Particle heat capacity (J/KgK) //exaggerated
        species->setHeatCapacity(1470); // Particle heat capacity (J/KgK)
        species->setThermalConductivity(0.5); // Particle heat conductivity
        species->setMassTransferCoefficient(1); //The mass transfer coefficient (m/s)
        species->setLatentHeatVaporization(2.256e6); // Latent heat of vaporisation of free water(J/kg)
        species->setLiquidDensity(1000); // density of liquid (kg/m^3)
        species->setEvaporationCoefficientA(0);   //The coefficients a and b for the grain
        species->setEvaporationCoefficientB(0); //are determined as 3.2 and -21.7 [29]
        species->setAmbientHumidity(0.1); // Relative humadity of drying air
        species->setAmbientEquilibriumMoistureContent(0.0); // Equilibrium moisture content under the condition of the drying air
        species->setAmbientVapourConcentration(0.0); // Vapour concentrations in the drying medium (kg/m^3)
        species->setAmbientTemperature(343); // Drying air temperature (K)
        //species->setContactAngle();
        double collisionTime = species->getCollisionTime(mass);
        // Set time step, final time and how often to output
        setTimeStep(0.05 * collisionTime);
        setTimeMax(8.0);
        setSaveCount(static_cast<int>(getTimeMax() / getTimeStep() / 100));
        // Set domain size
        setMin({-radius, -radius, 0});
        setMax({radius, radius, 4.0 * radius});

        // Add two particles
        HeatFluidCoupledParticle particle;
        particle.setSpecies(species);
        particle.setRadius(radius);
        particle.setPosition(Vec3D(0, 0, radius));
        particle.setLiquidVolume(6e-9); // Initial moisture volume (m^3)
        particle.setTemperature(298); // Initial temperature of the particle (K)
        particleHandler.copyAndAddObject(particle);
        //particle.setTemperature(298); // Initial temperature of the particle (K)
        particle.setTemperature(278); // Initial temperature of the particle (K)
        particle.setPosition(Vec3D(0, 0, 3.0 * radius));
        particleHandler.copyAndAddObject(particle);

        // add base wall
        InfiniteWall base;
        base.setSpecies(species);
        base.set(Vec3D(0, 0, -1), Vec3D(0, 0, 0));
        wallHandler.copyAndAddObject(base);
    }

    /**
     * outputs:
     *  - time
     *  - liquid film volumes of the two particles,
     *  - the total liquid bridge volume
     *  - temperatures of the two particles
     *  to the console
     */
    void printTime() const override {
        auto p0 = dynamic_cast<const HeatFluidCoupledParticle*>(particleHandler.getObject(0));
        auto p1 = dynamic_cast<const HeatFluidCoupledParticle*>(particleHandler.getLastObject());
        logger(INFO, "t % VF % % VB % T % %", getTime(),
               p0->getLiquidVolume(),
               p1->getLiquidVolume(),
               interactionHandler.getLiquidBridgeVolume(),
               p0->getTemperature(),
               p1->getTemperature()
        );
    }
};

/// Solve problem
int main() {
    EvaporationAndHeatTest problem;
    //problem.setWallsWriteVTK(true);
    //problem.setParticlesWriteVTK(true);
    //problem.setInteractionsWriteVTK(true);
    problem.solve();
    // check final liquid film and temperature values
    auto p0 = dynamic_cast<const HeatFluidCoupledParticle*>(problem.particleHandler.getObject(0));
    helpers::check(p0->getLiquidVolume(),0,1e-8,"LiquidFilmVolume");
    helpers::check(p0->getTemperature(),228.191,1e-3,"Temperature");
    return 0;
}
