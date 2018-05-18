// This file is part of MercuryDPM.
// 
// MercuryDPM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// MercuryDPM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with MercuryDPM.  If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 The Mercury Developers Team
// For the list of developers, see <http://www.MercuryDPM.org/Team>

#include "Mercury3D.h"
#include <Species/SinterFrictionSpecies.h>
#include <Species/NormalForceSpecies/ThermalSpecies.h>
#include <map>
#include <random>
#include <Boundaries/PeriodicBoundary.h>
//#include <Particles/ThermalParticle.h>

// define the species (type of particle interaction) we are going to use
typedef Species<ThermalSpecies<SinterNormalSpecies>, FrictionSpecies> ThermalSinterFrictionSpecies;

/**
 * Creates a cube of length l, periodic in x- and y-directions, with a lower wall in z-direction, constrained by gravity.
 * The system is filled repeatedly with enough particles to form a layer of height h.
 * These particles have a fixed temperature (backgroundTemperature), radius (mean: r, stdev: dr)
 * After filling, we wait until the system is relaxed (eneKin<0.02*enePot).
 * Then the top layer (all particles above h-2r, h being the height of the highest particle)
 * receives a certain amount of thermal energy (sinterEnergy).
 * The contact law is chosen such that all particles above the glass temperature start sinter (glassTemperature),
 * while the particles below that temperature are stable.
 * The thermal conductivity determines how quickly temperature is dissipated through the contacts (setThermalConductivity).
 * A new layer is added when the temperature drops below a given limit (maxTemperatureAtEndOfSinter).
 * This process is repeated until the cube is filled.
 *
 * SI-units are used.
 */
class DepositAndSinter : public Mercury3D
{
    // Parameters with default values, can be changed by the user in the main() function

    /// mean particle radius
    double r = 50e-6;
    /// standard deviation of particle radii
    double dr = 10e-6;
    /// cube length
    double l = 20 * r;
    /// layer height
    double h = 4 * r;
    /// system dimensions (1: chain of particles in z-direction, 2: x-z layer, 3: x-y-z layer)
    unsigned dim = 2;
    /// stores the species properties
    ThermalSinterFrictionSpecies* species;
    /// initial and boundary value for temperature
    double backgroundTemperature = 400;
    /// temperature above which sintering starts
    double glassTemperature = 420;
    /// sintering ends if the temperature drops below this temperature
    double maxTemperatureAtEndOfSinter = backgroundTemperature + 0.5*(glassTemperature-backgroundTemperature);
    /// how much power (in Watt/l^2) to add to the top layer
    /// (note, the number of particles in the top layer, and thus the temperature added, strongly depends on dim)
    double sinterEnergy = .2e-3;

public:

    /**
     * Setting default values (these can be changed by the user in main())
     */
    DepositAndSinter()
    {
        setName("DepositAndSinter");
        setGravity(Vec3D(0, 0, -9.8));
        //species parameters that do not depend on the particle radius
        species = speciesHandler.copyAndAddObject(ThermalSinterFrictionSpecies());
        species->setHeatCapacity(500);
        species->setThermalConductivity(20);
        species->setDensity(1000);
        //species->setSlidingFrictionCoefficient(0.1);
        species->setPenetrationDepthMax(0.5); // stops sintering if the overlap is larger than penetrationDepth*radius/2
        //setFileType(FileType::ONE_FILE);
        setXBallsAdditionalArguments(" -v0 -solidf -p 1 -cmode 8"); // modification to xballs output such that the particles are colored by temperature
        setSaveCount(200); //  output is written every saveCount-th time step
        setTimeMax(1e20); //do not set a time limit; instead, stop when the cube is filled
        setTimeStep(1); // this parameter is set in setupInitialConditions
        //setParticlesWriteVTK(true); //comment to turn off VTK output (needed to view the simulation in paraview)
    }

    /**
     * Setting all other simulation parameters that depend on the user-provided parameters
     */
    void setupInitialConditions() override
    {
        //if no species has been set yet, create a default species
        if (species->getLoadingStiffness() == 0)
        {
            // softness of particles the smaller, the stiffer
            double collisionTime = 5e-4;
            // determines dissipation of the particle contacts, how much relative velocity remains after one collision time
            double restitutionCoefficient = 0.001;
            // mass of smallest particle
            double mass = species->getMassFromRadius(r-2*dr);
            species->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, mass);
            species->setSlidingStiffness(2. / 7. * species->getLoadingStiffness());
            species->setSlidingDissipation(2. / 7. * species->getDissipation());
            species->setRollingStiffness(2. / 5. * species->getLoadingStiffness());
            species->setRollingDissipation(2. / 5. * species->getDissipation());
            species->setTorsionStiffness(2. / 5. * species->getLoadingStiffness());
            species->setTorsionDissipation(2. / 5. * species->getDissipation());
            // sintering only changes teh plastic overlap, not the actual overlap;
            // the  prefactor here determines the effect of of the plastic overlap on the actual overlap
            species->setUnloadingStiffnessMax(4.0 * species->getLoadingStiffness());
            // sintering law
            species->setSinterType(SINTERTYPE::TEMPERATURE_DEPENDENT_FRENKEL);
            // the rate of change of the plastic overlap is equal to fn/fa*sr, with
            // fn the normal force (due to gravity),
            // fa the sinterAdhesion (due to surface tension, which is dominant over van der Walls)
            // sr the sinterRate (a temperature-dependent proportionality constant; here we assume that sr=0 below the glass temperature)
            species->setSinterAdhesion(12e-9/r);
            species->setTemperatureDependentSinterRate([this] (double temperature) { return 3.0*(temperature>glassTemperature); });
        }

        //set time step based on the collision time
        double mass = species->getMassFromRadius(r);
        setTimeStep(0.02 * species->getCollisionTime(mass));

        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(l, l, l));

        logger(INFO, "Setting up periodic boundary conditions in x and y");
        if (dim > 1)
        {
            PeriodicBoundary p;
            p.set({1, 0, 0}, getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(p);
            if (dim == 3)
            {
                p.set({0, 1, 0}, getYMin(), getYMax());
                boundaryHandler.copyAndAddObject(p);
            }
        }

        insertBaseLayer();
        insertParticleLayer();
    }

    /**
     * This function is called after each time step. It allows modifications during the simulation.
     */
    void actionsAfterTimeStep() override
    {
        //stores if we are in the relaxation or in the sintering phase
        static bool relaxation = true;

        //the next two lines are to execute the following code every 50th time step
        static unsigned counter = 0;
        if (++counter < 50) return; else counter = 0;

        if (relaxation) {
            // relaxation mode:
            // once energy has died down, stop relaxation and start sintering
            if (getKineticEnergy() < 0.02 * getElasticEnergy())
            {
                relaxation = false;
                addSinterEnergy();
                logger(INFO, "Start sintering (t=%)",getTime());
            }
        } else {
            // sintering mode:
            // once temperature is below limit, insert new particle layer and start relaxation
            if (getMaximumTemperature()<maxTemperatureAtEndOfSinter) {
                relaxation = true;
                if (getCenterOfMass() < l / 2)
                {
                    insertParticleLayer();
                    logger(INFO, "Start relaxation (t=%)",getTime());
                } else {
                    setTimeMax(getTime());
                    logger(INFO, "Finish simulation (t=%)",getTime());
                }
            }
        }
    }

    // Inserts a layer of fixed particles on which the moving particles get deposited
    void insertBaseLayer() {
        ThermalParticle p;
        p.setSpecies(species);
        p.setRadius(r);
        p.setTemperature(backgroundTemperature);
        p.fixParticle();
        if (dim==1) {
            p.setPosition(Vec3D(l/2,l/2,-p.getRadius()));
            particleHandler.copyAndAddObject(p);
            //logger(ERROR,"p %",p);
        } else if (dim==2) {
            double nx = ceil(l/(2*r));
            for (double i=0.5; i<nx;++i)
            {
                p.setPosition(Vec3D(i/nx*l, l / 2, -p.getRadius()));
                particleHandler.copyAndAddObject(p);
            }
        } else {
            double nx = ceil(l/(2*r));
            double ny = ceil(l/(2*r));
            for (double i=0.5; i<nx;++i)
                for (double j=0.5; i<nx;++i)
                {
                    p.setPosition(Vec3D(i/nx*l,j/nx*l, -p.getRadius()));
                    particleHandler.copyAndAddObject(p);
                }
        }

    }

    // Call this function to insert a new particle layer
    void insertParticleLayer()
    {
        //set up random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        //set up normal distribution
        std::normal_distribution<> d(r, dr);

        //the volume to be added (assumes a volume fraction of 0.65)
        Mdouble addVolume;
        ThermalParticle p;
        p.setSpecies(species);
        p.setRadius(r);
        p.setTemperature(backgroundTemperature); //also divide by thermal density

        if (dim == 1)
        {
            addVolume = h / (2 * r) * p.getVolume() - 0.5 * p.getVolume();
        }
        else if (dim == 2)
        {
            addVolume = h / (2 * r) * l / (2 * r) * p.getVolume() - 0.5 * p.getVolume();
        }
        else
        {
            addVolume = 0.65 * h * l * l - 0.5 * p.getVolume();
        }

        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        Mdouble fillHeight = getZMin();
        while (addVolume > 0)
        {
            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
            if (dim == 1) x = l / 2;
            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
            if (dim < 3) y = l / 2;
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
            p.setPosition(Vec3D(x, y, z));
            // check if particle can be inserted
            if (checkParticleForInteraction(p))
            {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do
                {
                    p.setRadius(d(gen));
                } while (fabs(p.getRadius() - r) > 2 * dr); //reject too small or large radii
            }
            else
            {
                fillHeight += 0.01 * r; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO, " Inserted new layer, t=% N=%", getTime(), particleHandler.getNumberOfObjects());
    }

    // Call this function to insert thermal energy into the top particle layer
    void addSinterEnergy()
    {
        //add temperature above a certain height
        double h = getSurfaceHeight() - 2*r;
        //count particles in top layer
        double massHeatedParticles = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->getPosition().Z >= h)
            {
                massHeatedParticles += p->getMass();
            }
        }
        if (massHeatedParticles == 0)
        {
            logger(ERROR, "no particles to heat could be found; increasing height of top layer", massHeatedParticles, h);
        }
        //heat particles in top layer
        for (BaseParticle* p : particleHandler)
        {
            if (p->getPosition().Z > h)
            {
                ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
                tp->addTemperature(sinterEnergy / massHeatedParticles / species->getHeatCapacity());
            }
        }
    }

    /**
     * What to print to the console while the simulation is running
     */
    void printTime() const override
    {
        logger(INFO, "time=%, eneKin/enePot=%, eneTherm=%, maxTemp=%, h/l=%", getTime(),
               getKineticEnergy() / getElasticEnergy(), getThermalEnergy(),getMaximumTemperature(),getSurfaceHeight()/l);
    }

    /**
     * Prints the temperature into the data file, such that you can plot it in paraview
     */
    double getInfo(const BaseParticle& p) const override
    {
        auto tp = dynamic_cast<const ThermalParticle*>(&p);
        return tp->getTemperature();
    }

    /**
     * Returns the thermal energy present in the particles
     */
    Mdouble getThermalEnergy() const
    {
        Mdouble eneThermal = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneThermal += tp->getTemperature() * tp->getMass() * species->getHeatCapacity();
        }
        return eneThermal;
    }

    /**
     * Returns the maximum temperature in the particles
     */
    Mdouble getMaximumTemperature() const
    {
        Mdouble temperatureMax = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            if(tp->getTemperature()>temperatureMax) temperatureMax = tp->getTemperature();
        }
        return temperatureMax;
    }

    /**
     * Returns the center of mass of all particles
     */
    Mdouble getCenterOfMass() const
    {
        double momentum = 0.0;
        double mass = 0.0;
        for (const BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;
            mass += p->getMass();
            momentum += p->getMass() * p->getPosition().Z;
        }
        //logger(INFO,"mass %",mass);
        return momentum / mass;
    }

    /**
     * Returns the height of the highest particle
     */
    Mdouble getSurfaceHeight() const
    {
        double height = 0;
        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height;
    }

    /**
     * Returns the mean overlap of all particle contacts
     */
    Mdouble getMeanOverlap() const
    {
        double overlap = 0.0;
        for (const BaseInteraction* p : interactionHandler)
        {
            overlap += p->getOverlap();
        }
        return overlap / interactionHandler.getNumberOfObjects();
    }

};

/**
 * Main function; you can modify parameters here
 */
int main(int argc UNUSED, char* argv[] UNUSED)
{
    DepositAndSinter problem;
    problem.solve();
}
