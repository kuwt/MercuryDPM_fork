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
#include <Boundaries/PeriodicBoundary.h>
//#include <Particles/ThermalParticle.h>

// define the species (type of particle interaction) we are going to use
typedef Species<ThermalSpecies<SinterNormalSpecies>, FrictionSpecies> ThermalSinterSpecies;

/**
 * Restarts the DepositAndSinter simulation and tests the tensile strength of the sample:
 * First, gravity is turned off, and the base is removed.
 * Then the left-most/right-most particles (the 'walls') get a prescribed velocity, pulling the particles apart.
 */
class TensileTest : public Mercury3D
{
    // boolean switching between horizonal and vertical tensile testing
    bool horizontalTest_;
    // Determines the speed at which the test is carried out
    double strainRate = 1.0;
    // How far the sample is going to be stretched, relative to its original size
    double maxStrain = 0.5;

public:
    /**
     * Restarts from DepositAndSinter, then changes the simulation to a tensile test
     *
     * \todo TW the failure criteria for normal/frictional forces should be based on beam theory
     */
    TensileTest(bool horizontalTest) : horizontalTest_(horizontalTest)
    {
        // Restart from DepositAndSinter
        setName("DepositAndSinter");
        readRestartFile();
        setRestarted(false);

        // Name the simulation TensileTest
        if (horizontalTest) {
            setName("TensileTestHorizontal");
        } else {
            setName("TensileTestVertical");
        }

        // Set simulation time based on maxStrain
        setTimeMax(maxStrain/strainRate);

        // Add adhesive and frictional forces to simulate the failure behaviour after sintering
        auto species = (ThermalSinterSpecies*) speciesHandler.getLastObject();
        // adhesive force such that the max adhesive force is proportional to the plastic overlap.
        species->setCohesionStiffness(species->getLoadingStiffness());
        // 'high' frictional forces such that the normal force is the first to fail
        species->setSlidingFrictionCoefficient(10.0);
        species->setRollingFrictionCoefficient(0.5); //behaves odd for values >1
        species->setTorsionFrictionCoefficient(0); // set to zero to keep it simple

        // Set gravity to zero
        setGravity(Vec3D(0,0,0));
        
        // Un-fix base particles
        for (BaseParticle* p : particleHandler) {
            if (p->isFixed()) p->unfix();
        }

        // Set wall velocity based on strain rate
        const double velocity = getWallVelocity();
        // Set width of wall region
        const double wallWidth = 2.0*particleHandler.getMeanRadius();

        // Move left/right-most particles with a fixed velocity (use bottom/top particles for the vertical tensile test)
        if (horizontalTest)
        {
            //turn off periodic boundary conditions
            boundaryHandler.removeObject(0);
            for (BaseParticle* p : particleHandler)
            {
                if (p->getPosition().X - getXMin() < wallWidth) {
                    p->fixParticle();
                    p->setPrescribedVelocity([velocity](double){return Vec3D(-velocity, 0, 0);});
                } else if (getXMax() - p->getPosition().X < wallWidth) {
                    p->fixParticle();
                    p->setPrescribedVelocity([velocity](double){return Vec3D(velocity, 0, 0);});
                }
            }
        } else {
            for (BaseParticle* p : particleHandler)
            {
                if (p->getPosition().Z - getZMin() < wallWidth) {
                    p->fixParticle();
                    p->setPrescribedVelocity([velocity](double){return Vec3D(0,0,-velocity);});
                } else if (getZMax() - p->getPosition().Z > wallWidth) {
                    p->fixParticle();
                    p->setPrescribedVelocity([velocity](double){return Vec3D(0,0,velocity);});
                }
            }
        }

        // Change the look of the xballs output (temperature is not needed anymore)
        setXBallsAdditionalArguments("-v0 -solidf -p 1");
    }

    /**
     * For simplicity, I set all parameters in the constructor
     */
    void setupInitialConditions() override
    {
    }

    /**
     * What to print to the console while the simulation is running
     */
    void printTime() const override
    {
        logger(INFO, "time % ene % disp % force %", getTime(),
               getKineticEnergy() / getElasticEnergy(), 2.0*getWallVelocity()*getTime(), getWallForce());
    }

    void writeEneHeader(std::ostream& os) const override
    {
        os << "time\tene\tdisp\tforce\n";
    }

    void writeEneTimestep(std::ostream& os) const override
    {
        os  << getTime() << '\t'
            << getKineticEnergy() / getElasticEnergy() << '\t'
            << 2.0*getWallVelocity()*getTime() << '\t'
            << getWallForce() << '\n';
    }


    //returns wall velocity based on strain rate
    double getWallVelocity() const {
        return 0.5*(getXMax()-getXMin())*strainRate;
    }

    //returns wall normal force (tensile=positive)
    double getWallForce() const {
        double force = 0;
        if (horizontalTest_)
        {
            double center = 0.5 * (getXMin() + getXMax());
            for (BaseParticle* p : particleHandler)
            {
                if (p->isFixed())
                {
                    if (p->getPosition().X < center)
                    {
                        force += p->getForce().X;
                    }
                    else
                    {
                        force -= p->getForce().X;
                    }
                }
            }
        } else {
            double center = 0.5 * (getZMin() + getZMax());
            for (BaseParticle* p : particleHandler)
            {
                if (p->isFixed())
                {
                    if (p->getPosition().Z < center)
                    {
                        force += p->getForce().Z;
                    }
                    else
                    {
                        force -= p->getForce().Z;
                    }
                }
            }
        }
        return force;
    }
};

/**
 * Main function; you can modify parameters here
 */
int main(int argc UNUSED, char* argv[] UNUSED)
{
    logger(INFO,"Horizontal tensile test");
    TensileTest horizontalTest(true);
    horizontalTest.solve();

    logger(INFO,"Vertical tensile test");
    TensileTest verticalTest(false);
    verticalTest.solve();
}
