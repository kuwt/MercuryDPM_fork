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
#include "Species/SinterFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/SphericalWall.h"
//#include "Walls/InfiniteWall.h"
//#include "Walls/AxisymmetricIntersectionOfWalls.h"

/**
 * Single particle, indented slowly by spherical indenter.
 */
class Indenter : public Mercury3D
{
public:

    Indenter(int argc, char *argv[], Mdouble indenterRadius, Mdouble indentationVelocity, Mdouble indentationForce)
     : indentationVelocity_(indentationVelocity), indentationForce_(indentationForce)
    {
        logger(INFO, "Creating indenter with radius %, force %, velocity %", indenterRadius, indentationForce,
               indentationVelocity);
        logger.assert(indenterRadius > 0.0, "indenterRadius has to be positive");
        logger.assert(indentationVelocity > 0.0, "indentationVelocity has to be positive");
        logger.assert(indentationForce > 0.0, "indentationForce has to be positive");


        std::stringstream addToName;
        for (unsigned i=1; i<argc; i++)
        {
            if (!strcmp(argv[i], "-H"))
            {
                addToName << 'H' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-W"))
            {
                addToName << 'W' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-E"))
            {
                addToName << 'E' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-D"))
            {
                addToName << 'D' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-R"))
            {
                addToName << 'R' << atof(argv[i + 1]);
            }
        }
        setName("SinterBed1"+addToName.str());
        readRestartFile();
        setRestarted(false);
        setName("SinterBed2"+addToName.str());
        logger(INFO, "Name: %", getName());

        //Set species properties
        particleSpecies = dynamic_cast<SinterFrictionSpecies *>(speciesHandler.getObject(0));
        logger.assert(particleSpecies, "Particle species does not exist");
        indenterSpecies = dynamic_cast<SinterFrictionSpecies *>(speciesHandler.getObject(1));
        logger.assert(indenterSpecies, "Indenter species does not exist");
        indenterParticleSpecies = speciesHandler.getMixedObject(indenterSpecies, particleSpecies);

        particleSpecies->setSlidingFrictionCoefficient(100);
        particleSpecies->setRollingFrictionCoefficient(10);
        particleSpecies->setTorsionFrictionCoefficient(1);
        particleSpecies->setSinterRate(0);
        indenterParticleSpecies->setCohesionStiffness(0);

        Mdouble mass = particleSpecies->getMassFromRadius(particleHandler.getMeanRadius());
        logger(INFO, "old time step: %", getTimeStep());
        setTimeStep(particleSpecies->computeTimeStep(mass));
        logger(INFO, "new time step: %", getTimeStep());

        //Set wall-indenter properties
        indenter_ = wallHandler.copyAndAddObject(SphericalWall(indenterRadius, indenterSpecies));
        indenter_->setVelocity({0, 0, -indentationVelocity});
        initialBedHeight_ = getBedHeight();
        indenter_->setPosition({0, 0, initialBedHeight_ + indenterRadius});
    }

    void outputXBallsData(std::ostream &os) const override
    {
        //Adds one line to the particle data
        os << particleHandler.getNumberOfObjects() + 1 << " " << getTime() << " "
        << getXMin() << " " << getYMin() << " " << getZMin() << " "
        << getXMax() << " " << getYMax() << " " << getZMax() << " " << std::endl;
        //This outputs the particle data
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            outputXBallsDataParticle(i, 14, os);
        //Adds indenter
        os << indenter_->getPosition() << " "
        << indenter_->getVelocity() << " "
        << indenter_->getRadius() << "  0 0 0 0 0 0 0\n";
    }

    ///solve until system is relaxed
    bool isRelaxed() const
    {
        // Counter from 1 to 100, then reset.
        // Used to check only every 100 time steps
        static unsigned counter = 0;
        counter = (counter > 100) ? 1 : counter++;
        // Counter that is increased if the system is relaxed, reset to zero if not.
        // Simulation is discontinued if system is relaxed for 300 time steps
        static unsigned isStatic = 0;
        if (counter == 0)
        {
            if (getKineticEnergy() < 1e-5 * getElasticEnergy())
            {
                if (isStatic > 3)
                {
                    //logger(INFO, "t=%: System is relaxed", getTime());
                    return true;
                }
                else
                {
                    isStatic++;
                }
            }
            else
            {
                isStatic = 0;
            }
        }
        return false;
    }

    bool continueSolve() const override
    {
        static unsigned stage = 0;

        if (stage == 0)
        { //moving down
            if (indenter_->getForce().Z > indentationForce_)
            {
                logger(INFO, "t=%: Stage 0 (load) complete", getTime());
                stage++;
                indenter_->setVelocity({0, 0, 0});
            }
        }
        else if (stage == 1)
        { //stopping until system is relaxed
            if (isRelaxed())
            {
                logger(INFO, "t=%: Stage 1 (relax) complete", getTime());
                stage++;
                indenter_->setVelocity({0, 0, indentationVelocity_});
            }
        }
        else
        { //retracting until force on indenter is zero for a certain time interval
            if (indenter_->getForce().Z == 0)
            {
                if (isRelaxed())
                {
                    logger(INFO, "t=%: Stage 2 (unload) complete", getTime());
                    return false;
                }
            }
        }
        return true;
    }

    /** top of the particles (can be used as initial height of the indenter)) */
    Mdouble getBedHeight()
    {
        Mdouble bedHeight = getZMin();
        for (auto p : particleHandler)
        {
            if (!p->isFixed())
            {
                bedHeight = std::max(bedHeight, p->getPosition().Z + p->getRadius());
            }
        }
        logger(INFO, "Bed height: %, particle number: %", bedHeight, particleHandler.getNumberOfObjects());
        return bedHeight;
    }

    void printTime() const override
    {
        static int counter = 0;
        if (++counter == 1)
            writeEneHeader(std::cout);
        writeEneTimestep(std::cout);
    }

    /** creates custom ene header */
    void writeEneHeader(std::ostream &os) const override
    {
    }

    /** creates custom ene output */
    void writeEneTimestep(std::ostream &os) const override
    {
        Mdouble diameter = 2.0 * particleHandler.getObject(0)->getRadius();
        os << "time " << getTime()
        << " EneKin/EneEla " << std::left << std::setw(12) << getKineticEnergy() / getElasticEnergy()
        << " displacement/d " << (initialBedHeight_ + indenter_->getRadius() - indenter_->getPosition().Z) / diameter
        << " force/forceMax " << indenter_->getForce().Z / indentationForce_ << std::endl;
    }

private:
    SinterFrictionSpecies *particleSpecies; //pointer to the particle properties
    SinterFrictionSpecies *indenterSpecies; //pointer to the particle properties
    SinterFrictionMixedSpecies *indenterParticleSpecies; //pointer to the wall-particle collision properties
    SphericalWall *indenter_;
    Mdouble indentationVelocity_;
    Mdouble indentationForce_;
    Mdouble initialBedHeight_;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble timeMax = 0.01;
    Mdouble indentationDepth = 20e-6; //1um (corresponds to contact radius of sqrt(127*1)=11um)
    Mdouble indenterRadius = 0.5 * 127e-6;
    Mdouble indentationVelocity = indentationDepth / timeMax * 2.0;
    Mdouble indentationForce = 4e-3 * 1e-7; //4mN

    Indenter i(argc, argv, indenterRadius, indentationVelocity, indentationForce);
    i.setTimeStep(0.8 * i.getTimeStep());
    i.setTimeMax(timeMax);
    i.solve();

    helpers::writeToFile("SinterBed2.gnu",
                         "set xlabel 'displacement/d [um]'\n"
                          "set ylabel 'forceRel [mN]'\n"
                          "p 'SinterBed2.ene' u 6:8 w lp\n"
    );
    std::cout << "Execute 'gnuplot SinterBed2.gnu' to view output" << std::endl;
    std::cout << "Execute './SinterBed2.xballs -drwmax 4.01e-5 -cmode 2' to view xballs output" << std::endl;
}
