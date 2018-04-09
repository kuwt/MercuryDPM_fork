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
#include "Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"

class Bed : public Mercury3D
{
public:

    //set slave variables (i.e. compute stiffness, dissipation, timestep, wall and particle positions)
    Bed(int argc, char *argv[])
    {
        //create new simulation; as no specialised functions are needed,
        //Mercury3D is used directly instead of creating a derived class as usual.
        Mdouble radius = 2.0; //start with the largest particle size
        Mdouble domainHeight = 5.0; // in particle layers
        Mdouble domainWidth = 0.25; // in indenter diameter
        unsigned bedDimension = 3;
        //real: 1e7 times lower than real
        Mdouble elasticModulus = 3e2; //!!!! using overly soft particles
        //set other fixed parameters
        Mdouble restitutionCoefficient = 0.5; //overly dissipative particles
        Mdouble density = 1040; //1040 kg/m^3
        Mdouble sinterAdhesion = 0.0; //sintering not active
        Mdouble inverseSinterViscosity = 0.0; //sintering not active
        Mdouble gravity = 9.8;

        std::stringstream addToName;
        for (unsigned i=1; i<argc; i++)
        {
            if (!strcmp(argv[i], "-H"))
            {
                domainHeight = atof(argv[i + 1]); //in particle layers
                addToName << 'H' << domainHeight;
            }
            else if (!strcmp(argv[i], "-W"))
            {
                domainWidth = atof(argv[i + 1]);  //as fraction of the indenter diameter
                addToName << 'W' << domainWidth;
            }
            else if (!strcmp(argv[i], "-E"))
            {
                elasticModulus = atof(argv[i + 1]);  //in Pa
                addToName << 'E' << elasticModulus;
            }
            else if (!strcmp(argv[i], "-D"))
            {
                bedDimension = atoi(argv[i + 1]); //2D or 3D
                addToName << 'D' << bedDimension;
            }
            else if (!strcmp(argv[i], "-R"))
            {
                radius = atof(argv[i + 1]); //in nanometer
                addToName << 'R' << radius;
            }
        }
        setName("SinterBed0"+addToName.str());
        logger(INFO,"Name: %",getName());

        //scale into si units
        radius *= 1e-6;
        domainHeight *= 2.0*radius;
        domainWidth *= 127e-6;

        //set time-stepping and output parameters
        setFileType(FileType::ONE_FILE);
        setXBallsAdditionalArguments(" -v0 -solidf ");

        //set time-stepping and output parameters
        setSaveCount(200);
        setTimeMax(1e20);

        //set gravity
        setGravity(Vec3D(0, 0, -gravity));

        //add Species
        auto s = speciesHandler.copyAndAddObject(SinterFrictionSpecies());
        //density
        s->setDensity(density);
        //stiffness and dissipation
        Mdouble mass = s->getMassFromRadius(radius);
        //10% overlap k=4/3*3e9*sqrt(0.1)*4e-6=5000
        Mdouble stiffness = 4./3.*elasticModulus*sqrt(0.1)*radius;
        s->setStiffnessAndRestitutionCoefficient(stiffness, restitutionCoefficient, mass);
        //overly quick simulation, based on loading stiffness only!
        setTimeStep(5.0*s->computeTimeStep(mass));
        //plastic parameters
        s->setUnloadingStiffnessMax(3.0 * s->getLoadingStiffness());
        s->setPenetrationDepthMax(1.0);
        //friction
        s->setSlidingDissipation(2./7. * s->getDissipation());
        s->setSlidingStiffness(2./7. * s->getLoadingStiffness());
        s->setRollingDissipation(2./5. * s->getDissipation());
        s->setRollingStiffness(2./5. * s->getLoadingStiffness());
        s->setTorsionDissipation(2./5. * s->getDissipation());
        s->setTorsionStiffness(2./5. * s->getLoadingStiffness());
        s->setSinterAdhesion(sinterAdhesion);
        s->setInverseSinterViscosity(inverseSinterViscosity);
        ///\todo TW energy seems to behave erratically in mu_ro>=mu_sl; do they have to be in order?
        s->setSlidingFrictionCoefficient(0.5);
        s->setRollingFrictionCoefficient(0.2);
        s->setTorsionFrictionCoefficient(0.1);
        //output species statistics
        logger(INFO, "Added species");
        logger(INFO, "  fg=g*m=%", mass * getGravity().Z);
        logger(INFO, "  fad=ad*r=%", sinterAdhesion * radius);
        logger(INFO, "  ts=r^4*nu/ad=%",
               radius * radius * radius * radius / inverseSinterViscosity / sinterAdhesion);
        logger(INFO, "  tg=sqrt(d/g)=%", sqrt(2.0 * radius / getGravity().getLength()));
        logger(INFO, "  tc=sqrt(m/k)=%", sqrt(mass / s->getLoadingStiffness()));
        auto w = speciesHandler.copyAndAddObject(s); //wall species

        //set domain size
        setZMax(domainHeight); //typical thickness 20d
        setXMax(0.5*domainWidth); //2um indentation corresponds to contact radius of sqrt(127*2)=16um)
        setYMax(getXMax());
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(0);

        // the walls are set based on xMin, xMax, ..., zMax
        std::cout << "Inserting walls" << std::endl;
        InfiniteWall baseWall({0, 0, -1}, {0, 0, 0}, w);
        wallHandler.copyAndAddObject(baseWall);

//        // define side wall
//        AxisymmetricIntersectionOfWalls sideWall;
//        sideWall.setSpecies(w);
//        sideWall.setPosition({0, 0, 0});
//        sideWall.setOrientation({0, 0, 1});
//        sideWall.addObject({1, 0, 0}, {domainRadius, 0, 0});
//        wallHandler.copyAndAddObject(sideWall);
        auto b = PeriodicBoundary();
        b.set({1,0,0},getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(b);
        b.set({0,1,0},getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(b);

        // determine how many particles need to be inserted
        Mdouble domainVolume;
        Mdouble particleUnitVolume;
        if (bedDimension == 2)
        {
            domainVolume = domainWidth * domainHeight;
            particleUnitVolume = mathsFunc::square(radius) * 4.0;
        }
        else
        {
            domainVolume = domainWidth * domainWidth * domainHeight;
            particleUnitVolume = mathsFunc::cubic(radius) * 8.0;
        }
        Mdouble numberOfParticles = domainVolume / particleUnitVolume;
        logger(INFO, "Inserting % particles", std::floor(numberOfParticles));

        // Inserting particles; 50% of the particles are large (r_s=0.95*d), 50% small (r_s=1.05*d)
        ///small polydispersity
        Mdouble large = 1.05 * radius;
        Mdouble small = 0.95 * radius;
        Mdouble numberOfParticlesInserted = 0;
        Mdouble maxHeight = 0.0;
        BaseParticle p;
        p.setRadius(large);
        p.setSpecies(s);
        while (numberOfParticlesInserted < numberOfParticles)
        {
            //set random position
            Vec3D pos;
            pos.X = random.getRandomNumber(getXMin() + p.getRadius(), getXMax() - p.getRadius());
            if (bedDimension != 2)
                pos.Y = random.getRandomNumber(getYMin() + p.getRadius(), getYMax() - p.getRadius());
            pos.Z = random.getRandomNumber(getZMin() + p.getRadius(), maxHeight);
            p.setPosition(pos);

            //Mdouble r=sqrt(mathsFunc::square(pos.X)+mathsFunc::square(pos.Y));
            //if (r<(getXMax()-getXMin())/2.0-p.getRadius() && checkParticleForInteraction(p))
            if (checkParticleForInteraction(p))
            {
                particleHandler.copyAndAddObject(p);
                p.setRadius(random.getRandomNumber(small, large));
                ++numberOfParticlesInserted;
                std::cout << ".";
            }
            else
            {
                //slowly increase height so particles can get inserted
                maxHeight += .0002 * radius;
            }
        }
        logger(INFO, "Insertion complete");
    }

    void writeEneHeader(std::ostream &os) const override
    {
    }

    //override printTime function such that console output shows the state of relaxation (Ekin/Eela)
    void printTime() const override
    {
        writeEneTimestep(std::cout);
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
                    logger(INFO, "t=%: System is relaxed", getTime());
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

    ///solve until system is relaxed
    bool continueSolve() const override
    {
        return !isRelaxed();
    }

    void writeEneTimestep(std::ostream &os) const override
    {
        Mdouble overlap = 0.0;
        Mdouble plasticOverlap = 0.0;
        for (auto i:interactionHandler)
        {
            overlap += i->getOverlap();
            auto j = dynamic_cast<const SinterInteraction *>(i);
            if (j)
            {
                plasticOverlap += j->getPlasticOverlap();
            }
            else
            {
                auto j = dynamic_cast<const LinearPlasticViscoelasticInteraction *>(i);
                logger.assert(true, "Species is neither Sinter nor LinearPlasticViscoelastic");
                plasticOverlap += j->getMaxOverlap();
            }
        }
        unsigned N = interactionHandler.getNumberOfObjects();
        overlap /= N;
        plasticOverlap /= N;

        os
        << "time " << std::left << std::setw(12) << getTime()
        << " eneRatio " << std::left << std::setw(12) << getKineticEnergy() / getElasticEnergy()
        << " delta " << std::left << std::setw(12) << overlap
        << " delta0 " << std::left << std::setw(12) << plasticOverlap
        << std::endl;
    }

};

int main(int argc, char *argv[])
{
    //Description: Sintering.tex in Dropbox/Sintering
    // real time scales:
    // - collision time: sqrt(m/k)  = 50 us
    // - gravitational time sqrt(d/g) = 1 ms
    // - sinter time: 80 s

    //initialise bed
    Bed bed(argc,argv);

    //run
    bed.solve();

    helpers::writeToFile("SinterBed0.gnu",
                         "set logscale y\n"
                          "set xlabel 'time [s]'\n"
                          "set ylabel 'kin/ela Energy []'\n"
                          "p 'SinterBed0.ene' u 2:4 w lp\n"
    );
    logger(INFO, "Execute 'gnuplot SinterBed0.gnu' to view output");

}
