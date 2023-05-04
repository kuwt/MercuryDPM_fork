//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Species/Species.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include <string>
#include <sstream>
#include <cstring>

using mathsFunc::cubic;
using mathsFunc::square;
using constants::pi;

/**
 * Creates a 3D silo with steady state conditions (particles removed at teh base get re-inserted at the top)
 */
class Silo : public Mercury3D {
public:

    Silo(int argc, char *argv[]) {
        // Check if there are at least two arguments, and the first argument is "-r";
        // in that case, we restart an old code, using the second argument as fileName.
        if (argc>2 && !strcmp(argv[1],"-r"))
        {
            // Read the filename and restart from the command line
            std::string fileName = argv[2];
            // Restart
            logger(INFO, "Restarting %", fileName);
            setName(fileName);
            readRestartFile();
            // currently, the restarter overwrites the original data, fstat and restart files with new data;
            // uncomment this line to append to the original data and fstat files instead
            //Silo.setAppend(true);
        }
    }

    void setupInitialConditions() override {

        // make cross-checks of the input variables, ensuring that they are consistent with each other.
        logger.assert_always(outflowHeight_<=siloHeight_,"outflowHeight cannot be larger than siloHeight");
        logger.assert_always(outflowRadius_<=siloRadius_,"outflowRadius cannot be larger than siloRadius");

        Mdouble particleBulkVolume = particleNumber_*cubic(2.0*particleRadius_);
        Mdouble siloVolume = pi*(siloHeight_-outflowHeight_)*siloRadius_*siloRadius_ +
         pi/3.0*outflowHeight_*(square(siloRadius_)+siloRadius_*outflowRadius_+square(outflowRadius_));

        logger.assert_always(particleBulkVolume<=siloVolume,"particleBulkVolume % cannot be larger than siloVolume %", particleBulkVolume, siloVolume);
        logger(INFO,"Silo will be about %% full % %", particleBulkVolume/siloVolume*100, '%', particleBulkVolume, siloVolume);

        //setting system size
        setMin({-siloRadius_, -siloRadius_, 0});
        setMax({siloRadius_, siloRadius_, siloHeight_});

        setGravity(Vec3D(0, 0, -9.8));//setting gravity

        // set species (if not already set)
        if (speciesHandler.getNumberOfObjects()==0) {
            Mdouble maximumOverlap = 0.01 * (2.0*particleRadius_);
            Mdouble restitution = 0.1;
            //define the particle properties
            auto s = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
            s->setDensity(1e3);
            Mdouble mass = s->getMassFromRadius(particleRadius_);
            Mdouble maximumForce = std::fabs(siloHeight_/(2.0*particleRadius_)*mass*getGravity().Z);
            s->setStiffnessAndRestitutionCoefficient(maximumForce/maximumOverlap, restitution, mass);
            logger(INFO,"Set collision time to % such that overlaps are below %*particleRadius", s->getCollisionTime(mass), maximumOverlap/particleRadius_);
            s->setSlidingFrictionCoefficient(0.5);
            s->setSlidingStiffness(2.0/7.0*s->getStiffness());
            s->setSlidingDissipation(2.0/7.0*s->getDissipation());
            speciesHandler.copyAndAddObject(s);
            setTimeStep(0.02*s->getCollisionTime(mass));
        }


        //add walls

        //add outer cylindrical wall
        AxisymmetricIntersectionOfWalls w0({0,0,0}, {0,0,1}, {}, speciesHandler.getObject(0));
        w0.addObject({1,0,0},{siloRadius_,0,0});
        wallHandler.copyAndAddObject(w0);

        //add outer funnel wall
        AxisymmetricIntersectionOfWalls w1({0,0,0}, {0,0,1}, {}, speciesHandler.getObject(0));
        w1.addObject({outflowHeight_,0,-(siloRadius_-outflowRadius_)},{outflowRadius_,0,0});
        w1.addObject({1e-3,0,1},{outflowRadius_,0,0});
        wallHandler.copyAndAddObject(w1);
        
        // uncomment to add flat wall at base, thus preventing outflow
        // InfiniteWall w({0,0,-1}, {0,0,0}, speciesHandler.getObject(0));
        // wallHandler.copyAndAddObject(w);

        //add particles
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        //creating a random number whose magnitude is defined by 'polydispersity'
        //to give particles within the system a randomised distribution of sizes
        //within the desired bounds
        p.setRadius(particleRadius_);
        Mdouble zMax = 0;
        Mdouble dz = 0.0005*particleRadius_; //this number can be adjusted to make particle insertion quicker
        while (particleHandler.getNumberOfObjects() < particleNumber_) {
            // define the position at which the particle is to be inserted
            Mdouble x, y, rr, z = random.getRandomNumber(0, zMax);
            do {
                x = random.getRandomNumber(-siloRadius_, siloRadius_);
                y = random.getRandomNumber(-siloRadius_, siloRadius_);
                rr = x * x + y * y;
            } while (rr > mathsFunc::square(siloRadius_ - particleRadius_) || rr < mathsFunc::square(0.5 * siloRadius_ + particleRadius_));
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                double pol = random.getRandomNumber(-polydispersity_,polydispersity_);
                p.setRadius(particleRadius_*(1+pol));
                //std::cout << '+';
                //std::cout.flush();
            } else {
                //gradually increase the insertion domain to to insert particles as low as possible
                zMax += dz;
                //std::cerr << "warning: failed insertion attempt" << std::endl;
            }
        }
        logger(INFO,"All particles inserted below z=%*siloHeight", zMax/siloHeight_);
    }

    void printTime () const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", EneKin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                  << std::endl;
        std::cout.flush();
    }

    ///creating a reinsertion boundary
    void actionsAfterTimeStep() override {
        static unsigned nParticlesReinserted = 0;
        //take all particles below siloHeight_ z=-5r and reinsert them into the Silo above siloHeight_ zMax
        for (auto p: particleHandler) {
            if (p->getPosition().Z < -5.0*particleRadius_) {
                p->setVelocity({0, 0, 0});
                Mdouble zMax = 0;
                do {
                    zMax += 1e-3;
                    Mdouble x, y, rr, z = random.getRandomNumber(getZMax(), getZMax() + zMax);
                    do {
                        x = random.getRandomNumber(getXMin(), getXMax());
                        y = random.getRandomNumber(getYMin(), getYMax());
                        rr = x * x + y * y;
                    } while (rr > mathsFunc::square(getXMax() - p->getRadius()) ||
                             rr < mathsFunc::square(0.5 * getXMax() + p->getRadius()));
                    p->setPosition({x, y, z});
                } while (!checkParticleForInteraction(*p));
                nParticlesReinserted++;
                if (nParticlesReinserted%100==0) {
                    logger(INFO,"reinserted %% of all particles",(double)nParticlesReinserted/particleHandler.getNumberOfObjects()*100,'%');
                }
            }
        }
    }

    //scaling up the particle number
    void scaleUp(double factor)
    {
        particleRadius_ /= std::cbrt(factor);
        particleNumber_ *= factor;
    }

public: //setters and getters

    Mdouble getParticleRadius_() const
    {
        return particleRadius_;
    }

    void setParticleRadius_(Mdouble particleRadius_)
    {
        logger.assert_always(particleRadius_>0,"particleRadius should be positive");
        logger(INFO,"Setting particleRadius = %",particleRadius_);
        Silo::particleRadius_ = particleRadius_;
    }

    unsigned int getParticleNumber_() const
    {
        return particleNumber_;
    }

    void setParticleNumber_(unsigned int particleNumber_)
    {
        logger.assert_always(particleNumber_>0,"particleNumber should be 1 or larger");
        logger(INFO,"Setting particleNumber = %",particleNumber_);
        Silo::particleNumber_ = particleNumber_;
    }

    Mdouble getSiloRadius_() const
    {
        return siloRadius_;
    }

    void setSiloRadius_(Mdouble siloRadius_)
    {
        logger.assert_always(siloRadius_>0,"siloRadius should be positive");
        logger(INFO,"Setting siloRadius = %",siloRadius_);
        Silo::siloRadius_ = siloRadius_;
    }

    Mdouble getSiloHeight_() const
    {
        return siloHeight_;
    }

    void setSiloHeight_(Mdouble siloHeight_)
    {
        logger.assert_always(siloHeight_>0,"siloHeight should be positive");
        logger(INFO,"Setting siloHeight = %",siloHeight_);
        Silo::siloHeight_ = siloHeight_;
    }

    Mdouble getOutflowRadius_() const
    {
        return outflowRadius_;
    }

    void setOutflowRadius_(Mdouble outflowRadius_)
    {
        logger.assert_always(outflowRadius_>0,"outflowRadius should be positive");
        logger(INFO,"Setting outflowRadius = %",outflowRadius_);
        Silo::outflowRadius_ = outflowRadius_;
    }

    Mdouble getOutflowHeight_() const
    {
        return outflowHeight_;
    }

    void setOutflowHeight_(Mdouble outflowHeight_)
    {
        logger.assert_always(outflowHeight_>0,"outflowHeight should be positive");
        logger(INFO,"Setting outflowHeight = %",outflowHeight_);
        Silo::outflowHeight_ = outflowHeight_;
    }

    double getPolydispersity_() const
    {
        return polydispersity_;
    }

    void setPolydispersity_(double polydispersity_)
    {
        logger.assert_always(polydispersity_>=0,"polydispersity should be non-negative");
        logger(INFO,"Setting polydispersity = %",polydispersity_);
        Silo::polydispersity_ = polydispersity_;
    }

private:

    /// The radius possessed by particles
    Mdouble particleRadius_ = 5e-2;
    /// The number of particles within the system
    unsigned particleNumber_=1000;
    ///the (maximum) radius of the silo
    Mdouble siloRadius_=0.5;
    ///The height of the silo
    Mdouble siloHeight_=1.5;
    ///the radius of the outflow funnel
    Mdouble outflowRadius_=0.25;
    ///the height of the outflow funnel
    Mdouble outflowHeight_=0.0;
    ///The (fractional) degree of polydispersity (rMax/rMin-1) within the system
    double polydispersity_=0;
};

int main(int argc, char *argv[])
{
    //setting up the problem
    Silo Silo(argc,argv);
    Silo.setTimeMax(10.0); //the total duration of the simulation
    //Silo.setSaveCount(200); //setting how often to output data to file
    Silo.setSaveCount(200); //setting how often to output data to file
    //Silo.dataFile.setSaveCount(2000); //setting how often to output data to file
    Silo.fStatFile.setSaveCount(2000); //setting how often to output data to file
    Silo.restartFile.setSaveCount(10000); //setting how often to output data to file
    Silo.fStatFile.setFileType(FileType::NO_FILE);
    //Silo.dataFile.setFileType(FileType::NO_FILE);
    //Silo.scaleUp(8);
    Silo.setName("Silo3D"); //naming the output file
    Silo.setXBallsAdditionalArguments("-v0 -solidf");
    Silo.solve();
    return 0;
}
