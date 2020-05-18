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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Chute.h"

class ChutePeriodic : public Chute
{
public:
    
    void actionsBeforeTimeStep()
    {};

///creates flow particles in the whole chute
    BaseParticle create_inflow_particle()
    {
        BaseParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        P0.setPosition(Vec3D(random.getRandomNumber(getXMin() + 2.0 * P0.getRadius(), getXMax()),
                             random.getRandomNumber(getYMin() + 2.0 * P0.getRadius(), getYMax()),
                             random.getRandomNumber(getZMin() + 2.0 * P0.getRadius(), getInflowHeight())));
        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        return P0;
    }
    
    void setupInitialConditions()
    {
        
        Chute::setupInitialConditions();
        
        PeriodicBoundary b;
        b.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b);
        b.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b);
        
        numberOfParticles = 10000;
        
        add_flow_particles();
    }
    
    ///Add initial flow particles.
    void add_flow_particles()
    {
        logger(INFO, "Adding flowing particles");
        writeRestartFile();
        unsigned int N = particleHandler.getNumberOfObjects() + getChuteLength() * getChuteWidth() * getInflowHeight() /
                                                                mathsFunc::cubic(getInflowParticleRadius()) / 8;
        particleHandler.setStorageCapacity((N));
        Mdouble H = getInflowHeight();
        setZMax(1.2 * getInflowHeight());
        writeRestartFile();
        //try to find new insertable particles
        while (particleHandler.getSize() < N)
        {
            BaseParticle p = create_inflow_particle();
            if (checkParticleForInteraction(p))
            {
                particleHandler.copyAndAddObject(p);
            }
            else
            {
                setInflowHeight(getInflowHeight() + .0001 * getMaxInflowParticleRadius());
            }
        }
        logger(INFO, "InflowHeight= %", getInflowHeight());
        setHGridDistribution(HGridDistribution::EXPONENTIAL);
    }
    
    void setup()
    {
        // Problem parameters
        setName("chute_periodic");
        fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        autoNumber();
        setTimeMax(10);
        
        // Particle properties
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(2400.0);
        setInflowParticleRadius(.5e-3);
        species.setCollisionTimeAndRestitutionCoefficient(4e-4, 0.8, 1);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingStiffness(species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        setFixedParticleRadius(getInflowParticleRadius());
        setRoughBottomType(MULTILAYER);
        
        // Chute properties
        setChuteAngle(25.0);
        setChuteLength(20e-3);
        setChuteWidth(10e-3);
        setZMax(20e-3);
        setMaxFailed(1000);
        makeChutePeriodic();
        
        setTimeStep(4.0e-4 / 50);
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(150, getTimeMax(), getTimeStep()));
        logger(INFO, "dt = %", getTimeStep());
        
        //This set to colouring based of size and small vectors
        setXBallsColourMode(7);
        setXBallsVectorScale(1);
    }
    
    unsigned int numberOfParticles;
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    ChutePeriodic problem;
    problem.setup();
    logger.assert_always(argc > 1, "Argument needed!");
    problem.setChuteAngle(atof(argv[1]));
    logger(INFO, "Chute Angle: %", problem.getChuteAngle());
    problem.setZMax(5e-3);
    problem.solve();
    problem.write(std::cout, false);
    problem.writeRestartFile();
    
    ChutePeriodic problem2;
    problem2.setup();
    logger.assert_always(argc > 1, "Argument needed!");
   problem2.setChuteAngle(atof(argv[1]));
   
    logger(INFO, "Chute Angle: %", problem2.getChuteAngle());
    problem2.setZMax(10e-3);
    problem2.solve();
    problem2.write(std::cout, false);
    problem2.writeRestartFile();
    
    ChutePeriodic problem3;
    problem3.setup();
    logger.assert_always(argc > 1, "Argument needed!");
     problem3.setChuteAngle(atof(argv[1]));
    logger(INFO, "Chute Angle: %", problem3.getChuteAngle());
    problem3.setZMax(15e-3);
    problem3.solve();
    problem3.write(std::cout, false);
    problem3.writeRestartFile();
    
    ChutePeriodic problem4;
    problem4.setup();
    logger.assert_always(argc > 1, "Argument needed!");
    problem4.setChuteAngle(atof(argv[1]));
    logger(INFO, "Chute Angle: %" , problem4.getChuteAngle());
    problem4.setZMax(20e-3);
    problem4.solve();
    problem4.write(std::cout, false);
    problem4.writeRestartFile();
    
    ChutePeriodic problem5;
    problem5.setup();
    logger.assert_always(argc > 1, "Argument needed!");
    problem5.setChuteAngle(atof(argv[1]));
    logger(INFO, "Chute Angle: %", problem5.getChuteAngle());
    problem5.setZMax(25e-3);
    problem5.solve();
    problem5.write(std::cout, false);
    problem5.writeRestartFile();
}
