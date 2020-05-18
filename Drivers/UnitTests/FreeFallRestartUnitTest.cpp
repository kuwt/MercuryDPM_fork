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


///This code tests:
///1) Restarting
///2) Saving arcoss multiple files
///3) and accepting command line argument.
///4) Also tests restart reloading.

#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <Logger.h>

class FreeFall : public Mercury2D
{
public:
    
    void setupInitialConditions() override {
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
        wallHandler.copyAndAddObject(w);

        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);
        p.setPosition(Vec3D(p.getRadius(), 2.0*p.getRadius(), 0));
        particleHandler.copyAndAddObject(p);

        //particle reaches ground after r=gt^2/2 => t=sqrt(d/g)=1
        setGravity(Vec3D(0,-1,0));
    }

};


void runFreeFall(int argc, char* argv[])
{
    ///Start off my solving the default problem
    FreeFall dpm;
    dpm.setName("FreeFallRestart");
    dpm.setMax({1.0,1.5,0.0});
    auto species = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(4.0/constants::pi);
    species->setCollisionTimeAndRestitutionCoefficient(0.1,1,species->getMassFromRadius(0.5));
    dpm.setTimeStep(0.002);
    dpm.setTimeMax(2.1);
    dpm.setSaveCount(100);
    dpm.solve(argc, argv);
}

int main(int argc, char* argv[])
{
    //If code is being called with no arguments (e.g. as a selftest), it enters here and call the code again 3 times with different arguments.
    if (argc == 1)
    {
        logger(INFO, "Case 1: simple run");
        if (system("./FreeFallRestartUnitTest -name FreeFallRestart0"))
            logger(FATAL, "code did not run");

        logger(INFO, "Case 2: restart in the middle");
        //restarted at t=0.2

        if (system("./FreeFallRestartUnitTest -tmax 1.05 -name FreeFallRestart1"))
            logger(FATAL, "code did not run");
        
        if (system("./FreeFallRestartUnitTest -r FreeFallRestart1 -tmax 2.1"))
            logger(FATAL, "code did not run");
        
        
        logger(INFO, "Case 3: restart in the middle, using separate data files");

        if (system("./FreeFallRestartUnitTest -tmax 1.05 -name FreeFallRestart2 -fileTypeData 2"))
            logger(FATAL, "code did not run");

        if (system("./FreeFallRestartUnitTest -r FreeFallRestart2 -tmax 2.1"))
            logger(FATAL, "code did not run");
    }
    else
    {
        runFreeFall(argc, argv);
        return 0;
    }
    
    //To compare the code using gnuplot follow the instructions below.
    //gnuplot> plot 'free_fall_restart/free_fall_restart.ene' u 1:3 w l, 'free_fall_restart/free_fall_no_restart.ene' u 1:3 w l
    //../sc/fpdiff.py ./free_fall_restart.fstat ./free_fall_no_restart.fstat
    
    //Final stage now we check what we get.
    logger(INFO, "Finished running, now comparing");
    
    FreeFall dpm0;
    FreeFall dpm1;
    FreeFall dpm2;
    
    dpm0.readRestartFile("FreeFallRestart0.restart");
    dpm1.readRestartFile("FreeFallRestart1.restart");
    dpm2.readRestartFile("FreeFallRestart2.restart");
    
    auto p1 = dpm1.particleHandler.begin();
    auto p2 = dpm2.particleHandler.begin();

    for (auto p0 = dpm0.particleHandler.begin(); p0 != dpm0.particleHandler.end(); ++p0)
    {
        if (!(*p0)->getPosition().isEqualTo((*p1)->getPosition(),1e-6))
        {
            logger(FATAL, "Particles is not in the same place after restart. Before it was % and now it is %.",
                   (*p0)->getPosition(), (*p1)->getPosition());
        }
        if (!(*p0)->getPosition().isEqualTo((*p2)->getPosition(),1e-10))
        {
            logger(FATAL, "Particles velocities are not the same place. Before it was % and now it is %.",
                   (*p0)->getVelocity(), (*p1)->getVelocity());
        }
        ++p1;
        ++p2;
    }
    
    return 0;
}
