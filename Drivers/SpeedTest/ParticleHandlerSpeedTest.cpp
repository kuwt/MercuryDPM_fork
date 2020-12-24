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

#include <Boundaries/CubeInsertionBoundary.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <MercuryTime.h>


int main()
{
    // create a particle handler
    Mercury3D dpm;
    auto species = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    SphericalParticle particle(species);
    CubeInsertionBoundary insertionBoundary;
    insertionBoundary.set(particle,1e6,{0,0,0},{10000,10000,10000},{0,0,0},{0,0,0},.5,.5);
    insertionBoundary.setInitialVolume(1e5*constants::pi/6.);
    insertionBoundary.checkBoundaryBeforeTimeStep(&dpm);
    size_t n = dpm.particleHandler.getSize();
    logger(INFO,"Number of particles: %",n);
    logger(INFO,"Size of particles: % bytes, % doubles",sizeof(particle),sizeof(particle)/sizeof(double));

    // create a data vector
    std::vector<std::array<double,82>> dataVector(n);
    for (int i = 0; i < n; ++i) {
        dataVector[i][0] = dpm.particleHandler.getObject(i)->getPosition().X;
    }

    // create a small data vector
    std::vector<std::array<double,8>> smallDataVector(n);
    for (int i = 0; i < n; ++i) {
        smallDataVector[i][0] = dpm.particleHandler.getObject(i)->getPosition().X;
    }

    // create a particle vector
    std::vector<SphericalParticle> particleVector(n);
    for (int i = 0; i < n; ++i) {
        particleVector[i].setPosition(dpm.particleHandler.getObject(i)->getPosition());
    }

    Time timer;
    double sum=0;
    int repetitions = 1e8;

    for (int i = 0; i < repetitions; ++i) {
        sum += dpm.particleHandler.getObject(rand() % n)->getPosition().X;
    }
    double refTime = timer.toctic();
    logger(INFO, "Time to access randomly data in the particle handler: % s (% pct)", refTime, 100);

    for (int i = 0; i < repetitions; ++i) {
        sum += dpm.particleHandler.getObject(i % n)->getPosition().X;
    }
    double time = timer.toctic();
    int relTime = 100. * time / refTime;
    logger(INFO, "Time to access ordered data in a particle handler: % s (% pct)", time, relTime);

    for (int i = 0; i < repetitions; ++i) {
        sum += particleVector[rand() % n].getPosition().X;
    }
    time = timer.toctic();
    relTime = 100. * time / refTime;
    logger(INFO, "Time to access randomly data in a particle vector: % s (% pct)", time, relTime);

    for (int i = 0; i < repetitions; ++i) {
        sum += dataVector[rand() % n][0];
    }
    time = timer.toctic();
    relTime = 100. * time / refTime;
    logger(INFO, "Time to access randomly data in a vector of arrays: % s (% pct)", time, relTime);

    for (int i = 0; i < repetitions; ++i) {
        sum += smallDataVector[rand() % n][0];
    }
    time = timer.toctic();
    relTime = 100. * time / refTime;
    logger(INFO, "Time to access randomly data in a small vector of small arrays: % s (% pct)", time, relTime);

    logger(INFO,"Computation: %",sum);

    logger(INFO,"Conclusions:\n"
                " - 50\% gain: Ordered data access is much quicker than random\n"
                " - 20\% gain: Vectors are quicker than handlers\n"
                " -  0\% gain: Class type (SphericalParticle vs array<double,83>) makes no difference\n"
                " - 20\% gain: A lean array (8 doubles) is read quicker than a fat array (83 doubles) ");

}

