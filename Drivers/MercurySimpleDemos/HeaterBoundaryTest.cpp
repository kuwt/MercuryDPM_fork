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

#include<iostream>
#include <Species/LinearViscoelasticSpecies.h>

#include "Mercury3D.h"
#include "Boundaries/HeaterBoundary.h"
#include "Walls/InfiniteWall.h"
#include "MercuryTime.h"

/*!
 * A test of the HeaterBoundary.
 */

class HeaterBoundaryTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        const unsigned int N1 = static_cast<unsigned int>(pow(N, 0.33)) + 1;
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        for (unsigned int i = 0; i < N; ++i)
        {
            const unsigned int ix = (i % N1);
            //\todo IFCD: I added the static cast to double, since it seemed to me that integer division is not what's meant here. If you think the double is not needed, please remove both casts.
            const unsigned int iz = static_cast<unsigned int>( static_cast<double>(i) / N1 / N1);
            const unsigned int iy = (i - ix - N1 * N1 * iz) / N1;

            const double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            const double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            const double z = (getZMax() - getZMin()) * (iz + 1) / (N1 + 1);

            p0.setPosition(Vec3D(x, y, z));
//            p0.setVelocity(20 * Vec3D(0.1, 0.1, 0.1));
            p0.setVelocity(Vec3D(0,0,0));
            p0.setRadius(0.002);
            particleHandler.copyAndAddObject(p0);
        }

        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);

        HeaterBoundary b0;
        // b0.set3D(Vec3D(0.045,0.045,0.045), Vec3D(0.075,0.075,0.075), 1/particleHandler.getMass() );
        b0.set3D(Vec3D(0.0,0.0,0.0), Vec3D(0.1,0.1,0.1), 1/particleHandler.getMass() );
        boundaryHandler.copyAndAddObject(b0);
    }

    ///\brief Number of particles in the system.
    unsigned int N;
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    Time time;
    time.tic();

    HeaterBoundaryTest problem;
    LinearViscoelasticSpecies species;
    species.setDensity(2000);
    species.setDissipation(0.005);
    species.setStiffness(1e3);
    problem.speciesHandler.copyAndAddObject(species);

    problem.setXMax(0.1);
    problem.setYMax(0.1);
    problem.setZMax(0.1);
    problem.N = 1000;
    problem.setName("HeaterBoundaryTest");
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setTimeStep(4e-5);
//    logger(INFO, "dt / tc = % / % = 1/%", 
//            problem.getTimeStep(), species.getCollisionTime( 2000 * 4./3. * constants::pi * pow(0.002,3)),
//            1/(problem.getTimeStep() / species.getCollisionTime( 2000 * 4./3. * constants::pi * pow(0.002,3))));
    problem.setSaveCount(25);
    problem.setTimeMax(0.3);
    problem.setSystemDimensions(3);
    
    problem.setHGridMaxLevels(1);
    problem.setHGridCellOverSizeRatio(1.2);
    problem.setHGridUpdateEachTimeStep(false);

    problem.dataFile.setFileType(FileType::MULTIPLE_FILES);

    problem.solve();
    
    logger(INFO, "Total time to run this simulation: % s", time.toc());
}
