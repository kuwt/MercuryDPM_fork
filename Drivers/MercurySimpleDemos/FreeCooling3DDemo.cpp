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

//! [FCD_3D:headers]
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Mercury3D.h"
#include "MercuryTime.h"
//! [FCD_3D:headers]

/*
 * In this file 10^3 particles with the same velocity are placed in a tri-axial box. 
 * This makes them collide with the walls and eachother. 
 * Afterwards the same run is performed with hgrid on. 
 * It tests the working (and speedup) of the hgrid.
 */

//! [FCD_3D:class]
class FreeCooling3DDemoProblem : public Mercury3D
{
public:

    void setupInitialConditions() override {
        const unsigned int N1 = static_cast<unsigned int>(pow(N, 0.33)) + 1;
        //! [FCD_3D:particle]
        particleHandler.clear();
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        //! [FCD_3D:particle]
        //! [FCD_3D:placeparticles]
        for (unsigned int i = 0; i < N; ++i)
        {
            const unsigned int ix = (i % N1);
            const unsigned int iz = static_cast<unsigned int>( static_cast<double>(i) / N1 / N1);
            const unsigned int iy = (i - ix - N1 * N1 * iz) / N1;
            // set particle position
            const double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            const double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            const double z = (getZMax() - getZMin()) * (iz + 1) / (N1 + 1);
            p0.setPosition(Vec3D(x, y, z));
            // set random velocities for the particle
            p0.setVelocity(Vec3D(random.getRandomNumber(-2.0,2.0),  random.getRandomNumber(-2.0,2.0), random.getRandomNumber(-2.0,2.0)));
            p0.setRadius(0.0025);
            particleHandler.copyAndAddObject(p0);
        }
        //! [FCD_3D:placeparticles]

        //! [FCD_3D:CO_Mass]
        // Compute the center of mass velocity
        double particle_mass =  p0.getMass();
        double M_b = N*particle_mass; // mass of the bulk system
        Vec3D V_com = {0,0,0};

        for (int k = 0; k < particleHandler.getNumberOfObjects() ; k++){
            BaseParticle* p = particleHandler.getObject(k);
            V_com +=  (particle_mass*p->getVelocity())/M_b;
        }

        // Compute the reduced velocity for each particle
        for (int k = 0; k < particleHandler.getNumberOfObjects() ; k++){
            BaseParticle* p = particleHandler.getObject(k);
            p->setVelocity(p->getVelocity() - V_com);
        }
        //! [FCD_3D:CO_Mass]
        //! [FCD_3D:walls]
        PeriodicBoundary pb;
        pb.set(Vec3D(1,0,0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb);
        pb.set(Vec3D(0,1,0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(pb);
        pb.set(Vec3D(0,0,1), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(pb);
        //! [FCD_3D:walls]
    }

    //! [FCD_3D:aftertime]
    void actionsAfterTimeStep() override{
        if (getTime() > 4e5*getTimeStep()) {
            FC3D_Species.setDissipation(0.232);
        }
    }
    //! [FCD_3D:aftertime]

    ///\brief Number of particles in the system.
    //! [FCD_3D:datamembers]
    unsigned int N;
    LinearViscoelasticSpecies FC3D_Species;
    //! [FCD_3D:datamembers]
};
//! [FCD_3D:class]

//! [FCD_3D:main]
int main(int argc UNUSED, char* argv[] UNUSED)
{
    // Problem setup
    FreeCooling3DDemoProblem problem;
    //! [FCD_3D:species]
    LinearViscoelasticSpecies species;
    species.setDensity(2e3);
    species.setDissipation(0.0);
    species.setStiffness(1e3);
    problem.FC3D_Species = species;
    problem.speciesHandler.copyAndAddObject(species);
    //! [FCD_3D:species]
    //! [FCD_3DproblemSetup]
    problem.N = 1000;
    problem.setName("FreeCooling3DDemo");
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setTimeStep(5e-5);
    problem.setSaveCount(4000);
    problem.setTimeMax(50.0);
    problem.setMax(0.064,0.064,0.064);
    problem.setHGridMaxLevels(1);
    problem.setHGridCellOverSizeRatio(1.2);
    problem.setHGridUpdateEachTimeStep(false);
    //! [FCD_3DproblemSetup]
    //! [FCD_3D:solve]
    problem.setFileType(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    problem.solve();
    //! [FCD_3D:solve]
}
//! [FCD_3D:main]