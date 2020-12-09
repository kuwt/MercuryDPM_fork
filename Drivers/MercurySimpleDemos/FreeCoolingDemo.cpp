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

//! [FCD_2D:headers]
#include<iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Mercury2D.h"
//! [FCD_2D:headers]

/*
 * In this file 32^2 particles with the same velocity are placed in a bi-axial box. This makes them collide with
 * the walls and eachother. Afterwards the same run is performed with hgrid on. It tests the working
 * (and speedup) of the hgrid.
 */

//! [FCD_2D:class]
class FreeCoolingDemoProblem : public Mercury2D{
public:
    void setupInitialConditions() override {
        //! [FCD_2D:species]
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        FC2D_Species = species;
        species->setDensity(10000);
        species->setDissipation(0.0);
        species->setStiffness(1e3);
        //! [FCD_2D:species]
        //! [FCD_2D:particle]
        particleHandler.clear();
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        //! [FCD_2D:particle]
        //! [FCD_2D:placeparticles]
        int N1=static_cast<int>(sqrt(N))+1;
        for (int i=0;i<N;i++) {
            int ix = i % N1;
            int iy = i / N1;
            // set particle position
            double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            p0.setPosition(Vec3D(x, y, 0.0));
            // set random velocities for the particle
            p0.setVelocity(Vec3D(random.getRandomNumber(-0.001,0.001),  random.getRandomNumber(-0.001,0.001), 0.0));
            p0.setRadius(0.0002);
            p0.setSpecies(species);
            particleHandler.copyAndAddObject(p0);
        }
        //! [FCD_2D:placeparticles]

        //! [FCD_2D:CO_Mass]
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
        //! [FCD_2D:CO_Mass]
        //! [FCD_2D:walls]
        PeriodicBoundary pb;
        pb.set(Vec3D(1,0,0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb);
        pb.set(Vec3D(0,1,0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(pb);
        //! [FCD_2D:walls]
    }
    //! [FCD_2D:aftertime]
    void actionsAfterTimeStep() override{
        // Time to switch on dissipation
        if (getTime() > 4e5*getTimeStep())
        {
            FC2D_Species->setDissipation(0.25);
        }
    }
    //! [FCD_2D:aftertime]

    //! [FCD_2D:datamembers]
    int N = 1;
    LinearViscoelasticSpecies* FC2D_Species;

    //! [FCD_2D:datamembers]
};
//! [FCD_2D:class]

//! [FCD_2D:main]
int main()
{
    //! [FCD_2DproblemSetup]
    // Problem setup
    FreeCoolingDemoProblem problem;
    problem.setName("FreeCoolingDemo");
    problem.N=2000;
    problem.setGravity(Vec3D(0.0,0.0,0.0));
    problem.setTimeStep(5e-5);
    problem.setSaveCount(4000);
    problem.setTimeMax(100.0);
    problem.setMax({0.032,0.032,0.032});
    problem.setHGridMaxLevels(1);
    problem.setHGridCellOverSizeRatio(1.2);
    problem.setHGridUpdateEachTimeStep(false);
    //! [FCD_2DproblemSetup]
    //! [FCD_2D:solve]
    problem.setFileType(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    problem.solve();
    //! [FCD_2D:solve]
}
//![FCD_2D:main]

