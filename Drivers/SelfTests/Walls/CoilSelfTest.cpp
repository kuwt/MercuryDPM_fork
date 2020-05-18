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

//! [CST:headers]
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/Coil.h"
#include "Species/LinearViscoelasticSpecies.h"
//! [CST:headers]

//! [CST:class]
/*!
 * The CoilSelfTest is a class that tests the coil, by simulating the flow around a coil for a short time and then
 * comparing the restart file with an old restart file. If they are the same, the coil still works correctly.
 * This class can also be used to see how the coil works, and how it should be used in a bigger application.
 * This code is one of the advanced tutorials, therefore do not remove the ///![CST:xxx] lines.
 * \todo update documentation to simpler structure
 */
class CoilSelfTest : public Mercury3D
{
private:
    void setupInitialConditions() override
    {
        //! [CST:walls]
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        //front wall
        w.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w);

        //back wall
        w.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w);

        //bottom wall
        w.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w);

        //top wall
        w.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w);

        //left wall
        w.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w);

        //right wall
        InfiniteWallWithHole rightWall;
        rightWall.setSpecies(speciesHandler.getObject(0));
        rightWall.set(Vec3D(0, 0, 1), getZMax(), 1.0);
        wallHandler.copyAndAddObject(rightWall);
        //! [CST:walls]

        //! [CST:coil]
        // creation of the coil and setting of its properties
        Coil coil;
        coil.setSpecies(speciesHandler.getObject(0));
        // the syntax to set the coil geometry is as follows:
        // set(Start position, Length, Radius, Number of turns, Rotation speed, Thickness)
        coil.set(Vec3D(0, 0, 0), 1.0, 1.0 - particleRadius, 2.0, -1.0, 0.5 * particleRadius);
        auto pCoil = wallHandler.copyAndAddObject(coil);
        //! [CST:coil]

        //! [CST:particle]
        particleHandler.clear();
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(particleRadius);
        //! [CST:particle]
        /*
        //Single test case
        double distance;
                Vec3D normal;
        p0.setPosition(Vec3D(1.0,0.0,0.0));
        if(coil->getDistance_and_normal(p0, distance, normal))
                        std::cout<<"Collision, distance screw="<<distance<<std::endl;
                else
                        std::cout<<"No collision, distance screw="<<distance<<std::endl;
         */

        //! [CST:placeparticles]
        // Nx*Ny*Nz particles are created and placed on evenly spaced positions in
        // the domain [xMin_,xMax_]*[yMin_,yMax_]*[zMin_,zMax_] (unless their position
        // is already occupied by the coil).
        const unsigned int Nx = static_cast<unsigned int> (std::floor((getXMax() - getXMin()) / (2.1 * particleRadius)));
        const unsigned int Ny = static_cast<unsigned int> (std::floor((getYMax() - getYMin()) / (2.1 * particleRadius)));
        const unsigned int Nz = static_cast<unsigned int> (std::floor((getZMax() - getZMin()) / (2.1 * particleRadius)));
        Mdouble distance;
        Vec3D normal;

        for (unsigned int i = 0; i < Nx; i++)
        {
            for (unsigned int j = 0; j < Ny; j++)
            {
                for (unsigned int k = 0; k < Nz; k++)
                {
                    p0.setPosition(Vec3D(getXMin() + (getXMax() - getXMin()) * (0.5 + i) / Nx,
                                         getYMin() + (getYMax() - getYMin()) * (0.5 + j) / Ny,
                                         getZMin() + (getZMax() - getZMin()) * (0.5 + k) / Nz));
                    if (!pCoil->getDistanceAndNormal(p0, distance, normal)) //if there is no collision with the coil
                    {
                        particleHandler.copyAndAddObject(p0);
                    }
                    else
                    {
                        logger(DEBUG, "particle at position % could not be inserted", p0.getPosition());
                    }
                }
            }
        }
        //! [CST:placeparticles]
    }
    //! [CST:beforetime]
    ///After t = 1, the coil turns every time step.
    void actionsBeforeTimeStep() override
    {
        if (getTime() > 1)
        {
            coil->move_time(getTimeStep());
        }
    }
    //! [CST:beforetime]

public:
    ///\todo why is the coil public?
    // ! [CST:datamembers]
    Coil* coil;
    Mdouble particleRadius;
    //! [CST:datamembers]
};
//! [CST:class]

//! [CST:main]
int main(int argc UNUSED, char* argv[] UNUSED)
{
    
    // create CoilSelfTest object
    CoilSelfTest problem;

    // set some basic problem properties
    //! [CSTproblemSetup]
    problem.setName("CoilSelfTest");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0, -9.8, 0.0));

    // set problem geometry
    problem.setXMax(1.0);
    problem.setYMax(5.0);
    problem.setZMax(2.0);
    problem.setXMin(-1.0);
    problem.setYMin(-1.0);
    problem.setTimeMax(0.5);
    //! [CSTproblemSetup]
    problem.particleRadius = 0.2;

    //! [CST:species]
    LinearViscoelasticSpecies species;
    species.setDensity(1000);
    Mdouble tc = 0.05;
    Mdouble restitutionCoefficient = 0.8;

    Mdouble particleMass = pow(problem.particleRadius, 3) * constants::pi * 4.0 / 3.0 * species.getDensity();
    species.setCollisionTimeAndRestitutionCoefficient(tc, restitutionCoefficient, particleMass);
    problem.speciesHandler.copyAndAddObject(species);
    //! [CST:species]

    //! [CST:solve]
    problem.setTimeStep(0.02 * 0.05);
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000, problem.getTimeMax(),
                                                                                     problem.getTimeStep()));
    problem.solve();
    //! [CST:solve]
    
}
//! [CST:main]
