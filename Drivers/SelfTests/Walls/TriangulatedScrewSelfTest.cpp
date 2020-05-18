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

#include "Mercury3D.h"
#include "Walls/TriangulatedWall.h"
#include "Walls/Screw.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Logger.h"

/*!
 * \brief Tests the implementation of TriangulatedWall.
 * \details
 *  <img src="Walls/triangulatedWallsSelfTest.png" height="250px">
 */
class TriangulatedScrewSelfTest : public Mercury3D
{
public:
    TriangulatedScrewSelfTest ()
    {
        setName("TriangulatedScrewSelfTest");
        setMin({-1,-1,-1});
        setMax({1,1,1});

        setGravity(30*Vec3D(0,0,1));
        setTimeStep(1e-4);
        setTimeMax(1e3*getTimeStep()); //use 2.5e4 to get reasonable output
        //setTimeMax(getTimeStep());
        setSaveCount(100);
        setXBallsAdditionalArguments("-v0 -solidf");
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.setSaveCount(1e5);

        auto s = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        s->setDensity(1);
        s->setCollisionTimeAndRestitutionCoefficient(20.0*getTimeStep(),0.2,s->getMassFromRadius(radius));
        s->setSlidingFrictionCoefficient(0.1);
        s->setSlidingStiffness(2.0/7.0*s->getStiffness());
        s->setSlidingDissipation(2.0/7.0*s->getDissipation());
    }

    void setupInitialConditions() override
    {
        AxisymmetricIntersectionOfWalls a({0,0,0},{0,0,1},{},speciesHandler.getLastObject());
        a.addObject({1,0,0},{1,0,0});
        wallHandler.copyAndAddObject(a);

        PeriodicBoundary b;
        b.set({0,0,-1},getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(b);

        Mdouble extra = 1.2;
        Vec3D start = {0,0,-1*extra};
        Mdouble length = 2 * extra;
        Mdouble rad = 1;
        Mdouble nRotations = 1 * extra;
        Mdouble thickness = 0.1;
        Screw s(start, length, rad, nRotations, 0, thickness);
        s.writeVTK("Screw.vtk");
        TriangulatedWall w("Screw.vtk",speciesHandler.getLastObject());
        Mdouble omega = 0;
        w.setAngularVelocity({0,0,omega});
        w.setVelocity({0,0,length*omega});
        wallHandler.copyAndAddObject(w);

        //introduce particles randomly
        SphericalParticle p;
        p.setSpecies(speciesHandler.getLastObject());

        unsigned n = 1000;

        if (n==1) {
            p.setPosition({.5,.5,.9});
            p.setRadius(radius);
            particleHandler.copyAndAddObject(p);
        } else {
            for (unsigned i = 0; i < n; i++)
            {
                Vec3D r;
                r.X = random.getRandomNumber(getXMin(), getXMax());
                r.Y = random.getRandomNumber(getYMin(), getYMax());
                r.Z = random.getRandomNumber(getZMin(), getZMax());
                p.setPosition(r);
                p.setRadius(random.getRandomNumber(radius, 1.1 * radius));
                if (checkParticleForInteraction(p))
                    particleHandler.copyAndAddObject(p);
            }
        }

        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
    }

    Mdouble radius = 0.15;
};

int main()
{
    logger(INFO,"The simulated time interval is very small to reduce the length of this self test;"
     " to get nice output, increase timeMax to ~3 seconds,"
     " and add VTK output (uncomment the relevant lines in the main function).");

    TriangulatedScrewSelfTest t;
    //uncomment to get VTK output
    //t.setParticlesWriteVTK(true);
    t.setWallsWriteVTK(FileType::ONE_FILE);
    t.solve();
    return 0;
}
