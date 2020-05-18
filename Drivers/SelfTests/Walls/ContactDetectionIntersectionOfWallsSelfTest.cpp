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

#include <Walls/IntersectionOfWalls.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticReversibleAdhesiveSpecies.h>

/**
 * \brief Tests the contact detection between particles and IntersectionOfWalls.
 * \detail In particular, distinguishing face, edge and vertex contacts is tricky.
 * The most difficult case is when a face is less or equal in size to a particle, so this is tested here.
 **/
class ContactDetectionIntersectionOfWallsTest : public Mercury3D
{
public:
    void setupInitialConditions() override
    {
        setName("ContactDetectionIntersectionOfWallsSelfTest");
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::ONE_FILE);
        setTimeStep(1e-6);
        setTimeMax(getTimeStep());
        setMin(Vec3D(-.2,0,-.2));
        setMax(Vec3D(2.7,1,1.2));
        setXBallsAdditionalArguments("-w 1100 -s 0.9");

        LinearViscoelasticReversibleAdhesiveSpecies species;
        species.setDensity(1);
        species.setStiffness(1);
        species.setAdhesionStiffness(1);
        auto s0 = speciesHandler.copyAndAddObject(species);
        auto s1 = speciesHandler.copyAndAddObject(species);
        speciesHandler.getMixedObject(s0,s1)->setAdhesionForceMax(0.2);

        IntersectionOfWalls w;
        Vec3D minPos = Vec3D(0, 0, 0), maxPos = Vec3D(1, 1, 1);
        w.setSpecies(speciesHandler.getObject(1));
        w.addObject(Vec3D(1.0, 0.0, 0.0), minPos);
        w.addObject(Vec3D(0.0, 1.0, 0.0), minPos);
        w.addObject(Vec3D(0.0, 0.0, 1.0), minPos);
        w.addObject(Vec3D(-1.0, 0.0, 0.0), maxPos);
        w.addObject(Vec3D(0.0, -1.0, 0.0), maxPos);
        w.addObject(Vec3D(0.0, 0.0, -1.0), maxPos);
        w.addObject(Vec3D(1.0, 0.0, 1.0), Vec3D(0.9, 0.0, 0.0));
        wallHandler.copyAndAddObject(w);

        w.setPosition(Vec3D(1.5,0.6,0));
        wallHandler.copyAndAddObject(w);

        introduceParticlesAtWall();
    }

    void introduceParticlesAtWall()
    {

        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.01);

        Vec3D pos;
        Mdouble h = 2.0001*p.getRadius();
        pos.Y = 0.5;
        for (pos.X = getXMin(); pos.X <= getXMax(); pos.X += h)
        {
            for (pos.Z = getZMin(); pos.Z <= getZMax(); pos.Z += h)
            {
                Vec3D normal;
                Mdouble distance;
                p.setPosition(pos);
                bool touchesWall = false;
                for (auto w : wallHandler)
                {
                    //if touching the wall
                    if (w->getDistanceAndNormal(p, distance, normal))
                    {
                        touchesWall = true;
                        break;
                    }
                }
                if (touchesWall)
                {
                    particleHandler.copyAndAddObject(p);
                    logger(VERBOSE,"Inserted at %",p.getPosition());
                }
            }
        }
        std::cout << "Inserted particles: " << particleHandler.getNumberOfObjects() << std::endl;
    }

    void test()
    {
        //setWallsWriteVTK(FileType::ONE_FILE);
        //setParticlesWriteVTK(true);
        solve();
    }
};

int main(int argc, char* argv[])
{
    ContactDetectionIntersectionOfWallsTest().test();
    return 0;
}
