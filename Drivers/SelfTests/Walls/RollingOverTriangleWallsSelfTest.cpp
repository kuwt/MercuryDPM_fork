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

#include <Walls/Coil.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticBondedSpecies.h>
#include <Walls/TriangleWall.h>

/**
 * \brief Tests the contact detection between particles and a set of TriangleWalls.
 * \detail In particular, distinguishing face, edge and vertex contacts is tricky.
 * So here a particle is set to rollover a face, edge and vertex of a flat wall made from particles.
 **/
class RollingOverTriangleWalls : public Mercury3D
{
public:
    void setupInitialConditions() override
    {
        setName("RollingOverTriangleWalls");
        removeOldFiles();
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        setWallsWriteVTK(true);
        setParticlesWriteVTK(true);
        setSaveCount(2);

        setTimeStep(1e-3);
        setTimeMax(2.0);
        setGravity({0,0,-1});
        setMin(-1 * Vec3D(1, 1, 1));
        setMax(+1 * Vec3D(1, 1, 1));

        //use an adhesive species to visualise negative overlaps
        LinearViscoelasticBondedSpecies species;
        species.setDensity(6.0/constants::pi);
        species.setStiffness(2e3);
        species.setDissipation(0);
        species.setBondForceMax(0);
        auto s = speciesHandler.copyAndAddObject(species);

        //eight triangles forming a plane
        TriangleWall wall;
        wall.setSpecies(s);
        wall.setGroupId(1);
        double x[3]{-.7, 0, .7};
        double y = .9;
        double z[3]{-.7, 0, .7};
        /*
         *  1---4-----7    y
         *  |  /|  /  |    |
         *  | / 3-----6    |___x
         *  |/  |  /  |
         *  0---2-----5
         */
        Vec3D P[8] {
                {0,-1,0},
                {0,1,0},
                {1,-1,0},
                {1,0,0},
                {1,1,0},
                {2,-1,0},
                {2,0,0},
                {2,1,0},
        };

        wall.setVertices(P[0], P[1], P[4]);
        wallHandler.copyAndAddObject(wall);
        wall.setVertices(P[0], P[2], P[4]);
        wallHandler.copyAndAddObject(wall);
        wall.setVertices(P[2], P[3], P[6]);
        wallHandler.copyAndAddObject(wall);
        wall.setVertices(P[2], P[5], P[6]);
        wallHandler.copyAndAddObject(wall);
        wall.setVertices(P[3], P[4], P[7]);
        wallHandler.copyAndAddObject(wall);
        wall.setVertices(P[3], P[6], P[7]);
        wallHandler.copyAndAddObject(wall);

        //place a particle on top of the plane
        SphericalParticle particle;
        particle.setSpecies(s);
        particle.setRadius(0.5);
        double eq = particle.getMass() * -getGravity().getZ()/species.getStiffness();
        particle.setPosition({2,-0.5,0.5-eq});
        particle.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(particle);
    }

};

int main(int argc, char* argv[])
{
    RollingOverTriangleWalls().solve();
    logger(INFO,"Run \"p [0:2]  'RollingOverTriangleWalls.fstat' u 4:9\" to view force vs position");
    return 0;
}
