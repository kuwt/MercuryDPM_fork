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
class GetDistanceAndNormalForIntersectionOfWalls : public Mercury3D
{
public:
    void setupInitialConditions() override
    {
        setName("GetDistanceAndNormalForIntersectionOfWalls");
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::ONE_FILE);
        setTimeStep(1e-6);
        setTimeMax(getTimeStep());
        setMin(-.5*Vec3D(1,1,1));
        setMax(1.501*Vec3D(1,1,1));
        setXBallsAdditionalArguments("-w 1100 -s 0.9");

        //use an adhesive species to visualise negative overlaps
        LinearViscoelasticReversibleAdhesiveSpecies species;
        species.setDensity(1);
        species.setStiffness(1);
        species.setAdhesionStiffness(.1);
        species.setAdhesionForceMax(1);
        auto s = speciesHandler.copyAndAddObject(species);

        IntersectionOfWalls wall;
        wall.setSpecies(s);
        wall.addObject({0,0,1},{0,0,0});
        wall.addObject({1,0,0},{0,0,0});
        wall.addObject({-1,0,-1},{1,0,0});
        auto w = wallHandler.copyAndAddObject(wall);

        SphericalParticle particle;
        particle.setSpecies(s);
        particle.setRadius(0.2);
        auto p = particleHandler.copyAndAddObject(particle);

        Vec3D pos;
        Vec3D normal;
        Mdouble distance;
        Mdouble h = 0.01;
        std::ofstream file("GetDistanceAndNormalForIntersectionOfWalls.out");
        file << "x\tz\td\tnx\tnz\n";
        for (pos.X = getXMin(); pos.X <= getXMax(); pos.X += h) {
            for (pos.Z = getZMin(); pos.Z <= getZMax(); pos.Z += h) {
                p->setPosition(pos);
                w->getDistanceAndNormal(*p, distance, normal);
                file << pos.X
                     << '\t' << pos.Z
                     << '\t' << distance
                     << '\t' << normal.X
                     << '\t' << normal.Z
                     << '\n';
            }
        }
        file.close();

        helpers::writeToFile("GetDistanceAndNormalForIntersectionOfWalls.m",
                "clear variables\n"
                "raw=importdata('GetDistanceAndNormalForIntersectionOfWalls.out','\\t',1)\n"
                "n(1)=length(unique(raw.data(:,1)));\n"
                "n(2)=length(unique(raw.data(:,2)));\n"
                "x=reshape(raw.data(:,1),n);\n"
                "z=reshape(raw.data(:,2),n);\n"
                "d=reshape(raw.data(:,3),n);\n"
                "nx=reshape(raw.data(:,4),n);\n"
                "nz=reshape(raw.data(:,5),n);\n"
                "\n"
                "figure(1)\n"
                "contourf(x,z,d,20,'EdgeColor','none')\n"
                "hold on\n"
                "plot([0 0 1 0],[1 0 0 1],'k')\n"
                "hold off\n"
                "xlabel('x')\n"
                "ylabel('y')\n"
                "h=colorbar\n"
                "%xlabel(h, 'distance')\n"
                "title('distance')\n"
                "saveas(gcf,'distance.png')\n"
                "\n"
                "figure(2)\n"
                "m=8;\n"
                "quiver(x(1:m:end,1:m:end),z(1:m:end,1:m:end),nx(1:m:end,1:m:end),nz(1:m:end,1:m:end),0.5)\n"
                "hold on\n"
                "plot([0 0 1 0],[1 0 0 1],'k')\n"
                "hold off\n"
                "xlabel('x')\n"
                "ylabel('y')\n"
                "title('normal')\n"
                "saveas(gcf,'normal.png')\n"
                );
        logger(INFO,"Written file");
        exit(0);
    }
};

int main(int argc, char* argv[])
{
    GetDistanceAndNormalForIntersectionOfWalls().solve();
    return 0;
}
