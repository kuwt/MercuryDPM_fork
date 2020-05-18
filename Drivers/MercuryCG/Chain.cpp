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

#include <Walls/InfiniteWall.h>
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Species/LinearViscoelasticSpecies.h"

/** In this file a chain of 5 particles is settling under small gravity.
 * After that statistics are calculated.
 */
class Chain final : public Mercury3D
{
public:

    /**
     * Defines the problem setup
     */
    void setupInitialConditions() override
    {
        //define a particle species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());;
        species->setDensity(6.0 / constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(50e-4, 0.1, 1.0);


        //set gravity
        setGravity(Vec3D(0., 0., -1.));

        // set number and radius of particles
        unsigned n = 5;
        Mdouble r = 0.5;

        setXMax(r);
        setYMax(r);
        setZMax(n*2.0*r);
        setXMin(-r);
        setYMin(-r);
        setZMin(0);

        SphericalParticle p(species);
        p.setRadius(r);
        for (unsigned i = 0; i < n; ++i)
        {
            p.setPosition(Vec3D(0, 0, r*(1+2*i)));
            particleHandler.copyAndAddObject(p);
        }

        //set walls
        InfiniteWall w(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()),species);
        wallHandler.copyAndAddObject(w);
    }
};

int main()
{
    Chain chain;
    chain.setName("Chain");
    chain.setTimeStep(1e-4);
    chain.setTimeMax(1);
    chain.setSaveCount(250);
    //chain.restartFile.setName("ChainRenamed.restart");

    //define coarse-graining resolved in z
    CG<CGCoordinates::Z> cgZ;
    cgZ.setN(20);
    cgZ.setWidth(0.5);
    chain.cgHandler.copyAndAddObject(cgZ);

    chain.solve();

    logger(INFO,"Execute 'source Chain.sh' to get coarse-grained statistics of the last time step");
    helpers::writeToFile("Chain.sh","../MercuryCG/fstatistics Chain -stattype XZ -w 0.1 -h 0.05 -tmin 1 -tmax 30 -o Chain.XZ.stat\n"
            "../MercuryCG/MercuryCG Chain -stattype Z -w 0.1 -h 0.05 -tmin 1 -tmax 30 -o Chain.Z.stat\n"
            "../MercuryCG/MercuryCG Chain -stattype O -w 0.1 -h 0.05 -tmin 1 -tmax 30 -o Chain.O.stat\n");

    logger(INFO,"Run 'Chain.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("Chain.m","addpath('../MercuryCG/')\n"
            "data = loadstatistics('Chain.XZ.stat');\n"
            "colormap(1-gray)\n"
            "contourf(data.x,data.z,data.Density,20,'EdgeColor','none')\n"
            "c = colorbar\n"
            "c.Label.String = '\\rho';\n"
            "title('Density')\n"
            "xlabel('x')\n"
            "ylabel('z');\n"
            "axis equal\n"
            "%%\n"
            "particles=importdata('Chain.data',' ',12);\n"
            "z=particles.data(:,3);\n"
            "r=particles.data(:,7);\n"
            "a=linspace(0,2*pi,40);\n"
            "xCircle = sin(a);\n"
            "zCircle = cos(a);\n"
            "hold on;\n"
            "for i=1:length(x)\n"
            "  plot(x(i)+r(i)*xCircle,z(i)+r(i)*zCircle,'Color',.8*[1 1 1])\n"
            "end\n"
            "hold off");
    return 0;
}

