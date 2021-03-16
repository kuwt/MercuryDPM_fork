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
#include "Species/LinearViscoelasticSpecies.h"

/** A chain of 5 particles is settling under small gravity.
 * This is a copy of the class in Chain.cpp.
 */
class Chain final : public Mercury3D
{
public:

    /**
     * Defines the problem setup
     */
    void setupInitialConditions() override
    {
        //simple properties
        setName("Chain");
        setTimeStep(1e-4);
        setTimeMax(1);
        setSaveCount(250);

        //define a particle species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
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

/**
 * A simple DPM problem is coarse-grained in a select region to show selective CG works
 * @return
 */
int main()
{
    logger(INFO,"\nRun a simple DPM problem\n");
    Chain chain;
    chain.solve();

    logger(INFO,"\nCoarse-grain a select region (1<z<2)\n");
    std::string cmd = "./MercuryCG Chain -tMin 1 -timeAverage -z 1 2 -averageBeyondDomain 0";
    if (system(cmd.c_str())==-1) {
        logger(WARN,"system call failed");
    }

    logger(INFO,"\nRead in CG\n");
    std::ifstream statFile("Chain.stat");
    //ignore 2 header lines and coordinate
    statFile.ignore(2000,'\n');
    statFile.ignore(2000,'\n');
    statFile.ignore(2000,' ');
    statFile.ignore(2000,' ');
    //read in volume fraction
    Mdouble density = 0;
    statFile >> density;
    logger(INFO,"Density % (should be 1)",density);
}

