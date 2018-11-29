//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Walls/BasicIntersectionOfWalls.h>
#include <Walls/BasicUnionOfWalls.h>
#include "Mercury3D.h"
#include "Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h"
using constants::pi;
using mathsFunc::cubic;

class UnionOfWalls : public Mercury3D {

    void setupInitialConditions() override {
        setName("UnionOfWallsSelfTest");
        setTimeStep(1e-4);
        setTimeMax(1);
        setMin(-drumRadius*Vec3D(1,1,1));
        setMax(+drumRadius*Vec3D(1,1,1));
        setSaveCount(getTimeMax()/getTimeStep()/100);
        fStatFile.writeFirstAndLastTimeStep();
        restartFile.writeFirstAndLastTimeStep();
        setXBallsAdditionalArguments("-solidf");

        LinearPlasticViscoelasticSlidingFrictionSpecies species;
        species.setDensity(6.0/pi);
        species.setPlasticParameters(2e5,10e5,2e5,.5);
        species.setDissipation(25);
        species.setSlidingFrictionCoefficient(0.5);
        species.setSlidingDissipation(2./7.*species.getDissipation());
        species.setSlidingStiffness(2./7.*species.getLoadingStiffness());
        auto s = speciesHandler.copyAndAddObject(species);

        //define drum
        AxisymmetricIntersectionOfWalls drum;
        drum.setSpecies(s);
        drum.setPosition(Vec3D(0,0,0));
        drum.setAxis(Vec3D(0,1,0));
        drum.addObject(Vec3D(1,0,0),Vec3D(drumRadius,0,0));
        drum.setAngularVelocity(Vec3D(0,2.0*pi*100,0));
        drum.setVelocity(Vec3D(0,0,1));

        //define upper half
        InfiniteWall upperHalf;
        upperHalf.setSpecies(s);
        upperHalf.setPosition(Vec3D(0,0,0));
        upperHalf.set(Vec3D(0,0,1),Vec3D(0,0,0));

        //intersect drum and inflow
        BasicIntersectionOfWalls upperDrum;
        upperDrum.setSpecies(s);
        upperDrum.add(drum);
        upperDrum.add(upperHalf);

        //define lower half
        InfiniteWall lowerHalf;
        lowerHalf.setSpecies(s);
        lowerHalf.setPosition(Vec3D(0,0,0));
        lowerHalf.set(Vec3D(0,0,-1),Vec3D(0,0,0));

        //intersect drum and inflow
        BasicIntersectionOfWalls lowerDrum;
        lowerDrum.setSpecies(s);
        lowerDrum.add(drum);
        lowerDrum.add(lowerHalf);

        //union of lower and upper drum
        BasicUnionOfWalls unionDrum;
        unionDrum.setSpecies(s);
        unionDrum.add(lowerDrum);
        unionDrum.add(upperDrum);
        unionDrum.setAngularVelocity(Vec3D(0,2.0*pi*100,0));
        unionDrum.setVelocity(Vec3D(0,0,1));
        auto l = wallHandler.copyAndAddObject(unionDrum);

        BaseParticle particle;
        particle.setSpecies(s);
        particle.setRadius(0.5);
        particle.setPosition(Vec3D(0,0,-drumRadius+particle.getRadius()));
        //particle.setVelocity(Vec3D(0,0,100.0));
        auto p = particleHandler.copyAndAddObject(particle);
    }

    Mdouble drumRadius = 10.0;
};

int main()
{
    UnionOfWalls dpm;
    dpm.solve();
    return 0;
}
