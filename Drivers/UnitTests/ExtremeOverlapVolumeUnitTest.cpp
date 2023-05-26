//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Logger.h>

/*!
 * \brief This test checks the formula for computing the overlap volume between two spherical particles.
 * \details In this test, two particles are initialized and the position of the smaller particles is updated multiple
 * times, to then calculate and compare the overlap volume. The first and fourth overlap volumes are pre-calculated,
 * the second and third overlap volumes should equal half and fully the smaller particle volume.
 */
class ExtremeOverlapVolumeUnitTest : public DPMBase{
public:
	void setupInitialConditions() final
	{
        //Define the radii of the two particles
		Mdouble r0=0.01;
        Mdouble r1=0.001;

        //Define the bounding box (for viewing the output only; no walls are created)
		setYMax(2.0*(r0+r1));
		setYMin(0.0);
		setZMax(r0);
		setZMin(-r0);
		setXMax(r0);
		setXMin(-r0);

        //Create both particles and set their positions, radii, and velocities.
		SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
		p.setPosition(Vec3D(0.0,r0,0.0));
		p.setRadius(r0);
		particleHandler.copyAndAddObject(p);
		p.setRadius(r1);
		particleHandler.copyAndAddObject(p);

        // Small overlap. Should equal pre-calculated value.
        particleHandler.getObject(1)->setPosition(Vec3D(0.0,2.0*r0+0.9*r1,0.0));
        computeAllForces();
        helpers::check(interactionHandler.getObject(0)->getOverlapVolume(), 2.77675e-11, 1.0e-12,"overlapVolume");

        // Small particle center exactly on big particle surface. Overlap volume should equal half volume smallest particle.
        particleHandler.getObject(1)->setPosition(Vec3D(0.0,2.0*r0,0.0));
        computeAllForces();
        helpers::check(interactionHandler.getObject(0)->getOverlapVolume(), 0.5 * particleHandler.getObject(1)->getVolume(), 1.0e-10,"overlapVolume");

        // Small particle fully inside big particle, but still just touching big particle surface. Overlap volume should equal volume smallest particle.
        particleHandler.getObject(1)->setPosition(Vec3D(0.0,2.0*r0-r1,0.0));
        computeAllForces();
        helpers::check(interactionHandler.getObject(0)->getOverlapVolume(), particleHandler.getObject(1)->getVolume(), 1.0e-12,"overlapVolume");

        // Small particle fully inside big particle, and not touching the big particle surface anymore, causes the overlap volume calculations to return nonsensical values.
        // It should, however, equal to pre-calculated value.
        particleHandler.getObject(1)->setPosition(Vec3D(0.0,2.0*r0-1.1*r1,0.0));
        computeAllForces();
        helpers::check(interactionHandler.getObject(0)->getOverlapVolume(), 4.15244e-09, 1.0e-12,"overlapVolume");

        exit(0);
	}
};

int main(int argc, char *argv[])
{
	ExtremeOverlapVolumeUnitTest OverlapProblem;

    auto species = OverlapProblem.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    species->setDensity(2000);
    species->setStiffness(1e1);
    species->setDissipation(0.0);
    species->setSlidingFrictionCoefficient(1.0);
    species->setSlidingStiffness(2.0/7.0*species->getStiffness());
    species->setSlidingDissipation(2.0/7.0*species->getDissipation());
    
    OverlapProblem.setName("ExtremeOverlapUnitTest");
	OverlapProblem.solve(argc,argv);

    return 0;
}
