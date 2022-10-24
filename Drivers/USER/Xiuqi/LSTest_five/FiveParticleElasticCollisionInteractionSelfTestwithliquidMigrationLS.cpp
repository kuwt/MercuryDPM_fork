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
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Particles/LiquidFilmParticle.h"
#include "Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h"

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class FiveParticleElasticCollisionInteraction : public Mercury3D {
public:
    double radius_min, radius_max;

	void setupInitialConditions() override {
		setXMax(0.01);
        setYMax(0.01);
        setZMax(0.01);
		setGravity({0.0,0.0,-9.8});

        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies());
        species->setDensity(2650);
        species->setStiffness(1e4);
        species->setDissipation(0.002);
        species->setSlidingStiffness(0.4*species->getStiffness());
        species->setSlidingDissipation(0.4*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.1);
        species->setSlidingFrictionCoefficientStatic(0.1);

        Mdouble mass_min = 4.0/3.0 * constants::pi * pow(radius_min, 3.0) * species->getDensity();
        species->setLiquidBridgeVolumeMax(0.05 * mass_min/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.0);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);

		LiquidFilmParticle p;
		p.setSpecies(speciesHandler.getObject(0));

        for (unsigned int i = 0; i < 5; ++i){
            const double x = getXMax()/6.0 * (i + 1);
            const double z = getZMax()/6.0 * (i + 1);
            p.setPosition(Vec3D(x,0.0, z));
            p.setVelocity(Vec3D(random.getRandomNumber(-0.0005,0.0005),  0.0, 0.0));
            p.setRadius(random.getRandomNumber(radius_min, radius_max));
            Mdouble mass = 4.0/3.0 * constants::pi * pow(p.getRadius(), 3.0) * species->getDensity();
            p.setLiquidVolume(0.05 * mass/1000.0);
            particleHandler.copyAndAddObject(p);
        }


		wallHandler.clear();
		InfiniteWall w0;
		w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(-1, 0, 0), getMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0), getMax());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), getMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0), getMax());
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 0, 0, -1), getMin());
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 0, 0, 1), getMax());
		wallHandler.copyAndAddObject(w0);		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	FiveParticleElasticCollisionInteraction fiveParticleElasticCollisionInteractionProblem;
    fiveParticleElasticCollisionInteractionProblem.setName("FiveParticleElasticCollisionInteractionwithLiquidMigrationLSSelfTest");

    fiveParticleElasticCollisionInteractionProblem.radius_min = 0.0002;
    fiveParticleElasticCollisionInteractionProblem.radius_max = 0.0005;
    fiveParticleElasticCollisionInteractionProblem.setTimeMax(0.05);
    fiveParticleElasticCollisionInteractionProblem.setSaveCount(10);
    fiveParticleElasticCollisionInteractionProblem.setTimeStep(1e-5);
    fiveParticleElasticCollisionInteractionProblem.setFileType(FileType::ONE_FILE);
    //fiveParticleElasticCollisionInteractionProblem.getInteractionFile().setFileType(FileType::ONE_FILE);
    fiveParticleElasticCollisionInteractionProblem.setParticlesWriteVTK(true);
    fiveParticleElasticCollisionInteractionProblem.setWallsWriteVTK(FileType::ONE_FILE);


    fiveParticleElasticCollisionInteractionProblem.solve();

}
