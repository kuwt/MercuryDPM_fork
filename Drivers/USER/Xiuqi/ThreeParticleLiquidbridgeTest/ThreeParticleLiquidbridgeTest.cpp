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

#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Particles/LiquidFilmParticle.h"
#include "Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h"

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class ThreeParticleElasticCollisionInteraction : public Mercury2D {

	void setupInitialConditions() override {
		setXMax(0.01);
        setYMax(0.01);
        setZMax(0.01);
		setGravity({0.0,0.0,0.0});

		LiquidFilmParticle P0,P1,P2;
		P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
        P2.setSpecies(speciesHandler.getObject(0));
		P0.setPosition(Vec3D(0.00485,0.005,0.005));
        P1.setPosition(Vec3D(0.005,0.005,0.005));
		P2.setPosition(Vec3D(0.00515,0.005,0.005));
        P0.setLiquidVolume(pow(10, -12));
        P1.setLiquidVolume(pow(10, -12));
        P2.setLiquidVolume(pow(10, -12));

		P0.setVelocity(Vec3D(0.0,0.0,0.0));
		P1.setVelocity(Vec3D( 0.0,0.0,0.0));
        P2.setVelocity(Vec3D( 0.0,0.0,0.0));
	
		P0.setRadius(0.0001);
		P1.setRadius(0.0001);
        P2.setRadius(0.0001);
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
        particleHandler.copyAndAddObject(P2);
		
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
	ThreeParticleElasticCollisionInteraction threeParticleElasticCollisionInteractionProblem;
    threeParticleElasticCollisionInteractionProblem.setName("ThreeParticleLiquidbridgeTest");

    auto species = new LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies;
    threeParticleElasticCollisionInteractionProblem.speciesHandler.addObject(species);
    species->setDensity(2650);
    species->setStiffness(1.5e4);
    species->setSlidingStiffness(0.4*species->getStiffness());
    species->setSlidingDissipation(0.1);
    species->setSlidingFrictionCoefficient(0.4);
    species->setSlidingFrictionCoefficientStatic(0.4);

    species->setLiquidBridgeVolumeMax(2.0*pow(10, -12));
    species->setLiquidBridgeVolumeMin(0.0);
    species->setDistributionCoefficient(0.0);
    species->setSurfaceTension(0.07);
    species->setContactAngle(constants::pi/18);
    species->setViscosity(0.00085);

    threeParticleElasticCollisionInteractionProblem.setTimeMax(0.02);
    threeParticleElasticCollisionInteractionProblem.setSaveCount(10);
    threeParticleElasticCollisionInteractionProblem.setTimeStep(1e-5);

    threeParticleElasticCollisionInteractionProblem.setFileType(FileType::ONE_FILE);
    threeParticleElasticCollisionInteractionProblem.setParticlesWriteVTK(true);
    threeParticleElasticCollisionInteractionProblem.setInteractionsWriteVTK(true);
    threeParticleElasticCollisionInteractionProblem.setWallsWriteVTK(FileType::ONE_FILE);

	threeParticleElasticCollisionInteractionProblem.solve();

}
