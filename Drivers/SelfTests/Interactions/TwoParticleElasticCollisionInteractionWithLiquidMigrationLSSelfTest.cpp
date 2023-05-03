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
#include "Species/Species.h"
#include "Particles/LiquidFilmParticle.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationLSSpecies.h"

/// In this file two particles are symmetrically placed in a domain with opposite initial velocities, to test the liquid bridge force in normal
/// and tangential direction after colliding obliquely.

class TwoParticleElasticCollisionInteraction : public Mercury2D {
public:
    double radius;

	void setupInitialConditions() override {
		setXMax(0.01);
        setYMax(0.01);
        setZMax(0.01);
		setGravity({0.0,0.0,0.0});

        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationLSSpecies());
        species->setDensity(2650);
        species->setStiffness(1e4);
        species->setDissipation(0.002);
        species->setSlidingStiffness(0.4*species->getStiffness());
        species->setSlidingDissipation(0.4*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.4);
        species->setSlidingFrictionCoefficientStatic(0.4);
        Mdouble mass_species = species->getMassFromRadius(radius);
        //define the moisture level in mass ratio
        species->setLiquidBridgeVolumeMax(mass_species*0.05*2/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.0);
        //parameters in LiquidMigrationLSSpecies
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);

		LiquidFilmParticle P0,P1;
		P0.setSpecies(speciesHandler.getObject(0));
		P1.setSpecies(speciesHandler.getObject(0));
		P0.setPosition(Vec3D(0.007,0.005,0.005));
		P1.setPosition(Vec3D(0.003,0.005,0.005));
        P0.setRadius(radius);
        P1.setRadius(radius);

        Mdouble mass = species->getMassFromRadius(radius);
        logger(INFO, "Mass = %", mass);
        P0.setLiquidVolume(mass*0.05/1000.0);
        P1.setLiquidVolume(mass*0.05/1000.0);
	
		P0.setVelocity(Vec3D(-0.1,0.005,0.0));
		P1.setVelocity(Vec3D( 0.1,-0.005,0.0));

		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
		
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

    //test the force acting on particles when they get in contact and move apart
    void printTime() const override {
        Vec3D NormalForce = Vec3D(0,0,0);
        Vec3D TangentialForce = Vec3D(0,0,0);
        for (BaseInteraction * const i : interactionHandler) {
            SlidingFrictionInteraction *j = dynamic_cast<SlidingFrictionInteraction *>(i);
            NormalForce += i->getForce()-j->getTangentialForce();
            TangentialForce += j->getTangentialForce();
        }
        logger(INFO, "t = %", getTime());
        logger(INFO, "Normal force = %, Tangential force = %", NormalForce, TangentialForce);
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	TwoParticleElasticCollisionInteraction twoParticleElasticCollisionInteractionProblem;
    twoParticleElasticCollisionInteractionProblem.setName("TwoParticleElasticCollisionInteractionwithLiquidMigrationLSSelfTest");
    twoParticleElasticCollisionInteractionProblem.radius = 0.0002;

	twoParticleElasticCollisionInteractionProblem.setTimeMax(0.05);
	twoParticleElasticCollisionInteractionProblem.setSaveCount(10);
    twoParticleElasticCollisionInteractionProblem.setTimeStep(1e-5);
	twoParticleElasticCollisionInteractionProblem.fStatFile.setFileType(FileType::ONE_FILE);
	twoParticleElasticCollisionInteractionProblem.getInteractionFile().setFileType(FileType::ONE_FILE);
    twoParticleElasticCollisionInteractionProblem.setParticlesWriteVTK(true);
    twoParticleElasticCollisionInteractionProblem.interactionHandler.setWriteVTK(true);
    twoParticleElasticCollisionInteractionProblem.wallHandler.setWriteVTK(FileType::ONE_FILE);

	twoParticleElasticCollisionInteractionProblem.solve();

}