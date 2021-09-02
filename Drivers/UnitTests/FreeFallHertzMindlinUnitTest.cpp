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
#include <Species/HertzianViscoelasticMindlinSpecies.h>

/// This case does a single elastic particle falling on an infinite plane.
/// The k is chosen so that the maximum overlap with the wall is around 2% of the particles diameter;
/// whereas, the time step must be taken to ensure 50 steps with a collision.
class FreeFallHertzMindlinUnitTest : public DPMBase
{
public:
	
	void setupInitialConditions() override {
		setMin(radius*Vec3D(-1,0,-1));
		setMax(radius*Vec3D(1,2,1));

		InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(0,-1,0), Vec3D(0, getYMin(), 0));
		wallHandler.copyAndAddObject(w0);
		
		SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
		p0.setPosition(Vec3D(0.0,2.0*radius,0.0));
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(radius);
		particleHandler.copyAndAddObject(p0);
	}

public:
	Mdouble radius = 5e-3;
};

int main(int argc, char* argv[])
{
	logger(INFO, "Single particle bouncing vertically on the bottom plate");

	//The monosized 3 mm perfect-spherical tungsten particles are simulated.
	//Poisson Ratio is 0.28,
	//bulk density is 19250 Kg/m^3,
	//Young modulus 4.11E11 Pa,
	//restitution coefficient is 0.95,
	//Hertz-Mindlin model is used to describe the contact force between neighboring particles.
	Mdouble poissonRatio = 0.28;
	Mdouble density = 10000*6.0/constants::pi;
	Mdouble elasticModulus = 4.11e11; //100 times too soft
	Mdouble restitution = 0.5;

	// Make the problem and set the name
	FreeFallHertzMindlinUnitTest dpm;
	dpm.setName("FreeFallHertzMindlinUnitTest");
	dpm.setGravity(Vec3D(0,-10,0));
	dpm.setParticleDimensions(3);

	//Set the species of the particle and wall, and its properties
	HertzianViscoelasticMindlinSpecies species;
	species.setDensity(density);
    species.setEffectiveElasticModulusAndRestitutionCoefficient(elasticModulus, restitution);
	//https://en.wikipedia.org/wiki/Shear_modulus#References
    species.setEffectiveShearModulus(0.5 * elasticModulus / (1 + poissonRatio));
	species.setSlidingFrictionCoefficient(1.0);
	dpm.speciesHandler.copyAndAddObject(species);

	//set the parameters for the solver
	dpm.setSaveCount(40);
	dpm.fStatFile.setFileType(FileType::NO_FILE);
	dpm.setWallsWriteVTK(FileType::ONE_FILE);
	Mdouble relativeVelocity = 0.1;
	Mdouble tc = species.getCollisionTime(2.0*dpm.radius, species.getDensity(), relativeVelocity);
	logger(INFO,"Collision time %", tc);
	dpm.setTimeStep(0.02*tc);
	dpm.setTimeMax(0.1);

	//solve the system, the single particle will now bounce on the plate
	dpm.solve(argc, argv);

	helpers::writeToFile("FreeFallHertzMindlinUnitTest.gnu",
						 "set xlabel 't'; p 'FreeFallHertzMindlinUnitTest.ene' u 1:(sqrt($3)) t 'sqrt(E_kin)' w lp");
	logger(INFO,"type 'gnuplot FreeFallHertzMindlinUnitTest.gnu' to check energy conservation");
}

