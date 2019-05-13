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

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
using constants::pi;

/*
* This is our problem description. Everything is set up here.
* We inherit from Mercury3D, since this gives full flexibility.
* For more predefined problems (for instance, chutes), please see the
* documentation.
*/
class Drum : public Mercury3D
{
public:
	/* We define our own 'setupInitialConditions' function here,
     * which defines all the specifics about our simulation here.
     */
	void setupInitialConditions() override
	{
		//Particle radiusSmall
		Mdouble radiusSmall = 0.01;
		Mdouble radiusLarge = 0.02;
		//Rotation frequency
		Mdouble frequency = 0.5; //Hz

		//The first step: set any properties which are always true for your system.
		// (for instance, if gravity remains constant, set it here)setName("Drum");
		setGravity(Vec3D(0,0,-9.8));
		setTimeMax(4);
		setTimeStep(0.0002);
		//visualised length
		setXMax(0.2);
		//visualised height
		setZMax(0.2);
		//visualised width
		setYMax(0.1);
		setXMin(-getXMax());
		setZMin(-getZMax());
		setYMin(-getYMax());

		//Now, decide what Species you need for your system.
		LinearViscoelasticSlidingFrictionSpecies species;
		species.setDensity(2000);
		Mdouble massSmall = species.getDensity() * 4.0 / 3.0 * pi * radiusSmall * radiusSmall * radiusSmall;
		species.setCollisionTimeAndRestitutionCoefficient(25*getTimeStep(),0.1,massSmall);
		species.setSlidingStiffness(2./7.*species.getStiffness());
		species.setSlidingDissipation(2./7.*species.getDissipation());
		species.setSlidingFrictionCoefficient(1.0);
		speciesHandler.copyAndAddObject(species);

		//place a drum wall
		AxisymmetricIntersectionOfWalls drum;
		drum.setSpecies(speciesHandler.getObject(0));
		drum.setPosition(Vec3D(0,0,0));
		drum.setAxis(Vec3D(0,1,0));
		Vec3D position = Vec3D(getXMax(),0,0);
		Vec3D normal = Vec3D(1,0,0);
		drum.addObject(normal,position);
		drum.setAngularVelocity(Vec3D(0,2.0 * pi*frequency,0));
		wallHandler.copyAndAddObject(drum);

		//Place particles into the box
		SphericalParticle particle;
		particle.setSpecies(speciesHandler.getObject(0));
		//Insert a set of particles (90 small and 10 large):
		while (particleHandler.getNumberOfObjects()<90)
		{
			if (particleHandler.getNumberOfObjects()<89) {
				particle.setRadius(radiusLarge*random.getRandomNumber(0.9,1.1));
			} else {
				particle.setRadius(radiusSmall*random.getRandomNumber(0.9,1.1));
			}
			//Insert a particle:
			do
			{
				Mdouble x = random.getRandomNumber(getXMin(), getXMax());
				Mdouble z = random.getRandomNumber(getZMin(), getZMax());
				particle.setPosition(Vec3D(x,0,z));
			} while (!checkParticleForInteraction(particle));
			particleHandler.copyAndAddObject(particle);
		}
	}

};

int main(int argc, char **argv)
{
	Drum problem;
	problem.setName("DrumSize");
	problem.setSaveCount(100);
	problem.setWallsWriteVTK(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
	//problem.setParticlesWriteVTK(true);
	problem.setXBallsAdditionalArguments(" -v0 -solidf");
	problem.restartFile.setFileType(FileType::NO_FILE);
	problem.fStatFile.setFileType(FileType::NO_FILE);

	//solve the system, the single particle will now bounce on the plate
	problem.solve(argc, argv);
	return 0;
}
