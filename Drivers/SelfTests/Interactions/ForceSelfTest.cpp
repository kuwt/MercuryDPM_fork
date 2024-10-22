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
#include "Species/Species.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
//#include "Interactions/SlidingFrictionInteraction.h"
//#include "Interactions/LinearViscoelasticInteraction.h"

class ParticleParticleCollision : public DPMBase{
public:

	ParticleParticleCollision()
	{
		///\bug setting gravity here was a quick fix when I changed the default gravity to zero. Should gravity be 0 here?
		setGravity(Vec3D(0,-9.8,0));
		fStatFile.setFileType(FileType::NO_FILE);
		species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
	}

	void setupInitialConditions() override {
		static int count = -1;
		count++;

		// only do this for the first time
		if (!count) {
			setSystemDimensions(2);
			setParticleDimensions(3);
			species->setDensity(2000.0);

			unsigned int nParticle = 2;
			normal = Vec3D(-1.0,0.0,0.0);
			tangent = Vec3D(0.0,1.0,0.0);

			SphericalParticle P0;
			P0.setSpecies(speciesHandler.getObject(0));
			for (unsigned int j=0; j< nParticle; j++){
				P0.setRadius(random.getRandomNumber(0.0005,0.001));
				particleHandler.copyAndAddObject(P0);
				//P0.getSpecies()->computeMass(P0&);

			}

			//set random relative velocity
			initialNormalRelativeVelocity = random.getRandomNumber(-0.1,-0.05);
			initialTangentialRelativeVelocity = random.getRandomNumber( 0.05, 0.1);

			tc = random.getRandomNumber(0.0,1.0e-5);
			en = random.getRandomNumber(0.5,1.0);

			setXMax(0.005);
			setYMax(0.005);
			setZMax(0.005);

			m12 = (particleHandler.getObject(0)->getMass()*particleHandler.getObject(1)->getMass())/(particleHandler.getObject(0)->getMass()+particleHandler.getObject(1)->getMass());
			species->setCollisionTimeAndRestitutionCoefficient(tc,en,2.0*m12);
			setTimeMax(2.0*tc);
			setTimeStep(tc / 200);
			setSaveCount(1);
		}

		// do this for all solves
		particleHandler.getObject(0)->setVelocity(normal*initialNormalRelativeVelocity+tangent*initialTangentialRelativeVelocity);
		particleHandler.getObject(1)->setVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,0.0,0.0));

		particleHandler.getObject(0)->setPosition(Vec3D(0.0025-particleHandler.getObject(0)->getRadius(),0.0025,0.0025));
		particleHandler.getObject(1)->setPosition(Vec3D(0.0025+particleHandler.getObject(1)->getRadius(),0.0025,0.0025));
		particleHandler.getObject(0)->setOrientation({1,0,0,0});
		particleHandler.getObject(1)->setOrientation({1,0,0,0});

	}

	Vec3D getRelativeVelocity() {
		return particleHandler.getObject(0)->getVelocity() - particleHandler.getObject(1)->getVelocity() + Vec3D::cross(normal,particleHandler.getObject(0)->getAngularVelocity()*particleHandler.getObject(0)->getRadius() + particleHandler.getObject(1)->getAngularVelocity()*particleHandler.getObject(1)->getRadius());
	}

	void getRelativeVelocityComponents(Mdouble& normalRelativeVelocity, Mdouble& tangentialRelativeVelocity) {
		Vec3D relativeVelocity = getRelativeVelocity();
		normalRelativeVelocity = relativeVelocity.X / normal.X;
		tangentialRelativeVelocity = relativeVelocity.Y / tangent.Y;
	}

	///Calculates collision time for two copies of a particle of species 0
	Mdouble getCollisionTime(Mdouble mass){
		return species->getCollisionTime(mass);
	}


	Mdouble initialNormalRelativeVelocity, initialTangentialRelativeVelocity;
	Vec3D normal, tangent;
	Mdouble tc, en;
	Mdouble m12;
	LinearViscoelasticSlidingFrictionSpecies* species;
};

class WallParticleCollision : public ParticleParticleCollision{
public:

	void setupInitialConditions() override {
		static int count = -1;
		count++;

		// only do this for the first time
		if (!count) {
			setSystemDimensions(2);
			setParticleDimensions(3);
			species->setDensity(2000.0);


			InfiniteWall w0;
			w0.setSpecies(speciesHandler.getObject(0));
			w0.set(Vec3D(1, 0, 0), Vec3D(0.0025, 0, 0));
			wallHandler.copyAndAddObject(w0);

			SphericalParticle P0;
			P0.setSpecies(speciesHandler.getObject(0));
			particleHandler.copyAndAddObject(P0);
			setXMax(0.005);
			setYMax(0.005);
			setZMax(0.005);

			//set random masses

			particleHandler.getObject(0)->setRadius(random.getRandomNumber(0.0005,0.001));
			particleHandler.getObject(0)->getSpecies()->computeMass(particleHandler.getObject(0));

			//set random relative velocity
			normal = Vec3D(-1.0,0.0,0.0);
			tangent = Vec3D(0.0,1.0,0.0);
			initialNormalRelativeVelocity = random.getRandomNumber(-0.1,-0.05);
			initialTangentialRelativeVelocity = random.getRandomNumber( 0.05, 0.1);

			tc = random.getRandomNumber(0.0,1.0e-5);
			en = random.getRandomNumber(0.5,1.0);
			m12 = particleHandler.getObject(0)->getMass(); // wall counts as an infinite mass
			species->setCollisionTimeAndRestitutionCoefficient(tc,en,2.0*m12);
			setTimeMax(2.0*tc);
			setTimeStep(tc / 200);
			setSaveCount(1);
		}

		// do this for all solves
		particleHandler.getObject(0)->setVelocity(normal*initialNormalRelativeVelocity+tangent*initialTangentialRelativeVelocity);
		particleHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,0.0,0.0));

		particleHandler.getObject(0)->setPosition(Vec3D(0.0025-particleHandler.getObject(0)->getRadius(),0.0025,0.0025));
		particleHandler.getObject(0)->setOrientation({1,0,0,0});
	}

	Vec3D getRelativeVelocity() {
		return particleHandler.getObject(0)->getVelocity() + Vec3D::cross(normal,particleHandler.getObject(0)->getAngularVelocity()*particleHandler.getObject(0)->getRadius());
	}
	void getRelativeVelocityComponents(Mdouble& normalRelativeVelocity, Mdouble& tangentialRelativeVelocity) {
		Vec3D relativeVelocity = getRelativeVelocity();
		normalRelativeVelocity = relativeVelocity.X / normal.X;
		tangentialRelativeVelocity = relativeVelocity.Y / tangent.Y;
	}

};

void particleParticleTest()
{
    logger(INFO, "\ntesting particle-particle collisions ...\n\n", Flusher::NO_FLUSH);
    
    ParticleParticleCollision problem;
    
    problem.random.setRandomSeed(5);
    
    problem.setupInitialConditions();
    Mdouble normalRelativeVelocity, tangentialRelativeVelocity, analyticTangentialRelativeVelocity;
    
    logger(INFO, "5: without tangential forces");
    problem.setName("ForceSelfTest5");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           problem.initialTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    //problem.setAppend_to_files(true);
    
    logger(INFO, "6: with Coulomb friction");
    Mdouble mu = problem.random.getRandomNumber(0.0, 1.0);
    problem.species->setSlidingFrictionCoefficient(mu);
    problem.species->setSlidingDissipation(1e20);
    //problem.species->setSlidingStiffness(0.0);
    problem.setName("ForceSelfTest6");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = std::max(0.0, problem.initialTangentialRelativeVelocity +
                                                       mu * 3.5 * (1 + problem.en) *
                                                       problem.initialNormalRelativeVelocity);
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           analyticTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "7: with Coulomb friction, spring activated");
    problem.species->setSlidingStiffness(1.0);
    //problem.species->setSlidingDissipation(1);
    problem.setName("ForceSelfTest7");
    problem.solve();
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           analyticTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "8: with tangential viscous force");
    Mdouble et = problem.random.getRandomNumber(-1.0, 0.0);
    problem.species->setSlidingFrictionCoefficient(1e20);
    problem.species->setSlidingDissipation(-log(-et) / (2.0 * problem.tc) / 3.5 * 2.0 * problem.m12);
    problem.species->setSlidingStiffness(0.0);
    problem.setName("ForceSelfTest8");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = problem.initialTangentialRelativeVelocity *
                                         exp(-2.0 * 3.5 * problem.species->getSlidingDissipation() /
                                             (2.0 * problem.m12) * problem.getCollisionTime(2.0 * problem.m12));
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           analyticTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "9: with tangential elastic force");
    Mdouble et2 = problem.random.getRandomNumber(0.0, 1.0);
    problem.species->setSlidingFrictionCoefficient(1e20);
    problem.species->setSlidingDissipation(0.0);
    problem.species->setSlidingStiffness(
            problem.species->getStiffness() / 3.5 * mathsFunc::square(acos(-et2) / constants::pi));
    problem.setName("ForceSelfTest9");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = problem.initialTangentialRelativeVelocity *
                                         cos(sqrt(problem.species->getSlidingStiffness() / problem.m12 * 3.5) *
                                             problem.getCollisionTime(2.0 * problem.m12));
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           analyticTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
}

void wallParticleTest()
{
    logger(INFO, "\ntesting wall-particle collisions ...\n");
    
    srand(5);
    
    WallParticleCollision problem;
    problem.setupInitialConditions();
    
    Mdouble normalRelativeVelocity, tangentialRelativeVelocity, analyticTangentialRelativeVelocity;
    
    //problem.setAppend_to_files(true);
    
    logger(INFO, "0: without tangential forces");
    problem.setName("ForceSelfTest0");
    problem.solve();
    
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           problem.initialTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "1: with Coulomb friction");
    Mdouble mu = problem.random.getRandomNumber(0.0, 1.0);
    problem.species->setSlidingFrictionCoefficient(mu);
    problem.species->setSlidingDissipation(1e20);
    problem.species->setSlidingStiffness(0.0);
    problem.setName("ForceSelfTest1");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = std::max(0.0, problem.initialTangentialRelativeVelocity +
                                                       mu * 3.5 * (1 + problem.en) *
                                                       problem.initialNormalRelativeVelocity);
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           analyticTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "2: with Coulomb friction, spring activated");
    problem.species->setSlidingStiffness(1.0);
    problem.setName("ForceSelfTest2");
    problem.solve();
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           problem.initialTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "3: with tangential viscous force");
    Mdouble et = problem.random.getRandomNumber(-1.0, 0.0);
    problem.species->setSlidingFrictionCoefficient(1e20);
    problem.species->setSlidingDissipation(-log(-et) / (2.0 * problem.tc) / 3.5 * 2.0 * problem.m12);
    problem.species->setSlidingStiffness(0.0);
    problem.setName("ForceSelfTest3");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = problem.initialTangentialRelativeVelocity *
                                         exp(-2.0 * 3.5 * problem.species->getSlidingDissipation() /
                                             (2.0 * problem.m12) * problem.getCollisionTime(2.0 * problem.m12));
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           problem.initialTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
    
    logger(INFO, "4: with tangential elastic force");
    Mdouble et2 = problem.random.getRandomNumber(0.0, 1.0);
    problem.species->setSlidingFrictionCoefficient(1e20);
    problem.species->setSlidingDissipation(0.0);
    problem.species->setSlidingStiffness(
            problem.species->getStiffness() / 3.5 * mathsFunc::square(acos(-et2) / constants::pi));
    problem.setName("ForceSelfTest4");
    problem.solve();
    problem.getRelativeVelocityComponents(normalRelativeVelocity, tangentialRelativeVelocity);
    analyticTangentialRelativeVelocity = problem.initialTangentialRelativeVelocity *
                                         cos(sqrt(problem.species->getSlidingStiffness() / problem.m12 * 3.5) *
                                             problem.getCollisionTime(2.0 * problem.m12));
    logger(INFO, "tangentialRelativeVelocity: analytic=%, simulation=%\n",
           problem.initialTangentialRelativeVelocity, tangentialRelativeVelocity, Flusher::NO_FLUSH);
    logger(INFO, "normalRelativeVelocity: analytic=%, simulation=%",
           -problem.en * problem.initialNormalRelativeVelocity, normalRelativeVelocity);
}

int main()
{
    
    logger(INFO, "\nnote: analytic values are only asymptotically correct as delta->0");
    particleParticleTest();
    wallParticleTest();
    return 0;
}

