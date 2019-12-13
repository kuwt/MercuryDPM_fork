//
// Created by paolo on 1-10-19.
//

#include <iostream>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"

class prova : public Mercury3D{

    void setupInitialConditions() override{
        setDomain({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
        insertParticles();
        std::cout<<" nParticles: " << particleHandler.getSize() << std::endl;
        std::cout<<" Mass Fraction: " << particleHandler.getVolume() << std::endl;
    }

    void actionsAfterSolve() override {
        std::cout<<" interactions: " << interactionHandler.getSize() << std::endl;
    }

    void insertParticles()
    {
        int nParticlesInserted = 0;

        while (particleInsertionSuccessful(nParticlesInserted))
        {
            nParticlesInserted++;
        }
        logger(DEFAULT, "PARTICLE INSERTION TERMINATED SUCCESSFULLY\n");
    }

    bool particleInsertionSuccessful(int n) {

        int insertionFailCounter = 0;
        Mdouble rad, theta, phi;
        Vec3D particlePosition;
        SphericalParticle p0;

        //! setup of particle properties and initial conditions (besides position)
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(0.02*random(0.0198, 1.98));
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setGroupId(0);

        while (insertionFailCounter < 1000)
        {
            theta = constants::pi * random.getRandomNumber(0, 2.0);
            phi = acos(random.getRandomNumber(-1.0, 1.0));
            rad = p0.getRadius() + cbrt( random.getRandomNumber( 0, 1 ) ) * ( 0.5 * (getXMax()-getXMin()) - 2.01 * p0.getRadius() );

            particlePosition.X = rad * sin(phi) * cos(theta);
            particlePosition.Y = rad * sin(phi) * sin(theta);
            particlePosition.Z = rad * cos(phi);

            p0.setPosition(particlePosition);

            if (checkParticleForInteraction(p0))
            {
                particleHandler.copyAndAddObject(p0);
                return true;
            }

            insertionFailCounter++;
        }

        return false;
    }






};



int main(){
    prova problem;
    problem.random.randomise();

    problem.setTimeMax(0.01);

//! [T1:speciesProp]
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density
    species.setStiffness(258.5);//sets the spring stiffness.
    species.setDissipation(0.0); //sets the dissipation.
    problem.speciesHandler.copyAndAddObject(species);
    problem.setTimeStep(species.getCollisionTime(0.001)/50);
    problem.setName("ELIMI");
    problem.solve();

}