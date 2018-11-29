/* Hertzian2DUnitTest */
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Species/LinearPlasticViscoelasticSpecies.h"


struct ConstantRestitutionSelfTest : public Mercury3D {
    void actionsBeforeTimeStep() override {
        if (interactionHandler.getSize()!=0) {
            collisionTime = getTime();
            maxOverlap = std::max(maxOverlap,interactionHandler.getLastObject()->getOverlap());
        }
    }
    Mdouble collisionTime = 0;
    Mdouble maxOverlap = 0;
};

int main() 
{
    ConstantRestitutionSelfTest dpm;
    dpm.setName("ConstantRestitutionSelfTest");
    dpm.setSaveCount(2);

    auto species = dpm.speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSpecies());
    species->setDensity(1000);
    species->setConstantRestitution(true);
    Mdouble radius = 2.5e-3;
    const Mdouble stiffness = 24067/species->getMassFromRadius(radius);
    //species->setPlasticParameters(0.2*stiffness, stiffness, 0.873*stiffness, 0.05);
    species->setPlasticParameters(stiffness, stiffness, 0.0, 0.05);
    species->setRestitutionCoefficient(0.45, 1.0);
    dpm.setTimeStep(0.02*species->getCollisionTime(1.0));
    dpm.setTimeMax(60.0*dpm.getTimeStep());

    radius = 1e-4;
    Mdouble relativeVelocity = 1;

    BaseParticle particle;
    particle.setSpecies(species);
    particle.setRadius(radius);
    particle.setPosition(Vec3D(-radius,0,0));
    particle.setVelocity(Vec3D(0.5*relativeVelocity,0,0));
    auto particle1 = dpm.particleHandler.copyAndAddObject(particle);

    particle.setPosition(-particle.getPosition());
    particle.setVelocity(-particle.getVelocity());
    auto particle2 = dpm.particleHandler.copyAndAddObject(particle);

    dpm.setDomain(radius*Vec3D(-2,-1,-1),radius*Vec3D(2,1,1));

    dpm.solve();

    //Values that are now independent of the particle radius
    logger(INFO,"Max. Overlap %",dpm.maxOverlap);
    logger(INFO,"Collision time %",dpm.collisionTime);
    logger(INFO,"Restitution %",particle2->getVelocity().X/relativeVelocity);
    return 0;
}

