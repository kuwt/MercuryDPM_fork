#include "DPMBase.h"
#include "Walls/ArcWall.h"
#include "Species/LinearViscoelasticSpecies.h"


class ArcWallUnitTest : public DPMBase
{
public:
    void setupInitialConditions() override {
        setTimeMax(0);
        setTimeStep(1e-5);
        setMin(Vec3D(-2, -2, -2));
        setMax(Vec3D(2, 2, 2));
        setGravity(Vec3D(0, 0, 0));
    
        const Mdouble pipeRadius = 1;
        const Mdouble particleRadius = 0.1;
        
        auto species = speciesHandler.copyAndAddObject(new LinearViscoelasticSpecies);
        species->setHandler(&speciesHandler);
        species->setDensity(1);
        const Mdouble particleMass = species->getMassFromRadius(particleRadius);
        species->setStiffnessAndRestitutionCoefficient(1.0 * particleRadius, 0, particleMass);
    
        auto arcWall = new ArcWall();
        arcWall->setSpecies(species);
        arcWall->set(Vec3D(0, 0, 1),
                     Vec3D(0, 0, 0),
                     pipeRadius,
                     Vec3D(0, -1, 0),
                     120 * constants::degree);
        arcWall = wallHandler.copyAndAddObject(arcWall);
    
        auto p = new SphericalParticle;
        p->setSpecies(species);
        p->setRadius(particleRadius);
        
        Vec3D force;
        
        // Particle at the bottom of the pipe
        p->setPosition(Vec3D(0, -pipeRadius + 0.9*particleRadius, 0));
        interactionHandler.clear();
        p->setForce(Vec3D());
        computeForcesDueToWalls(p, arcWall);
        force = p->getForce();
        logger.assert_always((fabs(force.X) < 1e-9) && (force.Y > 0),
                             "The force on the particle at % should have no x-component and positive y-component, but it was %",
                             p->getPosition(),
                             force);
        
        // Particle at a 45 degree angle (testing direction)
        p->setPosition(Vec3D(-pipeRadius + 0.9*particleRadius, -pipeRadius + 0.9*particleRadius, 0));
        interactionHandler.clear();
        p->setForce(Vec3D());
        computeForcesDueToWalls(p, arcWall);
        force = p->getForce();
        logger.assert_always((force.X > 0) && (force.X == force.Y),
                             "The force on the particle at % should have equal and positive x-and y-components, but it was %",
                             p->getPosition(),
                             force);
        
        // Particle at a 90 degree angle (marginal case)
        p->setPosition(Vec3D(-pipeRadius + 0.9*particleRadius, 0, 0));
        interactionHandler.clear();
        p->setForce(Vec3D());
        computeForcesDueToWalls(p, arcWall);
        force = p->getForce();
        logger.assert_always((force.X > 0) && (fabs(force.Y) < 1e-9),
                             "The force on the particle at % should have positive x-component and no y-component, but it was %",
                             p->getPosition(),
                             force);
        
        // Cases where no interaction is expected
        
        // Particle not touching the pipe
        p->setPosition(Vec3D(0, -pipeRadius + 1.1*particleRadius, 0));
        interactionHandler.clear();
        p->setForce(Vec3D());
        computeForcesDueToWalls(p, arcWall);
        force = p->getForce();
        logger.assert_always(force.isZero(),
                             "No force was expected for a particle at % not touching the pipe, but the force was %",
                             p->getPosition(),
                             force);
        
        // Particle on the other side from the pipe
        p->setPosition(Vec3D(0, pipeRadius, 0));
        interactionHandler.clear();
        p->setForce(Vec3D());
        computeForcesDueToWalls(p, arcWall);
        force = p->getForce();
        logger.assert_always(force.isZero(),
                             "No force was expected for a particle at % on the far side of the axis from the pipe, but the force was %",
                             p->getPosition(),
                             force);
    }
};

int main(int argc, char *argv[])
{
    ArcWallUnitTest dpm;
    dpm.setName("ArcWallUnitTest");
    dpm.solve(argc, argv);
    return 0;
}
