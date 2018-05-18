/* Hertzian2DUnitTest */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Species/HertzianViscoelasticFrictionSpecies.h"

class Hertzian2DUnitTest : public Mercury2D
{
    public:

        Hertzian2DUnitTest()
        {

            setName("Hertzian2DSelfTest");
            setDomain({-1,-1,-1},{1,1,1});
            setTimeMax(4);
            setTimeStep(5e-5);
            setSaveCount(200);

            auto spec = new HertzianViscoelasticFrictionSpecies();
            spec->setDensity(1);
            spec->setElasticModulusAndRestitutionCoefficient(6e-3, 0.1);
            spec = speciesHandler.copyAndAddObject(spec);

            /* Collision between fixed and movable */
            auto pf = new BaseParticle;
            pf->setSpecies(spec);
            pf->setRadius(5e-2);
            pf->setPosition(Vec3D(0,0,0));
            pf->fixParticle();
            particleHandler.copyAndAddObject(pf);

            auto pm = new BaseParticle;
            pm->setSpecies(spec);
            pm->setRadius(5e-2);
            pm->setPosition(Vec3D(1,0,0));
            pm->setVelocity(Vec3D(-1,0,0));
            particleHandler.copyAndAddObject(pm);

            /* Collision between two movables */
            pm->setPosition(Vec3D(0,1,0));
            pm->setVelocity(Vec3D(1,0,0));
            particleHandler.copyAndAddObject(pm);
            
            pm->setPosition(Vec3D(1,1,0));
            pm->setVelocity(Vec3D(0,0,0));
            particleHandler.copyAndAddObject(pm);


        }

        ~Hertzian2DUnitTest() override {
        }
};

int main() 
{
    auto problem = new Hertzian2DUnitTest;
    problem->solve();
    return 0;
}

