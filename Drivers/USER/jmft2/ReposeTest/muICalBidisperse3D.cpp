/* muICalBidisperse3D
 * */

#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/PolydisperseInsertionBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <cassert>

#define MAX_STRLEN 1024
#define DEGREES (M_PI / 180.)
#define SPHERE_VOL (4.0 / 3.0 * constants::pi)




class muICalBidisperse3D : public Mercury3D
{
    public:
        muICalBidisperse3D(std::string parsfile, double thetaInDegrees)
        {
            std::ifstream paramsFile(parsfile);
            std::string name;
            double var;
            while (paramsFile >> name >> var)
            {
                pars[name] = var;
            }

            theta = thetaInDegrees * DEGREES;
            setName(parsfile.erase(parsfile.find_last_of('.')) + "-" + std::to_string(thetaInDegrees));

            std::cout << theta << std::endl;

            g = getParam("g", 1);

            setTimeMax(getParam("timeMax", 7));
            setTimeStep(getParam("timeStep", 1e-4));
            setSaveCount((unsigned int) getParam("saveCount", 1e4));

            setXMin(0);
            setXMax(getParam("length", 1.0));
            setYMin(0);
            setYMax(getParam("width", 1.0));
            setZMin(0);
            setZMax(getParam("height", 1.0));

            dataFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);

            smallRadius = getParam("smallRadius", 0.05);
            largeRadius = getParam("largeRadius", 0.1);
            basalRadius = getParam("basalRadius", 0.1);
            double smallRho = getParam("rho", 1.0);
            double smallSliding = getParam("sliding_small", 0) * DEGREES;
            double smallRolling = getParam("rolling_small", 0) * DEGREES;
            double smallTorsion = 0;  // CHANGEME
            double largeRho = getParam("rho", 1.0);
            double largeSliding = getParam("sliding_large", 0) * DEGREES;
            double largeRolling = getParam("rolling_large", 0) * DEGREES;
            double largeTorsion = 0;
            // double basalRho = getParam("basalRho");
            // double basalSliding = getParam("sliding_basal") * DEGREES;
            // double basalRolling = getParam("rolling_basal") * DEGREES;
            double basalRho = getParam("rho", 1.0);
            double basalSliding = getParam("sliding_large", 0) * DEGREES;
            double basalRolling = getParam("rolling_large", 0) * DEGREES;
            double basalTorsion = 0;

            double collisionTime = getParam("collisionTime", 5e-3);
            double restitutionCoefficient = getParam("restitutionCoefficient", 0.1);

            small = new LinearViscoelasticFrictionSpecies();
            small->setDensity(smallRho);
            small->setCollisionTimeAndRestitutionCoefficient(
                collisionTime,
                restitutionCoefficient,
                SPHERE_VOL * pow(smallRadius, 3) * smallRho
            );
            small->setSlidingFrictionCoefficient(tan(smallSliding));
            small->setSlidingStiffness(2.0 / 7.0 * small->getStiffness());
            small->setSlidingDissipation(2.0 / 7.0 * small->getDissipation());
            small->setRollingFrictionCoefficient(tan(smallRolling));
            small->setRollingStiffness(2.0 / 5.0 * small->getStiffness());
            small->setRollingDissipation(2.0 / 5.0 * small->getDissipation());
            small->setTorsionFrictionCoefficient(tan(smallTorsion));
            small = speciesHandler.copyAndAddObject(small);
            logger(
                INFO,
                "small species maximum collision velocity %",
                small->getMaximumVelocity(
                    smallRadius,
                    SPHERE_VOL * pow(smallRadius, 3) * smallRho
                )
            );
    
            smallPrototype = new SphericalParticle();
            smallPrototype->setSpecies(small);
            smallPrototype->setRadius(smallRadius);

            large = new LinearViscoelasticFrictionSpecies();
            large->setDensity(largeRho);
            large->setCollisionTimeAndRestitutionCoefficient(
                collisionTime,
                restitutionCoefficient,
                SPHERE_VOL * pow(largeRadius, 3) * largeRho
            );
            large->setSlidingFrictionCoefficient(tan(largeSliding));
            large->setSlidingStiffness(2.0 / 7.0 * large->getStiffness());
            large->setSlidingDissipation(2.0 / 7.0 * large->getDissipation());
            large->setRollingFrictionCoefficient(tan(largeRolling));
            large->setRollingStiffness(2.0 / 5.0 * large->getStiffness());
            large->setRollingDissipation(2.0 / 5.0 * large->getDissipation());
            large->setTorsionFrictionCoefficient(tan(largeTorsion));
            large = speciesHandler.copyAndAddObject(large);
            logger(
                INFO,
                "large species maximum collision velocity %",
                large->getMaximumVelocity(
                    largeRadius,
                    SPHERE_VOL * pow(largeRadius, 3) * largeRho
                )
            );
    
            largePrototype = new SphericalParticle();
            largePrototype->setSpecies(large);
            largePrototype->setRadius(largeRadius);
    
            basal = new LinearViscoelasticFrictionSpecies();
            basal->setDensity(basalRho);
            basal->setCollisionTimeAndRestitutionCoefficient(
                collisionTime,
                restitutionCoefficient,
                SPHERE_VOL * pow(basalRadius, 3) * basalRho
            );
            basal->setSlidingFrictionCoefficient(tan(basalSliding));
            basal->setSlidingStiffness(2.0 / 7.0 * basal->getStiffness());
            basal->setSlidingDissipation(2.0 / 7.0 * basal->getDissipation());
            basal->setRollingFrictionCoefficient(tan(basalRolling));
            basal->setRollingStiffness(2.0 / 5.0 * basal->getStiffness());
            basal->setRollingDissipation(2.0 / 5.0 * basal->getDissipation());
            basal->setTorsionFrictionCoefficient(tan(basalTorsion));
            basal = speciesHandler.copyAndAddObject(basal);
            logger(
                INFO,
                "basal species maximum collision velocity %",
                basal->getMaximumVelocity(
                    basalRadius,
                    SPHERE_VOL * pow(basalRadius, 3) * basalRho
                )
            );

            InfiniteWall bottomWall;
            bottomWall.setSpecies(basal);
            bottomWall.set(Vec3D(0, 0, -1), Vec3D(0, 0, 0));
            wallHandler.copyAndAddObject(bottomWall);

            /* Periodic boundaries */
            PeriodicBoundary xbounds, ybounds;
            xbounds.set(Vec3D(1, 0, 0), getXMin(), getXMax());
            ybounds.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(xbounds);
            boundaryHandler.copyAndAddObject(ybounds);

            /* Rough base */
            // TODO
        }

        void setupInitialConditions() override
        {
            char muIfileName[MAX_STRLEN];
            snprintf(muIfileName, MAX_STRLEN, "%s.muI", getName().c_str());
            auto muIfile = fopen(muIfileName, "w");
            setbuf(muIfile, nullptr);
            fprintf(muIfile, "time theta n depth mass xmom ke basfx basfy basfz\n");
            logger(INFO, "Started writing to .muI file\n");

            /* Gravity initially points downwards */
            setGravity(Vec3D(0, 0, -g));
            
            /* Random particles */
            insb = boundaryHandler.copyAndAddObject(new PolydisperseInsertionBoundary());
            insb->setGeometry(0,
                              getMin(),
                              getMax(),
                              Vec3D(0, 0, 0),
                              Vec3D(0, 0, 0));
            insb->addGenerandum(smallPrototype, 0.5, 0);
            insb->addGenerandum(largePrototype, 0.5, 0);
            insb->checkBoundaryBeforeTimeStep(this);
        }

        void actionsOnRestart() override
        {
            insb = dynamic_cast<PolydisperseInsertionBoundary*>(boundaryHandler.getObjectById(2));
        }
        
        Vec3D calculateBasalForce()
        {
            Vec3D basalForce;
            for (auto p : particleHandler)
            {
                if (p->isFixed())
                {
                    basalForce += p->getForce();
                }
            }
            basalForce += wallHandler.getObject(0)->getForce();
            
            return basalForce;
        }
        
        void actionsAfterTimeStep() override
        {
            if (not insb->isActivated())
                return;

            auto volume = particleHandler.getVolume();
            if (volume > getTotalVolume() * 0.60)
            {
                insb->deactivate();
                std::cout << "Deactivating the insertion boundary" << std::endl;

                setGravity(Vec3D(g * cos(theta), 0, -g * sin(theta)));
            }
        }

        void writeOutputFiles() override
        {
            Mercury3D::writeOutputFiles();
            // TODO
        }
    
        void printTime() const override
        {
            DPMBase::printTime();
            std::cout << particleHandler.getNumberOfObjects() << "  " << particleHandler.getVolume() << std::endl;
        }
        
        void actionsAfterSolve() override
        {
            dataFile.setFileType(FileType::MULTIPLE_FILES);
            writeDataFile();
            fStatFile.setFileType(FileType::MULTIPLE_FILES);
            writeFStatFile();
        }

    private:
        std::map<std::string, double> pars;
        RNG generator;
        LinearViscoelasticFrictionSpecies* small;
        LinearViscoelasticFrictionSpecies* large;
        LinearViscoelasticFrictionSpecies* basal;
        SphericalParticle* smallPrototype;
        SphericalParticle* largePrototype;
        PolydisperseInsertionBoundary* insb;
        bool notYetRemovedInsb;

        double g;
        double theta;  // slope angle in radians
        double smallRadius;
        double largeRadius;
        double basalRadius;
    
        double getParam(const std::string& key, double def)
        {
            try {
                return pars.at(key);
            }
            catch (std::out_of_range& e) {
                return def;
            }
        }
    
        double getParam(const std::string& key)
        {
            try {
                return pars.at(key);
            }
            catch (std::out_of_range& e) {
                std::cerr << "Could not get key " << key << std::endl;
                throw e;
            }
        }
    
};


int main(const int argc, char** argv)
{
    assert (argc > 2);
    auto problem = new muICalBidisperse3D(argv[1], std::stod(argv[2]));
    argv[2] = argv[0];
    problem->solve(argc - 2, argv + 2);
    delete problem;
    return 0;
}
