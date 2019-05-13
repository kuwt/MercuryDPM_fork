#include "Mercury3D.h"
#include "Boundaries/DeletionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "Math/JonnyTools.h"
#include "File.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>

typedef struct {
    double timeStep, timeMax;
    int saveEvery;
    double rho, collisionTime, restitutionCoefficient;
    double particleRadius, dispersity_min, dispersity_max;
    double boundbox, veleps;

    double frictionCoefficientSlide;
    double frictionCoefficientRoll;
    double g, thetaMax, thetaDot;
} RHT_Parameters;

void rhtpars_write(char* str, size_t size, RHT_Parameters pars)
{
    snprintf(str, size,
            "timeStep %f timeMax %f saveEvery %d\n"
            "rho %f collisionTime %f restitutionCoefficient %f\n"
            "particleRadius %f dispersity %f to %f\n"
            "boundbox %f veleps %f\n"
            "beta_s %f beta_r %f\n"
            "g %f thetaMax %f thetaDot %f\n",
            pars.timeStep, pars.timeMax, pars.saveEvery,
            pars.rho, pars.collisionTime, pars.restitutionCoefficient,
            pars.particleRadius, pars.dispersity_min, pars.dispersity_max,
            pars.boundbox, pars.veleps,
            atan(pars.frictionCoefficientSlide)*180/M_PI,
            atan(pars.frictionCoefficientRoll)*180/M_PI,
            pars.g, pars.thetaMax, pars.thetaDot);
}

void rhtpars_fwrite(FILE * file, RHT_Parameters pars)
{
    char str[1024];
    rhtpars_write(str, 1024, pars);
    fputs(str, file);
}

class ReposeHeapTest : public Mercury3D
{
    public:
        ReposeHeapTest(RHT_Parameters pars_in)
        {
            fprintf(stderr, "constructor for ReposeHeapTest\n");
            pars = pars_in;

            setTimeMax(pars.timeMax);
            setTimeStep(pars.timeStep);
            setSaveCount(pars.saveEvery);
            generator.randomise();

            dataFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);

            /* Define species */
            LinearViscoelasticFrictionSpecies* speciesP;
            speciesP = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
            speciesP->setDensity(pars.rho);
            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    pars.collisionTime,
                    pars.restitutionCoefficient, 
                    (4.0/3.0)*constants::pi*pow(pars.particleRadius, 3)*pars.rho
            );

            speciesP->setSlidingFrictionCoefficient(pars.frictionCoefficientSlide); 
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(pars.frictionCoefficientRoll);
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(pars.frictionCoefficientRoll);
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());

            /* Plane */
            setGravity(Vec3D(0,0,-pars.g));
            wall = wallHandler.copyAndAddObject(InfiniteWall());
            wall->setSpecies(speciesP);
            wall->set(Vec3D(0,0,-1), Vec3D(0,0,0));
            

            /* Deletion boundaries */
            /*
            DeletionBoundary delb;
            delb.set(Vec3D(1,0,0), 4*pars.boundbox);
            boundaryHandler.copyAndAddObject(delb);
            delb.set(Vec3D(-1,0,0), 4*pars.boundbox);
            boundaryHandler.copyAndAddObject(delb);
            */

            SphericalParticle p0;
            p0.setSpecies(speciesP);
            for (double x = -.5*pars.boundbox; x < .5*pars.boundbox; x += 2*pars.particleRadius*pars.dispersity_max)
                for (double y = -.5*pars.boundbox; y < .5*pars.boundbox; y += 2*pars.particleRadius*pars.dispersity_max)
                    for (double z = pars.particleRadius; z < pars.boundbox+pars.particleRadius; z += 2*pars.particleRadius*pars.dispersity_max)
                    {
                        p0.setRadius(pars.particleRadius * 
                                generator.getRandomNumber(pars.dispersity_min, pars.dispersity_max));
                        p0.setPosition(Vec3D(x,y,z));
                        p0.setVelocity(Vec3D(
                                generator.getRandomNumber(-pars.veleps,pars.veleps), 
                                generator.getRandomNumber(-pars.veleps,pars.veleps), 
                                generator.getRandomNumber(-pars.veleps,pars.veleps)
                                    ));
                        particleHandler.copyAndAddObject(p0);
                    }

            startingEnergy = fabs( getGravitationalEnergy() ) ;
            fprintf(stderr, "startingEnergy = %f\n", startingEnergy);
            finishedCollapse = false;
            finishedCollapseAt = 0;
            startedRolling = false;
            startedRollingAt = 0;

            theta = 0;
        }

        ~ReposeHeapTest() {}

        void actionsAfterTimeStep()
        {
            if (!finishedCollapse)
            {
                double last_hra = heapReposeAngle;
                heapReposeAngle = getHeapReposeAngle();
                double dhradt = (getTime() == 0) ? 0 : (heapReposeAngle - last_hra)/getTimeStep();

                fprintf(heapAngleFile, "%f %f %f\n", 
                        getTime(), heapReposeAngle, dhradt);

                double ene_kin = getKineticEnergy();
                // double ene_thm = getThermalEnergy();
                // double eneratio = ene_kin / ene_thm - 1; // bulk over therm

                // if (fabs(dhradt) < 1e-4 && (eneratio < 1e-4 || ene_kin < 1e-20) && getTime() > 1)
                if (fabs(dhradt) < 1e-4 && (ene_kin < 1e-5 * getTotalMass()) && getTime() > 1)
                {
                    // outputXBallsData(dataFile.getFstream());
                    finishedCollapse = true;
                    finishedCollapseAt = getTime();
                    finishedCollapseAt_ind = dataFile.getCounter();
                    fprintf(stderr, "finished collapse at time %f\n", finishedCollapseAt);
                    fprintf(stderr, "heap repose angle is %f deg\n", 
                            heapReposeAngle * 180. / M_PI);
                    fclose(heapAngleFile);

                }
                else
                    setTimeMax(getTimeMax() + getTimeStep());
            }

            if (finishedCollapse) 
            {
                theta = fmin(pars.thetaMax, pars.thetaDot * (getTime() - finishedCollapseAt));
            }

            if (finishedCollapse && !startedRolling)
            {
                if (getGravitationalEnergy() < 0) // This means that we have started to slide, 
                                                  // but this is a sufficient condition, 
                                                  // not a necessary condition.
                {
                    startedRolling = true;
                    startedRollingAt = getTime();
                    startedRollingAt_ind = dataFile.getCounter();
                    staticReposeAngle = heapReposeAngle + theta;
                    fprintf(stderr, "started rolling at time %f, theta=%f deg, thetas = %f deg\n",
                            startedRollingAt, theta * 180 / M_PI, staticReposeAngle * 180 / M_PI);
                }
                else
                    setTimeMax(getTimeMax() + getTimeStep());
            }

            if (startedRolling)
            {
                setTimeMax(getTimeMax() - 1);
            }

            wall->set(Vec3D(-pars.g*sin(theta),0,-pars.g*cos(theta)), Vec3D(0,0,0));
        }

        double getHeapReposeAngle() 
        {
            double n = particleHandler.getNumberOfObjects();
            double* xs = (double*)calloc(n, sizeof(double));
            double* zs = (double*)calloc(n, sizeof(double));
            int i = 0;
            for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
            {
                BaseParticle* p = *it;
                xs[i] = p->getPosition().X;
                zs[i] = p->getPosition().Z;
                // fprintf(stdout, "i %d z %f\n", i, zs[i]);
                i++;
            }
            qsort(xs, n, sizeof(double), qsort_cmp);
            qsort(zs, n, sizeof(double), qsort_cmp);

            double perc = 0.99;
            double length = 0.5 * (getPercentile(xs, n, perc) - getPercentile(xs, n, 1-perc));
            double height = getPercentile(zs, n, 1);
            // fprintf(stdout, "length %f height %f\n", length, height);
            free(xs); free(zs);
            return atan(height / length);
        }

    // private:
        RHT_Parameters pars;
        RNG generator;
        InfiniteWall* wall;

        double theta;
        double startingEnergy;
        bool finished;
        bool finishedCollapse;
        double finishedCollapseAt;
        unsigned int finishedCollapseAt_ind;
        bool startedRolling;
        double startedRollingAt;
        unsigned int startedRollingAt_ind;
        double heapReposeAngle;
        double staticReposeAngle;

        FILE * heapAngleFile;
};

