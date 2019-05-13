#include "Mercury3D.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
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
#define MAX_STRLEN 1024

#define BOTTOMWALL
#define ROUGHBASE
// #define CYLDELB

#ifdef CYLDELB
#include "Boundaries/CylindricalDeletionBoundary.h"
#endif

// #define DELB

typedef struct {
    double timeStep, timeMax;
    int saveEvery;
    double rho, collisionTime, restitutionCoefficient;
    double particleRadius, dispersity_min, dispersity_max;
    double radiusDisc, radiusRelease, flowRate, veleps;

    double frictionCoefficientSlide;
    double frictionCoefficientRoll;
    double g;
} AT_Parameters;

void atpars_write(char* str, size_t size, AT_Parameters pars)
{
    snprintf(str, size,
            "timeStep %f timeMax %f saveEvery %d\n"
            "rho %f collisionTime %f restitutionCoefficient %f\n"
            "particleRadius %f dispersity %f to %f\n"
            "radiusDisc %f radiusRelease %f flowRate %f veleps %f\n"
            "beta_s %f beta_r %f\n"
            "g %f\n",
            pars.timeStep, pars.timeMax, pars.saveEvery,
            pars.rho, pars.collisionTime, pars.restitutionCoefficient,
            pars.particleRadius, pars.dispersity_min, pars.dispersity_max,
            pars.radiusDisc, pars.radiusRelease, pars.flowRate, pars.veleps,
            atan(pars.frictionCoefficientSlide)*180/M_PI,
            atan(pars.frictionCoefficientRoll)*180/M_PI,
            pars.g);
}

class AvalancheTest : public Mercury3D
{
    public:
        AvalancheTest(AT_Parameters pars_in)
        {
            fprintf(stderr, "constructor for AvalancheTest\n");
            pars = pars_in;

            setTimeMax(pars.timeMax);
            setTimeStep(pars.timeStep);
            setSaveCount(pars.saveEvery);

            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);

            /* Define species */
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

            setGravity(Vec3D(0,0,-pars.g));

#ifdef BOTTOMWALL
            /* Plane */
            InfiniteWall wall;
            wall.setSpecies(speciesP);
            wall.set(Vec3D(0,0,-1), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(wall);
#endif

#ifdef ROUGHBASE 
            /* Rough base */
            int Nbasal = pow(pars.radiusDisc / pars.particleRadius, 2);
            for (int i = 0; i < Nbasal; i++)
            {
                SphericalParticle pbasal;
                pbasal.setSpecies(speciesP);
                pbasal.setRadius(pars.particleRadius * generator.getRandomNumber
                        (pars.dispersity_min, pars.dispersity_max));
                // uniformly random position inside a disc 
                double r = pars.radiusDisc * sqrt(generator.getRandomNumber(0,1));
                double theta = generator.getRandomNumber(0,2*M_PI);
                pbasal.setPosition(Vec3D(
                            r * cos(theta), r * sin(theta), 0));
                pbasal.setVelocity(Vec3D(0,0,0));
                pbasal.fixParticle();
                particleHandler.copyAndAddObject(pbasal);
            }
#endif

            /* Insertion and deletion boundaries */
            /*
            CubeInsertionBoundary insb;
            BaseParticle insp;
            insp.setSpecies(speciesP);
            insb.set(&insp, 1, Vec3D(-pars.radiusRelease,-pars.radiusRelease,pars.radiusDisc), 
                                Vec3D(+pars.radiusRelease,+pars.radiusRelease,pars.radiusDisc*2),
                                Vec3D(0,0,0), Vec3D(0,0,0), 
                                pars.particleRadius * pars.dispersity_min,
                                pars.particleRadius * pars.dispersity_max);
            boundaryHandler.copyAndAddObject(insb);
            */

#ifdef CYLDELB
            CylindricalDeletionBoundary cyldelb;
            cyldelb.set(pars.radiusDisc);
            boundaryHandler.copyAndAddObject(cyldelb);

            /* Need this extra large particle to make sure the largest particle
             * never gets deleted. I *think* this will stop segfaults from the
             * deletion boundary. */
            SphericalParticle phack;
            phack.setSpecies(speciesP);
            phack.setRadius(pars.particleRadius * pars.dispersity_max * 2);
            phack.setPosition(Vec3D(0,0, -50*pars.particleRadius));
            phack.setVelocity(Vec3D(0,0,0));
            phack.fixParticle();
            particleHandler.copyAndAddObject(phack);
#endif

#ifdef DELB
            DeletionBoundary delb;
            delb.set(Vec3D(+1,0,0), pars.radiusDisc);
            boundaryHandler.copyAndAddObject(delb);
            delb.set(Vec3D(-1,0,0), pars.radiusDisc);
            boundaryHandler.copyAndAddObject(delb);
            delb.set(Vec3D(0,+1,0), pars.radiusDisc);
            boundaryHandler.copyAndAddObject(delb);
            delb.set(Vec3D(0,-1,0), pars.radiusDisc);
            boundaryHandler.copyAndAddObject(delb);
#endif

            reposeAngle = 0;
            nextInsertionTime = (4. * M_PI * pars.rho * pow(pars.particleRadius,3)) / (3. * pars.flowRate);
            fprintf(stderr, "insertion once every %f\n", nextInsertionTime);
            fprintf(stderr, "2a/v = %f\n", 2*pars.particleRadius*pars.dispersity_max / sqrt(pars.g * pars.radiusDisc));
            step = 0;
        }

        ~AvalancheTest() {
            fclose(heapAngleFile); 
        }

        void actionsAfterTimeStep()
        {
            reposeAngle = getHeapReposeAngle();
            if (step++ % pars.saveEvery == 0)
            {
                fprintf(heapAngleFile, "%f %d %f\n",
                        getTime(), particleHandler.getNumberOfObjects(), 
                        reposeAngle * 180 / M_PI);
                fflush(heapAngleFile);
            }

            /* Introduce a new particle, or not, depending on the flow rate. */
            if (getTime() > nextInsertionTime)
            {
                SphericalParticle pnew;
                pnew.setSpecies(speciesP);
                pnew.setRadius(
                        pars.particleRadius * generator.getRandomNumber(
                            pars.dispersity_min, pars.dispersity_max));
                /* Assume that the angle will never exceed 45 degrees... */
                pnew.setPosition(Vec3D(
                            generator.getRandomNumber(-1,1) * pars.radiusRelease,
                            generator.getRandomNumber(-1,1) * pars.radiusRelease,
                            pars.radiusDisc
                            ));
                pnew.setVelocity(Vec3D(0,0,-sqrt(pars.g * pars.radiusDisc)));
                // pnew.setVelocity(Vec3D(0,0,0));
                particleHandler.copyAndAddObject(pnew);
                nextInsertionTime += (4. * M_PI * pars.rho * pow(pars.particleRadius,3)) / (3. * pars.flowRate);

            }

        }

        double getHeapReposeAngle() 
        {
            double n = particleHandler.getNumberOfObjects();
            double* rs = (double*)calloc(n, sizeof(double));
            double* zs = (double*)calloc(n, sizeof(double));
            int i = 0; int actual_n = n;
            for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
            {
                BaseParticle* p = *it;
                Vec3D pos = p->getPosition();
                double rcoord = sqrt(pos.X*pos.X + pos.Y*pos.Y);
                // We only count this particle if it is unfixed, but unmoving, and inside the
                // domain that we care about..
                if (!(p->isFixed()) 
                        && rcoord < pars.radiusDisc 
                        && pos.Z > 0
                        && pos.Z < pars.radiusDisc 
                        && p->getVelocity().getLength() < 0.1 * sqrt(pars.g*pars.radiusDisc) // what should this number be..?
                )
                {
                    rs[i] = rcoord;
                    zs[i] = pos.Z;
                    i++;
                }
                else
                {
                    actual_n--;
                }
            }
            // fprintf(stderr, "actual_n is %d\n", actual_n);
            n = actual_n;
            rs = (double*)realloc(rs, sizeof(double)*n);
            zs = (double*)realloc(zs, sizeof(double)*n);

            double height, spread, theta;
            if (n > 0)
            {
                qsort(rs, n, sizeof(double), qsort_cmp);
                qsort(zs, n, sizeof(double), qsort_cmp);

                double perc = 0.99;
                spread = getPercentile(rs, n, perc);
                // spread = pars.radiusDisc;
                height = getPercentile(zs, n, perc);
                // fprintf(stdout, "length %f height %f\n", length, height);
                theta = atan(height/spread);
            }
            else
                theta = std::numeric_limits<double>::quiet_NaN();
            free(rs); free(zs);
            return theta;
        }

    // private:
        AT_Parameters pars;
        RNG generator;
        LinearViscoelasticFrictionSpecies* speciesP;
        double reposeAngle;
        double nextInsertionTime;
        unsigned long int step;

        FILE * heapAngleFile;
};

