#include "Mercury2D.h"

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
// #define ROUGHBASE
#define DELB

typedef struct 
{
    double timeStep, timeMax;
    int saveEvery;
    double rho, g;
    double particleRadius, dispersity;
    double baseRadius, releaseRadius, releaseHeight, releaseRate;
    // releaseRate is given as an area flux
    double collisionTime, restitutionCoefficient, betaslide, betaroll;
} at2d_pars_t;

class AvalancheTest2D : public Mercury2D
{
    public:
        AvalancheTest2D(at2d_pars_t pars_in)
        {
            fprintf(stderr, "constructor for AvalancheTest2D\n");
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
                    constants::pi*pow(pars.particleRadius, 2)*pars.rho
            );

            speciesP->setSlidingFrictionCoefficient(pars.betaslide); 
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(pars.betaroll);
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(0);
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());

            setGravity(Vec3D(0,-pars.g,0));

#ifdef BOTTOMWALL
            /* Plane */
            InfiniteWall wall;
            wall.setSpecies(speciesP);
            wall.set(Vec3D(0,-1,0), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(wall);
#endif

#ifdef DELB
            /* Deletion boundaries */
            DeletionBoundary delb;
            delb.set(Vec3D(1,0,0), pars.baseRadius);
            boundaryHandler.copyAndAddObject(delb);
            delb.set(Vec3D(-1,0,0), pars.baseRadius);
            boundaryHandler.copyAndAddObject(delb);
#endif
        }

        void setupInitialConditions() 
        {
            char haffn[MAX_STRLEN];
            snprintf(haffn, MAX_STRLEN, "%s.haf", getName().c_str());
            heapAngleFile = fopen(haffn, "w");
            setbuf(heapAngleFile, NULL);
            step = 0;
        }

        void actionsAfterTimeStep()
        {
            /* Calculate the heap repose angle and write it to the heapAngleFile
             */
            if (step++ % pars.saveEvery == 0)
            {
                reposeAngle = getHeapReposeAngle();
                fprintf(heapAngleFile, "%f %d %f\n",
                        getTime(), particleHandler.getNumberOfObjects(),
                        reposeAngle * 180/M_PI);
            }

            /* Introduce new particles maybe. */
            double probinsert = pars.releaseRate * pars.timeStep / (M_PI * pars.rho * pow(pars.particleRadius,2));
            if (generator.getRandomNumber(0,1) < probinsert) 
            {
                SphericalParticle pnew;
                pnew.setSpecies(speciesP);
                pnew.setRadius(pars.particleRadius * (1 + pars.dispersity*generator.getRandomNumber(-1,1)));
                pnew.setPosition(Vec3D(
                            pars.releaseRadius * generator.getRandomNumber(-1,1),
                            pars.releaseHeight,
                            0));
                pnew.setVelocity(Vec3D(0,-pars.restitutionCoefficient*sqrt(pars.g * pars.releaseHeight),0));
                if (checkParticleForInteraction(pnew))
                    particleHandler.copyAndAddObject(pnew);
            }
        }

        double getHeapReposeAngle()
        {
            unsigned int n = particleHandler.getNumberOfObjects();
            double* rs = (double*) calloc(n, sizeof(double));
            double* ys = (double*) calloc(n, sizeof(double));
            int i = 0; int actual_n = n;
            for (std::vector<BaseParticle*>::iterator it = particleHandler.begin();
                    it != particleHandler.end(); ++it)
            {
                const BaseParticle* p = *it;
                Vec3D pos = p->getPosition();
                double r = fabs(pos.X);
                if (       
                           !(p->isFixed())
                        && r > pars.releaseRadius 
                        && pos.Y > 0 
                     // && pos.Y < pars.releaseHeight
                        && p->getVelocity().getLength()
                             < pow(pars.restitutionCoefficient,2)*sqrt(pars.g * pars.releaseHeight ) 
                )
                {
                    rs[i] = r;
                    ys[i] = pos.Y;
                    i++;
                }
                else
                    actual_n--;
            }
            n = actual_n;
            rs = (double*)realloc(rs, sizeof(double)*n);
            ys = (double*)realloc(ys, sizeof(double)*n);

            double height;
            double theta;
            if (n > 0)
            {
                qsort(rs, n, sizeof(double), qsort_cmp);
                qsort(ys, n, sizeof(double), qsort_cmp);
                double perc = 0.99;
                height = getPercentile(ys, n, perc);
                theta = atan(height / (pars.baseRadius - pars.releaseRadius));
            }
            else
                theta = 0;

            fprintf(stderr, "t = %f, n = %d, height = %f, theta = %f deg\n", 
                    getTime(), n, height, theta*180/M_PI);
            free(rs); free(ys);
            return theta;
        }

    private:
        at2d_pars_t pars;
        RNG generator;
        LinearViscoelasticFrictionSpecies* speciesP;
        double reposeAngle;
        unsigned long int step;
        FILE * heapAngleFile;
};
