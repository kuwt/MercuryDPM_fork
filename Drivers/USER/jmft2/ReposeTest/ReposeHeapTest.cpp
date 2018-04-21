/* ReposeHeapTest - A simple test to find the static and dynamic angles of
 * repose of a heap of particles by measuring the length and height of a
 * *single* release of particles. 
 *
 * A single release and collapse will not give the same results as a slow,
 * continuous flow of particles (a la Matt)... see AvalancheTest. 
 * */

#include "ReposeHeapTest.h"

int main (const int argc, const char** argv)
{
    /*
    if (argc <= 1)
    {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        exit(-1);
    }
    */

    RHT_Parameters pars;
    pars.timeMax = 50;
    pars.timeStep = 1e-3;
    pars.saveEvery = 20;
    pars.rho = 1;
    pars.collisionTime = 1e-2;
    pars.restitutionCoefficient = 0.2;
    pars.frictionCoefficientSlide = tan(34 * M_PI / 180);
    pars.frictionCoefficientRoll = tan(12 * M_PI / 180);
    pars.particleRadius = 0.05;
    pars.dispersity_min = 0.9;
    pars.dispersity_max = 1.1;

    pars.boundbox = 1.2;
    pars.veleps = 1e-5;

    pars.g = 1;
    pars.thetaMax = 45 * M_PI / 180;
    pars.thetaDot = 1 * M_PI / 180;

    ReposeHeapTest problem(pars);
    problem.setName("ReposeHeapTest-heapReposeAngle");

    // problem.heapAngleFile = fopen("heapangle.txt", "w");

    problem.solve();

    fprintf(stdout, "restitution coefficient = %f\n", pars.restitutionCoefficient);
    fprintf(stdout, "betaslide  = %f deg\n", atan(pars.frictionCoefficientSlide) * 180 / M_PI );
    fprintf(stdout, "betaroll  = %f deg\n", atan(pars.frictionCoefficientRoll) * 180 / M_PI );
    fprintf(stdout, "particle radius: [%f:%f]\n", 
            pars.particleRadius*pars.dispersity_min, 
            pars.particleRadius*pars.dispersity_max);
    fprintf(stdout, "finished collapse at time %f (counter %d)\n", 
            problem.finishedCollapseAt, problem.finishedCollapseAt_ind);
    fprintf(stdout, "heap repose angle is %f deg\n", 
            problem.heapReposeAngle * 180 / M_PI);
    fprintf(stdout, "started rolling at time %f (counter %d)\n", 
            problem.startedRollingAt, problem.startedRollingAt_ind);
    fprintf(stdout, "static repose angle is %f deg\n",
            problem.staticReposeAngle * 180 / M_PI);

    return 0;
}

