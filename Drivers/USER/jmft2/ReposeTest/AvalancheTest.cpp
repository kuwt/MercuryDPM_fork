/* AvalancheTest - Slowly release particles above a circular disc to form a
 * cone of particles. The slope of the cone will go between the dynamic and
 * static angles of repose as avalanches occur. 
 * */

#include "AvalancheTest.h"

int main (const int argc, const char** argv)
{
    char label[] = "larger";
    AT_Parameters pars;
    pars.rho = 1;
    pars.g = 1;
    pars.radiusDisc = 1;
    pars.timeMax = 12000;

    pars.timeStep = 5e-4;
    pars.saveEvery = 200;
    pars.collisionTime = 1e-2;
    pars.restitutionCoefficient = 0.2;
    pars.frictionCoefficientSlide = tan(34 * M_PI / 180);
    pars.frictionCoefficientRoll = tan(10 * M_PI / 180);
    pars.particleRadius = 0.04;
    pars.dispersity_min = 0.9;
    pars.dispersity_max = 1.1;

    pars.radiusRelease = 0.01;
    pars.flowRate = 8e-4;
    pars.veleps = 1e-5;

    AvalancheTest problem(pars);
    problem.generator.randomise();

    char problemname[MAX_STRLEN];
    if (strlen(label) > 0)
        snprintf(problemname, MAX_STRLEN, "AvalancheTest-%s", label);
    else 
        strcpy(problemname, "AvalancheTest");
    char haffn[MAX_STRLEN];
    snprintf(haffn, MAX_STRLEN, "%s.haf", problemname);
    problem.setName(problemname);
    problem.heapAngleFile = fopen(haffn, "w");

    problem.solve();

    fprintf(stdout, "restitution coefficient = %f\n", pars.restitutionCoefficient);
    fprintf(stdout, "betaslide = %f deg\n", atan(pars.frictionCoefficientSlide) * 180 / M_PI );
    fprintf(stdout, "betaroll = %f deg\n", atan(pars.frictionCoefficientRoll) * 180 / M_PI );
    fprintf(stdout, "particle radius: [%f:%f]\n", 
            pars.particleRadius*pars.dispersity_min, 
            pars.particleRadius*pars.dispersity_max);

    return 0;
}

