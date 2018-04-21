/* AvalancheTest2D */
#include "AvalancheTest2D.h"

#include <string.h>

/* Define some useful presets */
typedef enum {
    BALLOTINI, SAND, GRIT, FRICTIONLESS    
} material;

void at2d_pars_preset(at2d_pars_t *pars, char* label, material mat);

int main (const int argc, const char** argv)
{
    char label[MAX_STRLEN];
    at2d_pars_t pars;
    if (argc != 2) 
    {
        fprintf(stderr, "Usage: %s [ballotini|sand|grit|frictionless]\n", argv[0]);
        exit(-1);
    }
    else
    {
        if (strcmp(argv[1], "ballotini") == 0)
            at2d_pars_preset(&pars, label, BALLOTINI);
        if (strcmp(argv[1], "sand") == 0)
            at2d_pars_preset(&pars, label, SAND);
        if (strcmp(argv[1], "grit") == 0)
            at2d_pars_preset(&pars, label, GRIT);
        if (strcmp(argv[1], "frictionless") == 0)
            at2d_pars_preset(&pars, label, FRICTIONLESS);
    }

    AvalancheTest2D problem(pars);

    char problemname[MAX_STRLEN];
    if (strlen(label) > 0)
        snprintf(problemname, MAX_STRLEN, "AvalancheTest2D-%s", label);
    else 
        strcpy(problemname, "AvalancheTest2D");
    problem.setName(problemname);

    problem.solve();

    fprintf(stdout, "restitution coefficient = %f\n", pars.restitutionCoefficient);
    fprintf(stdout, "betaslide = %f deg\n", pars.betaslide * 180 / M_PI );
    fprintf(stdout, "betaroll = %f deg\n", pars.betaroll * 180 / M_PI );
    fprintf(stdout, "particle radius: [%f:%f]\n", 
            pars.particleRadius*(1-pars.dispersity), 
            pars.particleRadius*(1+pars.dispersity));

    return 0;
}

void at2d_pars_preset(at2d_pars_t *pars, char* label, material mat)
{
    pars->timeStep = 8e-5;
    pars->timeMax = 12000;
    pars->saveEvery = 2000;

    pars->rho = 1;
    pars->g = 1;

    pars->baseRadius = 1;
    pars->particleRadius = 0.005;
    pars->releaseRadius = 0.05;
    pars->releaseRate = 0.004;
    pars->releaseHeight = 1.5;

    pars->collisionTime = 4e-3;

    switch(mat)
    {
        case BALLOTINI:
            pars->dispersity = 0.1;
            pars->restitutionCoefficient = 0.66;
            pars->betaslide = tan(21.8 * M_PI / 180.);
            pars->betaroll = tan(0.0 * M_PI / 180.);
            strncpy(label, "ballotini", MAX_STRLEN);
            break;

        case SAND:
            pars->dispersity = 0.25;
            pars->restitutionCoefficient = 0.66;
            pars->betaslide = tan(34.0 * M_PI / 180.);
            pars->betaroll = tan(5.0 * M_PI / 180.);
            strncpy(label, "sand", MAX_STRLEN);
            break;

        case GRIT:
            pars->dispersity = 0.2;
            pars->restitutionCoefficient = 0.66;
            pars->betaslide = tan(34.0 * M_PI / 180.);
            pars->betaroll = tan(25.0 * M_PI / 180.);
            strncpy(label, "grit", MAX_STRLEN);
            break;

        case FRICTIONLESS:
            pars->dispersity = 0.1;
            pars->restitutionCoefficient = 0.66;
            pars->betaslide = tan(0.0 * M_PI / 180.);
            pars->betaroll = tan(0.0 * M_PI / 180.);
            strncpy(label, "frictionless", MAX_STRLEN);
            break;

        default:
            break;
    }
}
