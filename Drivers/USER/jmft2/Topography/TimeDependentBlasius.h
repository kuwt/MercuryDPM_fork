#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    double timeStep, timeMax;
    int saveEvery;
    double theta;
    double length, height, vel;
    double reservoirTemperature;

    double rho, g;
    double collisionTime, restitutionCoefficient;

    double baseRadius, baseDispersity, baseConc;
    double base_betaslide, base_betaroll;

    double particleRadius;
    double dispersity;
    double betaslide, betaroll;
} TimeDependentBlasius_pars_t;

/* Read a parameter file, with filename parsfn, and store the parameters into
 * pars. */
void TimeDependentBlasius_pars_read(TimeDependentBlasius_pars_t* pars, const char* parsfn);
/* Write the parameters pars to the file handler f. */
void TimeDependentBlasius_pars_write(FILE* f, TimeDependentBlasius_pars_t pars);
/* Return a sensible set of default parameters. */
TimeDependentBlasius_pars_t TimeDependentBlasius_pars_default();

void TimeDependentBlasius_pars_read(TimeDependentBlasius_pars_t* pars, const char* parsfn)
{
    FILE * parsfile = fopen(parsfn, "r");
    if (parsfile == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", parsfn);
        exit(-1);
    }

    char* str = NULL;
    size_t len = 0;
    ssize_t read;

    char par_type[255];
    char par_value[255];

    while (read = getline(&str, &len, parsfile) != -1)
    {
        /* Ignore empty lines or lines which begin with #, for they are
         * comments. */
        if (str[0] == '#' || str[0] == '\n')
            continue;

        /* Parse the line. The first word will be a parameter name and the
         * second will be a parameter value. The type of parameter value is a
         * string; depending on which parameter it is, we may have to use atof
         * or atoi to convert it into a number. */
        sscanf(str, "%s %s", par_type, par_value);

        if (!strcmp(par_type, "timeStep"))
            pars->timeStep = atof(par_value);
        if (!strcmp(par_type, "timeMax"))
            pars->timeMax = atof(par_value);
        if (!strcmp(par_type, "saveEvery"))
            pars->saveEvery = atoi(par_value);
        if (!strcmp(par_type, "theta"))
            pars->theta = atof(par_value) * M_PI / 180.;
        if (!strcmp(par_type, "length"))
            pars->length = atof(par_value);
        if (!strcmp(par_type, "height"))
            pars->height = atof(par_value);
        if (!strcmp(par_type, "vel"))
            pars->vel = atof(par_value);
        if (!strcmp(par_type, "reservoirTemperature"))
            pars->reservoirTemperature = atof(par_value);

        if (!strcmp(par_type, "base_betaslide"))
            pars->base_betaslide = atof(par_value) * M_PI / 180.;
        if (!strcmp(par_type, "base_betaroll"))
            pars->base_betaroll = atof(par_value) * M_PI / 180.;
        if (!strcmp(par_type, "baseRadius"))
            pars->baseRadius = atof(par_value);
        if (!strcmp(par_type, "baseDispersity"))
            pars->baseDispersity = atof(par_value);
        if (!strcmp(par_type, "baseConc"))
            pars->baseConc = atof(par_value);

        if (!strcmp(par_type, "particleRadius"))
            pars->particleRadius = atof(par_value);
        if (!strcmp(par_type, "dispersity"))
            pars->dispersity = atof(par_value);

        if (!strcmp(par_type, "betaslide") || !strcmp(par_type, "beta"))
            pars->betaslide = atof(par_value) * M_PI / 180.;
        if (!strcmp(par_type, "betaroll"))
            pars->betaroll = atof(par_value) * M_PI / 180.;

        if (!strcmp(par_type, "collisionTime"))
            pars->collisionTime = atof(par_value);
        if (!strcmp(par_type, "restitutionCoefficient"))
            pars->restitutionCoefficient = atof(par_value);
        
        if (!strcmp(par_type, "rho"))
            pars->rho = atof(par_value);
        if (!strcmp(par_type, "g"))
            pars->g = atof(par_value);
	}

    fprintf(stderr, "Read the following parameters:\n");
    TimeDependentBlasius_pars_write(stderr, *pars);

    fclose(parsfile);
}

void TimeDependentBlasius_pars_write(FILE* f, TimeDependentBlasius_pars_t pars)
{
    fprintf(f, "timeStep               %g\n", pars.timeStep);
    fprintf(f, "timeMax                %g\n", pars.timeMax);
    fprintf(f, "saveEvery              %d\n", pars.saveEvery);
    fprintf(f, "theta                  %.2f\n", pars.theta * 180 / M_PI);
    fprintf(f, "xmin                   %lf\n", pars.xmin);
    fprintf(f, "xmax                   %lf\n", pars.xmax);

    fprintf(f, "length                 %g\n", pars.length);
    fprintf(f, "height                 %g\n", pars.height);
    fprintf(f, "vel                    %g\n", pars.vel);
    fprintf(f, "reservoirTemperature   %g\n", pars.reservoirTemperature);

    fprintf(f, "base_betaslide         %.2f\n", pars.base_betaslide * 180 / M_PI);
    fprintf(f, "base_betaroll          %.2f\n", pars.base_betaroll * 180 / M_PI);
    fprintf(f, "baseRadius             %g\n", pars.baseRadius);
    fprintf(f, "baseDispersity         %g\n", pars.baseDispersity);
    fprintf(f, "baseConc               %g\n", pars.baseConc);

    fprintf(f, "particleRadius         %g\n", pars.particleRadius);
    fprintf(f, "dispersity             %g\n", pars.dispersity);
    fprintf(f, "betaslide              %.2f\n", pars.betaslide * 180 / M_PI);
    fprintf(f, "betaroll               %.2f\n", pars.betaroll * 180 / M_PI);

    fprintf(f, "collisionTime          %g\n", pars.collisionTime);
    fprintf(f, "restitutionCoefficient %.3f\n", 
                    pars.restitutionCoefficient);
    fprintf(f, "rho                    %lf\n", pars.rho);
    fprintf(f, "g                      %lf\n", pars.g);
}

/* 'Sensible' default values. */
TimeDependentBlasius_pars_t TimeDependentBlasius_pars_default()
{
    TimeDependentBlasius_pars_t pars;
    pars.timeStep = 2e-4; // 1/20 of collision time
    pars.timeMax = 50;   // 5 seconds
    pars.saveEvery = 200; // 250 fps (so a total of 250*5 = 1250 frames)
    pars.xmin = -1.0; // 10cm to the left
    pars.xmax = +1.0; // 10cm to the right

    pars.length = 2.0;
    pars.height = 1.0; 
    pars.vel    = 1.0;
    pars.reservoirTemperature = 0.0;

    pars.base_betaslide  = 21.8 * M_PI/180.;
    pars.base_betaroll   =  0.0 * M_PI/180.;
    pars.baseRadius      =  0.01;
    pars.baseDispersity  =  0.10;
    pars.baseConc        =  1.0;

    /* Ballotini */
    pars.particleRadius  =  0.01; // so that delta = d/H = 1/50 
    pars.dispersity      =  0.10;
    pars.betaslide       = 21.8 * M_PI/180.;
    pars.betaroll        =  0.0 * M_PI/180.;

    pars.theta           = 16.0 * M_PI/180; 
    pars.collisionTime   = 4e-3;
    pars.restitutionCoefficient = 0.66;

    pars.rho = 1;
    pars.g = 1;

    return pars;
}

