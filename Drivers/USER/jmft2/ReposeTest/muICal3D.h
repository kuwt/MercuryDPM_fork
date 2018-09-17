#include <cmath>

typedef struct {
    double timeStep, timeMax;
    int saveEvery;
    double length, width, height;
    double theta;

    double rho, g;
    double collisionTime, restitutionCoefficient; 
    double particleRadius, dispersity;
    double betaslide, betaroll, betators;
    double baseRadius, baseDispersity, baseConc;
    double base_betaslide, base_betaroll, base_betators;
} muICal3D_pars_t;

void muICal3D_pars_read(muICal3D_pars_t* pars, const char* parsfn);
void muICal3D_pars_write(FILE* f, muICal3D_pars_t pars);
void muICal3D_pars_preset(muICal3D_pars_t* pars);

void muICal3D_pars_read(muICal3D_pars_t* pars, const char* parsfn) 
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

    while ((read = getline(&str, &len, parsfile)) != -1)
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
        else if (!strcmp(par_type, "timeMax"))
            pars->timeMax = atof(par_value);
        else if (!strcmp(par_type, "saveEvery"))
            pars->saveEvery = atoi(par_value);
        else if (!strcmp(par_type, "length"))
            pars->length = atof(par_value);
        else if (!strcmp(par_type, "width"))
            pars->width = atof(par_value);
        else if (!strcmp(par_type, "height"))
            pars->height = atof(par_value);
        else if (!strcmp(par_type, "theta"))
            pars->theta = atof(par_value) * M_PI / 180.;

        else if (!strcmp(par_type, "rho"))
            pars->rho = atof(par_value);
        else if (!strcmp(par_type, "g"))
            pars->g = atof(par_value);
        else if (!strcmp(par_type, "collisionTime"))
            pars->collisionTime = atof(par_value);
        else if (!strcmp(par_type, "restitutionCoefficient"))
            pars->restitutionCoefficient = atof(par_value);

        else if (!strcmp(par_type, "particleRadius"))
            pars->particleRadius = atof(par_value);
        else if (!strcmp(par_type, "dispersity"))
            pars->dispersity = atof(par_value);
        else if (!strcmp(par_type, "betaslide"))
            pars->betaslide = atof(par_value) * M_PI / 180.;
        else if (!strcmp(par_type, "betaroll"))
            pars->betaroll = atof(par_value) * M_PI / 180.;
        else if (!strcmp(par_type, "betators"))
            pars->betators = atof(par_value) * M_PI / 180.;

        else if (!strcmp(par_type, "baseRadius"))
            pars->baseRadius = atof(par_value);
        else if (!strcmp(par_type, "baseDispersity"))
            pars->baseDispersity = atof(par_value);
        else if (!strcmp(par_type, "baseConc"))
            pars->baseConc = atof(par_value);
        else if (!strcmp(par_type, "base_betaslide"))
            pars->base_betaslide = atof(par_value) * M_PI / 180.;
        else if (!strcmp(par_type, "base_betaroll"))
            pars->base_betaroll = atof(par_value) * M_PI / 180.;
        else if (!strcmp(par_type, "base_betators"))
            pars->base_betators = atof(par_value) * M_PI / 180.;

        else
            fprintf(stderr, "Warning: Ignored parameter %s, value %s\n",
                    par_type, par_value);
	}

    fprintf(stderr, "Read the following parameters:\n");
    muICal3D_pars_write(stderr, *pars);

    fclose(parsfile);
}

void muICal3D_pars_write(FILE* f, muICal3D_pars_t pars)
{
    fprintf(f, "timeStep %g\n", pars.timeStep);
    fprintf(f, "timeMax %g\n", pars.timeMax);
    fprintf(f, "saveEvery %d\n", pars.saveEvery);

    fprintf(f, "length %lf\n", pars.length);
    fprintf(f, "width %lf\n", pars.width);
    fprintf(f, "height %lf\n", pars.height);
    fprintf(f, "theta %lf\n", pars.theta * 180 / M_PI);

    fprintf(f, "rho %f\n", pars.rho);
    fprintf(f, "g %f\n", pars.g);
    fprintf(f, "collisionTime %lf\n", pars.collisionTime);
    fprintf(f, "restitutionCoefficient %lf\n", 
                    pars.restitutionCoefficient);

    fprintf(f, "particleRadius %g\n", pars.particleRadius);
    fprintf(f, "dispersity %lf\n", pars.dispersity);
    fprintf(f, "betaslide %lf\n", pars.betaslide * 180 / M_PI);
    fprintf(f, "betaroll %lf\n",  pars.betaroll  * 180 / M_PI);
    fprintf(f, "betators %lf\n",  pars.betators  * 180 / M_PI);

    fprintf(f, "baseRadius %g\n", pars.baseRadius);
    fprintf(f, "baseDispersity %lf\n", pars.baseDispersity);
    fprintf(f, "baseConc %lf\n", pars.baseConc);
    fprintf(f, "base_betaslide %lf\n", pars.base_betaslide * 180 / M_PI);
    fprintf(f, "base_betaroll %lf\n",  pars.base_betaroll  * 180 / M_PI);
    fprintf(f, "base_betators %lf\n",  pars.base_betators  * 180 / M_PI);

}

// void muICal3D_pars_preset(muICal3D_pars_t *pars, char* label, material mat)
void muICal3D_pars_preset(muICal3D_pars_t *pars)
{
    pars->timeStep = 8e-5;
    pars->timeMax = 500;
    pars->saveEvery = 1250;

    pars->length    = 1;
    pars->width     = 0.2;
    pars->height    = 0.25;
    pars->theta     = 0. * M_PI / 180.;

    pars->rho = 1;
    pars->g = 1;
    pars->collisionTime = 4e-3;
    pars->restitutionCoefficient = 0.66;

    pars->particleRadius =  0.005;
    pars->dispersity     =  0.1;
    pars->betaslide      = 21.8 * M_PI / 180.;
    pars->betaroll       =  0.0 * M_PI / 180.;
    pars->betators       =  0.0 * M_PI / 180.;
    pars->baseRadius     =  0.005;
    pars->baseDispersity =  0.1;
    pars->baseConc       =  1.0;
    pars->base_betaslide = 21.8 * M_PI / 180.;
    pars->base_betaroll  =  0.0 * M_PI / 180.;
    pars->base_betators  =  0.0 * M_PI / 180.;

}
