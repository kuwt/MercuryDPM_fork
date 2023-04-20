/* ReposeHeapTest-multiple - Arrays of parameters! */

#include "ReposeHeapTest.h"
#include <pthread.h>
const size_t MAX_STRLEN = 1024;

void sendmail(const char* user, const char* subject, const char* msg, const char* mailbuf_fn)
{
    FILE * mailbuf = fopen(mailbuf_fn, "w");
    fprintf(mailbuf, "%s", msg);
    fclose(mailbuf);
    char command[MAX_STRLEN];
    snprintf(command, MAX_STRLEN, "cat %s | mail -s \"%s\" %s",
            mailbuf_fn, subject, user);
    int result = system(command);
}



typedef struct {
    RHT_Parameters pars;
    ReposeHeapTest* instance;
    char problemname[MAX_STRLEN];
    double theta_h;
    double theta_s;
} ptwrap_s;

/*
void* ptwrap_runsims(void* threadargs)
{
    ptwrap_s* args = (ptwrap_s*) threadargs;
    fprintf(stderr, "AAAAA! %p, %f\n", args->instance, args->pars.frictionCoefficientRoll);
    args->instance->solve();

    args->theta_h = args->instance->heapReposeAngle;
    args->theta_s = args->instance->staticReposeAngle;

    char pars_write[MAX_STRLEN]; rhtpars_write(pars_write, MAX_STRLEN, args->pars);
    char msg[MAX_STRLEN];
    snprintf(msg, MAX_STRLEN, 
            "Exultation! %s is finished.\n"
            "%s\n"
            "beta_s %f deg, beta_r %f deg, theta_h %f deg, theta_s %f deg\n", 
            args->problemname, 
            pars_write,
            atan(args->pars.frictionCoefficientSlide)  * 180 / M_PI,
            atan(args->pars.frictionCoefficientRoll)  * 180 / M_PI,
            args->theta_h * 180 / M_PI,
            args->theta_s * 180 / M_PI
            );
    sendmail("jmft2", args->problemname, msg, "/tmp/mailbuf.txt");

    fprintf(stderr, "Exultation! %s is finished.\n", args->problemname);
}
*/


int main (const int argc, const char** argv)
{
    if (argc == 1)
    {
        fprintf(stderr, "Usage: %s seriesconffile\n", argv[0]);
        exit(-1);
    }

    FILE * seriesconffile = fopen(argv[1], "r");
    /* First line gives the number of simulations in this series and the name of the series. */
    size_t N; char seriesname[MAX_STRLEN];
    fscanf(seriesconffile, "%zd %s\n", &N, seriesname);

    /*
    char seriesname[] = "elasticparticles4";
    size_t N = 40;
    */
    RHT_Parameters pars[N];
    ReposeHeapTest** instances = (ReposeHeapTest**) malloc(N * sizeof(ReposeHeapTest*));
    for (int i = 0; i < N; i++) 
    {
        pars[i].timeMax = 50;
        pars[i].rho = 1;
        pars[i].dispersity_min = 0.9;
        pars[i].dispersity_max = 1.1;
        /*
        pars[i].particleRadius = 0.05;
        pars[i].collisionTime = 1e-2;
        pars[i].restitutionCoefficient = 0.6;
        pars[i].frictionCoefficientSlide = tan(34 * M_PI / 180);
        pars[i].frictionCoefficientRoll = tan((16 + 0.4 * i) * M_PI / 180);
        */
        double betaslide_d, betaroll_d;

        char buf[MAX_STRLEN];
        double dispersity;
        char* result = fgets(buf, MAX_STRLEN, seriesconffile);
        sscanf(buf, "%lf %lf %lf %lf %lf %lf",
                &(pars[i].particleRadius),
                &dispersity,
                &(pars[i].collisionTime),
                &(pars[i].restitutionCoefficient),
                &betaslide_d,
                &betaroll_d
        );
        pars[i].dispersity_min = 1 - dispersity;
        pars[i].dispersity_max = 1 + dispersity;

        pars[i].frictionCoefficientSlide = tan(betaslide_d * M_PI / 180);
        pars[i].frictionCoefficientRoll = tan(betaroll_d * M_PI / 180);
        // pars[i].timeStep = 5e-4;
        pars[i].timeStep = pars[i].collisionTime / 40;
        pars[i].saveEvery = 20;

        pars[i].boundbox = 1.2;
        pars[i].veleps = 1e-5;

        pars[i].g = 1;
        pars[i].thetaMax = 45 * M_PI / 180;
        pars[i].thetaDot = 1 * M_PI / 180;

        // rhtpars_fwrite(stdout, pars[i]);
        instances[i] = new ReposeHeapTest(pars[i]);
        instances[i]->generator.randomise();

        char problemname[MAX_STRLEN];
        snprintf(problemname, MAX_STRLEN, "ReposeHeapTest-%s-run%d out of %zu", seriesname, i, N);
        instances[i]->setName(problemname);
        fprintf(stderr, "%d %p %s\n", i, instances[i], problemname);
        char problem_haf[MAX_STRLEN];
        snprintf(problem_haf, MAX_STRLEN, "%s.haf", problemname);
        instances[i]->heapAngleFile = fopen(problem_haf, "w");
        setbuf(instances[i]->heapAngleFile, NULL);
    }
    double theta_h[N];
    double theta_s[N];

    for (int i = 0; i < N; i++)
    {
        instances[i]->solve();
        theta_h[i] = instances[i]->heapReposeAngle;
        theta_s[i] = instances[i]->staticReposeAngle;
        char pars_write[MAX_STRLEN]; rhtpars_write(pars_write, MAX_STRLEN, pars[i]);
        char msg[MAX_STRLEN];
        snprintf(msg, MAX_STRLEN, 
                "Exultation! %s is finished.\n"
                "%s\n"
                "beta_s %f deg, beta_r %f deg, theta_h %f deg, theta_s %f deg\n", 
                instances[i]->getName().c_str(), 
                pars_write,
                atan(pars[i].frictionCoefficientSlide)  * 180 / M_PI,
                atan(pars[i].frictionCoefficientRoll)  * 180 / M_PI,
                theta_h[i] * 180 / M_PI,
                theta_s[i] * 180 / M_PI
               );
        sendmail("jmft2", instances[i]->getName().c_str(), msg, "/tmp/mailbuf.txt");
    }

    /*
    for (int i = 0; i < N; i++)
    {
        ReposeHeapTest instance(pars[i]);
        instance.generator.randomise();

        char problemname[MAX_STRLEN];
        snprintf(problemname, MAX_STRLEN, "ReposeHeapTest-%s-run%d", seriesname, i);
        instance.setName(problemname);
        fprintf(stderr, "%s\n", problemname);

        instance.solve();
        theta_h[i] = instance.heapReposeAngle;
        theta_s[i] = instance.staticReposeAngle;

        char pars_write[MAX_STRLEN]; rhtpars_write(pars_write, MAX_STRLEN, pars[i]);
        char msg[MAX_STRLEN];
        snprintf(msg, MAX_STRLEN, 
                "Exultation! %s is finished.\n"
                "%s\n"
                "beta_s %f deg, beta_r %f deg, theta_h %f deg, theta_s %f deg\n", 
                problemname, 
                pars_write,
                atan(pars[i].frictionCoefficientSlide)  * 180 / M_PI,
                atan(pars[i].frictionCoefficientRoll)  * 180 / M_PI,
                theta_h[i] * 180 / M_PI,
                theta_s[i] * 180 / M_PI
               );
        sendmail("jmft2", problemname, msg, "/tmp/mailbuf.txt");
    }
    */

    /*
    pthread_t threads[N];
    ptwrap_s thread_data[N];
    for (int i = 0; i < N; i++)
    {
        thread_data[i].pars = pars[i];
        thread_data[i].instance = instances[i];
    }

    for (int i = 0; i < N; i++)
    {
        // pthread_create(&threads[i], NULL, ptwrap_runsims, (void*) &thread_data[i]);
    }
    */

    /* Send a summary email after everything is done. */
    char msgsum[MAX_STRLEN * N];
    for (int i = 0; i < N; i++)
    {
        char msgsum_part[MAX_STRLEN];
        /*
        snprintf(msgsum_part, MAX_STRLEN, "%d: beta_s %f deg, beta_r %f deg, theta_h %f deg, theta_s %f deg\n", 
                i, 
                atan(pars[i].frictionCoefficientSlide)  * 180 / M_PI,
                atan(pars[i].frictionCoefficientRoll)  * 180 / M_PI,
                theta_h[i] * 180 / M_PI,
                theta_s[i] * 180 / M_PI
               );
        */
        snprintf(msgsum_part, MAX_STRLEN, "%f %f %f %f %f\n",
                pars[i].restitutionCoefficient,
                atan(pars[i].frictionCoefficientSlide)  * 180 / M_PI,
                atan(pars[i].frictionCoefficientRoll)  * 180 / M_PI,
                theta_h[i] * 180 / M_PI,
                theta_s[i] * 180 / M_PI
        );
        fprintf(stdout, "%s", msgsum_part);
        strncat(msgsum, msgsum_part, MAX_STRLEN*N);
    }
    char subject[MAX_STRLEN];
    snprintf(subject, MAX_STRLEN, "ReposeHeapTest-%s summary", seriesname);
    sendmail("jmft2", subject, msgsum, "/tmp/mailbuf.txt");

    free(instances);
}
