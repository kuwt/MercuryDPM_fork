//
// Created by paolo on 6-10-19.
//
//Ciao!

#ifndef MERCURY_BASEDOUBLEPOROSITYSTRESSSTRAIND_H
#define MERCURY_BASEDOUBLEPOROSITYSTRESSSTRAIND_H

#endif //MERCURY_BASEDOUBLEPOROSITYSTRESSSTRAIND_H

#include <iostream>
#include <Mercury3D.h>
#include <Particles/SphericalParticle.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/HertzianViscoelasticFrictionSpecies.h>

class BaseDoublePorosityStressStrainD : public Mercury3D {

public:

    BaseDoublePorosityStressStrainD();

    BaseDoublePorosityStressStrainD(Mdouble pD, Mdouble lS, Mdouble rTS, Mdouble sFC, Mdouble bS, Mdouble sRPTS, std::string n);

    BaseDoublePorosityStressStrainD(Mdouble pD, Mdouble tC, Mdouble r, Mdouble rTS, Mdouble sFC, Mdouble bS, int nP, std::string n);

    BaseDoublePorosityStressStrainD(Mdouble pD, Mdouble lS, Mdouble rTS, Mdouble sFC, Mdouble bS, std::string n);

    ~BaseDoublePorosityStressStrainD() override;

    std::vector<Vec3D> clusterPositions;
    std::vector<Mdouble> clusterRadii;


private:

    //! Initializing variables, boxSize, species, radii, walls, particles, timeStep, output, name
    void setupInitialConditions() override;

    //! Makes data analysis and writes to isotropic compression file
    void actionsBeforeTimeStep() override;

    //! Damps velocities
    void actionsAfterTimeStep() override;

    //! Prints time and values of interest
    void printTime() const override;

    void dampVelocities();

    //! Creating stress strain control walls
    void createWalls();

    void insertParticles();

    bool particleInsertionSuccessful(int n);

    void setSpecies();

    //! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis();

    //! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps();

    //! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity();

    //! computes stress
    void computeStress();

    /*
     * -------------------------- VARIABLES ------------------------------------
     */

    // Simulation
    std::string name;
    int stage;
    Mdouble t0;
    Mdouble isocOutputTimeInterval;
    Mdouble boxSize;
    int counterCompression;
    int counterRelaxation;
    Mdouble relaxationTime;
    Mdouble compressionTime;

    //Species
    Mdouble loadingStiffness;
    Mdouble relativeTangentialStiffness;
    Mdouble slidingFrictionCoefficient;
    Mdouble particleDispersity;

    // Domain
    Mdouble isoStrainDot;
    Mdouble strainRatePerTimeStep;
    Mdouble volumeBox;
    Mdouble initialSolidFraction;

    // Walls
    StressStrainControlBoundary w0;
    StressStrainControlBoundary* w1;

    // Particles
    Mdouble radiusParticle;
    Mdouble volumeParticle;
    Mdouble massParticle;
    Mdouble collisionTime;
    int nParticles;
    Mdouble totalParticleVolume;
    Mdouble dampingCoefficient;
    Mdouble restitutionCoefficient = 0.5;

    // Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;

    // Data analysis
    Mdouble  meanForceOnInteraction;
    Mdouble  meanCoordinationNumber;
    Mdouble  maxRelativeOverlap;
    Mdouble  meanRelativeOverlap;
    Mdouble  minRelativeOverlap;
    Mdouble  voidRatio;
    int floaters;
    Mdouble  e;
    Mdouble rMin;
    Mdouble rMax;

    // Output
    //!\brief cluster data file.
    std::vector<Mdouble> eMValues =
            {
                    0.35172,
                    0.35116,
                    0.3506,
                    0.35003,
                    0.34946,
                    0.34888,
                    0.3483,
                    0.34772,
                    0.34713,
                    0.34654,
                    0.34594,
                    0.34534,
                    0.34473,
                    0.34412,
                    0.34351,
                    0.34289,
                    0.34227,
                    0.34164,
                    0.34101,
                    0.34037,
                    0.33973,
                    0.33909,
                    0.33844,
                    0.33779,
                    0.33713,
                    0.33647,
                    0.33581,
                    0.33514,
                    0.33446,
                    0.33379,
                    0.33311,
                    0.33242,
                    0.33173,
                    0.33103,
                    0.33034,
                    0.32963,
                    0.32893,
                    0.32821,
                    0.3275,
                    0.32678,
                    0.32605,
                    0.32533,
                    0.32459,
                    0.32386,
                    0.32312,
                    0.32237,
                    0.32162,
                    0.32087,
                    0.32011,
                    0.31935,
                    0.31858,
                    0.31781,
                    0.31704,
                    0.31626,
                    0.31547,
                    0.31469,
                    0.3139,
                    0.3131,
                    0.3123,
                    0.3115,
                    0.31069,
                    0.30987,
                    0.30906,
                    0.30824,
                    0.30741,
                    0.30658,
                    0.30575,
                    0.30491,
                    0.30407,
                    0.30322,
                    0.30237,
                    0.30151,
                    0.30066,
                    0.29979,
                    0.29892,
                    0.29805,
                    0.29718,
                    0.2963,
                    0.29541,
                    0.29452,
                    0.29363,
                    0.29273,
                    0.29183,
                    0.29093,
                    0.29002,
                    0.2891,
                    0.28818,
                    0.28726,
                    0.28633,
                    0.2854,
                    0.28447,
                    0.28353,
                    0.28259,
                    0.28164,
                    0.28069,
                    0.27973,
                    0.27877,
                    0.2778,
                    0.27684,
                    0.27586
            };

    std::vector<Mdouble> eMHalves =
            {    0.35144,
                 0.35088,
                 0.350315,
                 0.349745,
                 0.34917,
                 0.34859,
                 0.34801,
                 0.347425,
                 0.346835,
                 0.34624,
                 0.34564,
                 0.345035,
                 0.344425,
                 0.343815,
                 0.3432,
                 0.34258,
                 0.341955,
                 0.341325,
                 0.34069,
                 0.34005,
                 0.33941,
                 0.338765,
                 0.338115,
                 0.33746,
                 0.3368,
                 0.33614,
                 0.335475,
                 0.3348,
                 0.334125,
                 0.33345,
                 0.332765,
                 0.332075,
                 0.33138,
                 0.330685,
                 0.329985,
                 0.32928,
                 0.32857,
                 0.327855,
                 0.32714,
                 0.326415,
                 0.32569,
                 0.32496,
                 0.324225,
                 0.32349,
                 0.322745,
                 0.321995,
                 0.321245,
                 0.32049,
                 0.31973,
                 0.318965,
                 0.318195,
                 0.317425,
                 0.31665,
                 0.315865,
                 0.31508,
                 0.314295,
                 0.3135,
                 0.3127,
                 0.3119,
                 0.311095,
                 0.31028,
                 0.309465,
                 0.30865,
                 0.307825,
                 0.306995,
                 0.306165,
                 0.30533,
                 0.30449,
                 0.303645,
                 0.302795,
                 0.30194,
                 0.301085,
                 0.300225,
                 0.299355,
                 0.298485,
                 0.297615,
                 0.29674,
                 0.295855,
                 0.294965,
                 0.294075,
                 0.29318,
                 0.29228,
                 0.29138,
                 0.290475,
                 0.28956,
                 0.28864,
                 0.28772,
                 0.286795,
                 0.285865,
                 0.284935,
                 0.284,
                 0.28306,
                 0.282115,
                 0.281165,
                 0.28021,
                 0.27925,
                 0.278285,
                 0.27732,
                 0.27635,
                 0.1 //This last value is added by me and is here just because eMValues and eMHalves have different lenght
            };

    std::vector<Mdouble> targetValues =
            {
                    500000,
                    515151.5152,
                    530303.0303,
                    545454.5455,
                    560606.0606,
                    575757.5758,
                    590909.0909,
                    606060.6061,
                    621212.1212,
                    636363.6364,
                    651515.1515,
                    666666.6667,
                    681818.1818,
                    696969.697,
                    712121.2121,
                    727272.7273,
                    742424.2424,
                    757575.7576,
                    772727.2727,
                    787878.7879,
                    803030.303,
                    818181.8182,
                    833333.3333,
                    848484.8485,
                    863636.3636,
                    878787.8788,
                    893939.3939,
                    909090.9091,
                    924242.4242,
                    939393.9394,
                    954545.4545,
                    969696.9697,
                    984848.4848,
                    1000000,
                    1015151.5152,
                    1030303.0303,
                    1045454.5455,
                    1060606.0606,
                    1075757.5758,
                    1090909.0909,
                    1106060.6061,
                    1121212.1212,
                    1136363.6364,
                    1151515.1515,
                    1166666.6667,
                    1181818.1818,
                    1196969.697,
                    1212121.2121,
                    1227272.7273,
                    1242424.2424,
                    1257575.7576,
                    1272727.2727,
                    1287878.7879,
                    1303030.303,
                    1318181.8182,
                    1333333.3333,
                    1348484.8485,
                    1363636.3636,
                    1378787.8788,
                    1393939.3939,
                    1409090.9091,
                    1424242.4242,
                    1439393.9394,
                    1454545.4545,
                    1469696.9697,
                    1484848.4848,
                    1500000,
                    1515151.5152,
                    1530303.0303,
                    1545454.5455,
                    1560606.0606,
                    1575757.5758,
                    1590909.0909,
                    1606060.6061,
                    1621212.1212,
                    1636363.6364,
                    1651515.1515,
                    1666666.6667,
                    1681818.1818,
                    1696969.697,
                    1712121.2121,
                    1727272.7273,
                    1742424.2424,
                    1757575.7576,
                    1772727.2727,
                    1787878.7879,
                    1803030.303,
                    1818181.8182,
                    1833333.3333,
                    1848484.8485,
                    1863636.3636,
                    1878787.8788,
                    1893939.3939,
                    1909090.9091,
                    1924242.4242,
                    1939393.9394,
                    1954545.4545,
                    1969696.9697,
                    1984848.4848,
                    2000000
            };

};














