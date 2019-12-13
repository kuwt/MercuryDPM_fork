//
// Created by paolo on 4-10-19.
//
//Ciao!

#ifndef MERCURY_DOUBLEPOROSITYSTRESSSTRAIN_H
#define MERCURY_DOUBLEPOROSITYSTRESSSTRAIN_H

#endif //MERCURY_DOUBLEPOROSITYSTRESSSTRAIN_H


#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
//#include "Boundaries/BaseClusterInsertionBoundary.h"
#include <Boundaries/StressStrainControlBoundary.h>
#include "BaseCluster.h"

class DoublePorosityStressStrain : public Mercury3D{

public:

    DoublePorosityStressStrain();

    DoublePorosityStressStrain(std::vector<Mdouble> radC,std::vector<Vec3D> posC, Mdouble rLSIntra, Mdouble ExpuSMintra, Mdouble rExpCSIntra, Mdouble rLSInter, Mdouble rTS, Mdouble sFC, Mdouble bS, std::string n);

    ~DoublePorosityStressStrain() override;


private:


    void setupInitialConditions() override;

    void actionsAfterTimeStep() override;

    void printTime() const override;

    void setSpecies();

    void makeWalls();

    void makeIsoCFile();

    void dampVelocities();

    //! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis();

    //! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps();

    //! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity();

    //! computes stress
    void computeStress();

    void writeIsoCFile();

    void writeToIsoCFileTarget();

    //! Geometry
    std::vector<Mdouble> radiiCluster;
    std::vector<Vec3D> positionsCluster;
    Mdouble radiusParticle;
    Mdouble boxSize;
    StressStrainControlBoundary w0;
    Mdouble isoStrainDot;

    //! Species and cluister
    LinearPlasticViscoelasticFrictionSpecies* species;
    std::vector<LinearPlasticViscoelasticFrictionSpecies*> speciesVector;
    int nClusters;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;
    Mdouble densityParticle;
    Mdouble slidingFrictionCoefficient;
    Mdouble loadingStiffnessIntra;
    Mdouble unLoadingStiffnessMaxIntra;
    Mdouble loadingStiffnessInter;
    Mdouble cohesionStiffness;
    Mdouble relativeTangentialStiffness;

    //BaseClusterInsertionBoundary* cluster1;
    int nParticles;
    Mdouble sizeDispersityParticleIntra;

    //! Name and output
    std::ostringstream name;
    std::ofstream isoCFile;
    int fileOutputCount;
    Mdouble fileOutputTime;

    //! Data of interest
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlapIntra_;
    Mdouble meanRelativeOverlapIntra_;
    Mdouble maxRelativeOverlapIntra_;
    Mdouble minRelativeOverlapInter_;
    Mdouble meanRelativeOverlapInter_;
    Mdouble maxRelativeOverlapInter_;
    Mdouble meanCoordinationNumber_;
    Mdouble volumeBox;
    Mdouble voidRatio;
    Mdouble e;
    Mdouble em;
    Mdouble eM;
    Mdouble totalParticleVolume;
    Mdouble inertia;
    Mdouble NC;

    //! Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble stressZZIntra;
    Mdouble totalStress;
    Mdouble uniStress;
    Mdouble uniStressIntra;


    //! Computation
    Mdouble t0;
    Mdouble compressionTime;
    Mdouble relaxationTime;
    int stage;
    Mdouble dampingCoefficient;
    bool isStrainRateControlled;
    Mdouble strainRatePerIteration;
    int counterCompression;
    Mdouble relaxPoint1;
    Mdouble relaxPoint2;
    int counterPositiveForce;
    int counterTotalForces;

    std::ofstream resultsFile;
    std::vector<Mdouble> results;

    std::vector<Mdouble> eValues = {
            0.96,
            0.95919,
            0.95837,
            0.95755,
            0.95672,
            0.95588,
            0.95504,
            0.95419,
            0.95334,
            0.95248,
            0.95161,
            0.95074,
            0.94986,
            0.94898,
            0.94809,
            0.94719,
            0.94629,
            0.94538,
            0.94446,
            0.94354,
            0.94261,
            0.94168,
            0.94074,
            0.93979,
            0.93884,
            0.93788,
            0.93692,
            0.93595,
            0.93497,
            0.93399,
            0.933,
            0.93201,
            0.93101,
            0.93,
            0.92899,
            0.92797,
            0.92694,
            0.92591,
            0.92487,
            0.92383,
            0.92278,
            0.92172,
            0.92066,
            0.91959,
            0.91852,
            0.91744,
            0.91635,
            0.91526,
            0.91416,
            0.91305,
            0.91194,
            0.91083,
            0.9097,
            0.90857,
            0.90744,
            0.9063,
            0.90515,
            0.90399,
            0.90283,
            0.90167,
            0.9005,
            0.89932,
            0.89813,
            0.89694,
            0.89575,
            0.89454,
            0.89333,
            0.89212,
            0.8909,
            0.88967,
            0.88844,
            0.8872,
            0.88595,
            0.8847,
            0.88344,
            0.88218,
            0.88091,
            0.87963,
            0.87835,
            0.87706,
            0.87576,
            0.87446,
            0.87316,
            0.87184,
            0.87052,
            0.8692,
            0.86787,
            0.86653,
            0.86519,
            0.86384,
            0.86248,
            0.86112,
            0.85975,
            0.85837,
            0.85699,
            0.85561,
            0.85421,
            0.85282,
            0.85141,
            0.85,
    };



};
