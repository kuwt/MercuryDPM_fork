//
// Created by paolo on 4-10-19.
//
//Ciao!

#ifndef MERCURY_DOUBLEPOROSITYSTRESSSTRAINC_H
#define MERCURY_DOUBLEPOROSITYSTRESSSTRAINC_H

#endif //MERCURY_DOUBLEPOROSITYSTRESSSTRAINC_H


#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
//#include "Boundaries/BaseClusterInsertionBoundary.h"
#include <Boundaries/StressStrainControlBoundary.h>
#include "BaseCluster.h"

class DoublePorosityStressStrainC : public Mercury3D{

public:

    DoublePorosityStressStrainC();

    DoublePorosityStressStrainC(std::vector<Mdouble> radC,std::vector<Vec3D> posC, Mdouble rappCompress, Mdouble rExpCSIntra, Mdouble sFC, Mdouble relKt, Mdouble bSX, Mdouble bSY, Mdouble bSZ, std::string n);

    ~DoublePorosityStressStrainC() override;


private:


    void setupInitialConditions() override;

    void actionsAfterTimeStep() override;

    void printTime() const override;

    void setSpecies();

    void setSpecies2();

    void makeWalls();

    void makeIsoCFile();

    void makeOverlFile();

    void dampVelocities();

    //! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis();

    //! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps();

    //! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity();

    //! computes stress
    void computeStress();

    //! Computes void dimension
    void computeVoidDimensions();

    bool particleInsertionSuccessful(Mdouble r);

    void writeOverlFile();

    void writeIsoCFile();

    void writeToIsoCFileTarget();

    void makeContactModelFile();

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
    Mdouble compressionRatio;

    //BaseClusterInsertionBoundary* cluster1;
    int nParticles;
    Mdouble sizeDispersityParticleIntra;

    //! Name and output
    std::ostringstream name;
    std::ofstream isoCFile;
    std::ofstream overlFile;
    std::ofstream cmFile;
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
    // pore size / radius particle
    Mdouble relativePoreSize;

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
    Mdouble stressE0pre;
    Mdouble stressE0;


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
    int counterTractionForce;
    int counterTotalForces;
    bool primoValore;

    std::ofstream resultsFile;
    std::vector<Mdouble> results;
    std::ofstream stressFile;
    std::vector<Mdouble> stress;

    std::vector<Mdouble> eValues = {
            0.72,
            0.71879,
            0.71758,
            0.71636,
            0.71515,
            0.71394,
            0.71273,
            0.71152,
            0.7103,
            0.70909,
            0.70788,
            0.70667,
            0.70545,
            0.70424,
            0.70303,
            0.70182,
            0.70061,
            0.69939,
            0.69818,
            0.69697,
            0.69576,
            0.69455,
            0.69333,
            0.69212,
            0.69091,
            0.6897,
            0.68848,
            0.68727,
            0.68606,
            0.68485,
            0.68364,
            0.68242,
            0.68121,
            0.68,
            0.67879,
            0.67758,
            0.67636,
            0.67515,
            0.67394,
            0.67273,
            0.67152,
            0.6703,
            0.66909,
            0.66788,
            0.66667,
            0.66545,
            0.66424,
            0.66303,
            0.66182,
            0.66061,
            0.65939,
            0.65818,
            0.65697,
            0.65576,
            0.65455,
            0.65333,
            0.65212,
            0.65091,
            0.6497,
            0.64848,
            0.64727,
            0.64606,
            0.64485,
            0.64364,
            0.64242,
            0.64121,
            0.64,
            0.63879,
            0.63758,
            0.63636,
            0.63515,
            0.63394,
            0.63273,
            0.63152,
            0.6303,
            0.62909,
            0.62788,
            0.62667,
            0.62545,
            0.62424,
            0.62303,
            0.62182,
            0.62061,
            0.61939,
            0.61818,
            0.61697,
            0.61576,
            0.61455,
            0.61333,
            0.61212,
            0.61091,
            0.6097,
            0.60848,
            0.60727,
            0.60606,
            0.60485,
            0.60364,
            0.60242,
            0.60121,
            0.6,
    };



};
