//
// Created by paolo on 1-10-19.
//
//Ciao!

#ifndef MERCURY_SINGLECLUSTERSTRESSSTRAIN_H
#define MERCURY_SINGLECLUSTERSTRESSSTRAIN_H

#endif //MERCURY_SINGLECLUSTERSTRESSSTRAIN_H


#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include "BaseCluster.h"
#include <Boundaries/StressStrainControlBoundary.h>

class SingleClusterStressStrain : public Mercury3D{

public:

    SingleClusterStressStrain();

    SingleClusterStressStrain(int nP, Mdouble rP, Mdouble sDP, Mdouble dP, Mdouble rC, Mdouble pDM,
                             Mdouble sFC, Mdouble lS, Mdouble uSM, Mdouble cS, Mdouble rTS);

    ~SingleClusterStressStrain() override;


private:


    void setupInitialConditions() override;

    void actionsAfterTimeStep() override;

    void actionsAfterSolve() override;

    void printTime() const override;

    //void computeExternalForces(BaseParticle * CI) override;

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

    //! increasing cohesive forces
    void increaseCohesiveForces();

    void writeIsoCFile();

private:

    //! Geometry
    Mdouble radiusParticle;
    Mdouble boxSize;
    StressStrainControlBoundary w0;
    Mdouble isoStrainDot;

    //! Species and cluister
    LinearPlasticViscoelasticFrictionSpecies* species;
    LinearPlasticViscoelasticFrictionSpecies* speciesSimulation;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;
    Mdouble densityParticle;
    Mdouble slidingFrictionCoefficient;
    Mdouble loadingStiffness;
    Mdouble unLoadingStiffnessMax;
    Mdouble cohesionStiffness;
    Mdouble relativeTangentialStiffness;
    Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(radiusParticle*(1-sizeDispersityParticle), 3) / 3;
    Mdouble collisionTimeSmallestMass = sqrt(smallestMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );
    Mdouble rCluster;
    Mdouble rForMassFraction;

    int nParticles;
    Mdouble sizeDispersityParticle;

    //! Name and output
    std::ostringstream name;
    std::ofstream isoCFile;
    int fileOutputCount;
    Mdouble fileOutputTime;

    //! Data of interest
    Mdouble meanForceOnInteraction_;
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlap_;
    Mdouble meanRelativeOverlap_;
    Mdouble maxRelativeOverlap_;
    Mdouble meanCoordinationNumber_;
    int floaters;
    Mdouble volumeBox;
    Mdouble voidRatio;
    Mdouble e;
    Mdouble em;
    Mdouble eM;
    Mdouble totalParticleVolume;
    Mdouble inertia;

    //! Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;


    //! Computation
    Mdouble initialPositionWall;
    bool start;
    Mdouble t0;
    int stage;
    Mdouble compressionTime;
    Mdouble relaxationTime;
    Mdouble dampingCoefficient;
    bool relax;
    bool isStrainRateControlled;
    Mdouble strainRatePerIteration;
    int counterTotalForces;
    int counterPositiveForce;

};
