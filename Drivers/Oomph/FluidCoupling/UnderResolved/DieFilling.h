//
// Created by mitchel on 1/23/20.
//

#ifndef MERCURYDPM_DIEFILLING_H
#define MERCURYDPM_DIEFILLING_H

#include "../../../../Kernel/Oomph/FluidCoupling/UnderResolved/UnderResolvedCoupling.h"

class DieFilling : public UnderResolvedCoupling
{
public:
    DieFilling(const double& xMin, const double& xMax,
                 const double& yMin, const double& yMax,
                 const double& zMin, const double& zMax,
                 const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                 oomph::TimeStepper* time_stepper_pt = &oomph::Mesh::Default_TimeStepper) :
            UnderResolvedCoupling(xMin,xMax,yMin,yMax,zMin,zMax,nx,ny,nz,time_stepper_pt)
    {
        //pinBC();
        //setBC();
    };
    
    void pinBC() override;
    void setBC() override;
    
    void setupInitialConditions() override;
    
    void setDieLength(double dieLength) {dieLength_ = dieLength;}
    void setDieWidth(double dieWidth) {dieWidth_ = dieWidth;}
    void setDieDepth(double dieDepth) {dieDepth_ = dieDepth;}
    
    double getDieLength() {return dieLength_;}
    double getDieWidth() {return dieWidth_;}
    double getDieDepth() {return dieDepth_;}
    
private:
    // Die dimensions
    double dieLength_;
    double dieWidth_;
    double dieDepth_;
    
    Mdouble scaleFactorPSizeMCC = 1.0;
    
    Mdouble rhoSteel = 7850;
    Mdouble rSteel = 0.5;
    Mdouble dSteel = 10e-4;
    Mdouble slidingFrictionSteel = 0.2;
    Mdouble rollingFrictionSteel = 0.0;
    Mdouble rhoGlass = 2650;
    Mdouble rGlass = 0.8;
    Mdouble dGlass = 10e-2;
    Mdouble slidingFrictionGlass = 0.4;
    Mdouble rollingFrictionGlass = 0.0;
    Mdouble torsionFrictionGlass = 0.0;
    
    Mdouble rhoMCC = 1582;
    Mdouble rMCC = 0.5;
    Mdouble dMCC = 100e-6*scaleFactorPSizeMCC; // MCC agglomerates are ~200 mu m
    Mdouble slidingFrictionMCC = 0.2;
    Mdouble rollingFrictionMCC = 0.0;
    Mdouble torsionFrictionMCC = 0.0;
    
    //Mdouble mercStiffness = pGlass->getLoadingStiffness();
    Vec3D frictionVal = {slidingFrictionGlass, rollingFrictionGlass, torsionFrictionGlass};
    
    LinearPlasticViscoelasticFrictionSpecies* pSteel;
    LinearPlasticViscoelasticFrictionSpecies* pGlass;
    LinearPlasticViscoelasticFrictionSpecies* pMCC;

};

#endif //MERCURYDPM_DIEFILLING_H
