//
// Created by mitchel on 1/23/20.
//

#ifndef MERCURYDPM_DIEFILLING_H
#define MERCURYDPM_DIEFILLING_H

#include "../../../../Kernel/Oomph/FluidCoupling/UnderResolved/UnderResolvedCoupling.h"

template <class ELEMENT>
class DieFilling : public UnderResolvedCoupling<ELEMENT>
{
public:
    DieFilling(const double& xMin, const double& xMax,
                 const double& yMin, const double& yMax,
                 const double& zMin, const double& zMax,
                 const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                 oomph::TimeStepper* time_stepper_pt = &oomph::Mesh::Default_TimeStepper) :
            UnderResolvedCoupling<ELEMENT>(xMin,xMax,yMin,yMax,zMin,zMax,nx,ny,nz,time_stepper_pt)
    {
        pinBC();
        setBC();

        setupInitialConditions();
    };
    
    void setupInitialConditions() override;
    void setBC() override;
    void pinBC() override;
    
    void setDieLength(double dieLength) {dieLength_ = dieLength;}
    void setDieWidth(double dieWidth) {dieWidth_ = dieWidth;}
    void setDieDepth(double dieDepth) {dieDepth_ = dieDepth;}
    
    double getDieLength() {return dieLength_;}
    double getDieWidth() {return dieWidth_;}
    double getDieDepth() {return dieDepth_;}
    
private:
    double dieLength_;
    double dieWidth_;
    double dieDepth_;
    
    Mdouble scaleFactorPSizeMCC = 1.0;
    
    Mdouble rhoSteel = 7850;
    Mdouble rSteel = 0.5;
    Mdouble dSteel = 4e-3;
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
    
    Vec3D frictionVal = {slidingFrictionGlass, rollingFrictionGlass, torsionFrictionGlass};
    
    LinearPlasticViscoelasticFrictionSpecies* pSteel;
    LinearPlasticViscoelasticFrictionSpecies* pGlass;
    LinearPlasticViscoelasticFrictionSpecies* pMCC;

};

template <class ELEMENT>
void DieFilling<ELEMENT>::setupInitialConditions()
{
    logger(INFO,"Call to MercuryDieFilling::setupInitialConditions()");

    double tc = 50.* DPMBase::getTimeStep();

    this->speciesHandler.clear();

    LinearPlasticViscoelasticFrictionSpecies steel;
    steel.setHandler(&this->speciesHandler);
    steel.setDensity(rhoSteel);
    Mdouble mSteel = steel.getMassFromRadius(dSteel/2.0);
    steel.setCollisionTimeAndRestitutionCoefficient(tc,rSteel,mSteel);
    steel.setUnloadingStiffnessMax(5.0*steel.getLoadingStiffness());
    steel.setPenetrationDepthMax(0.05*dSteel);
    steel.setSlidingFrictionCoefficient(slidingFrictionSteel);
    steel.setSlidingStiffness(2.0 / 7.0 * steel.getLoadingStiffness());
    steel.setSlidingDissipation(2.0 / 7.0 * steel.getDissipation());
    steel.setRollingFrictionCoefficient(rollingFrictionSteel);
    steel.setRollingStiffness(0.2 * steel.getLoadingStiffness());
    steel.setRollingDissipation(0.2 * steel.getDissipation());
    pSteel = this->speciesHandler.copyAndAddObject(steel);

    LinearPlasticViscoelasticFrictionSpecies glass;
    glass.setHandler(&this->speciesHandler);
    glass.setDensity(rhoGlass);
    Mdouble mGlass = glass.getMassFromRadius(dGlass/2.0);
    glass.setCollisionTimeAndRestitutionCoefficient(tc,rGlass,mGlass);
    glass.setUnloadingStiffnessMax(5.0*glass.getLoadingStiffness());
    glass.setPenetrationDepthMax(0.05*dGlass);
    glass.setSlidingFrictionCoefficient(slidingFrictionGlass);
    glass.setSlidingStiffness(2.0 / 7.0 * glass.getLoadingStiffness());
    glass.setSlidingDissipation(2.0 / 7.0 * glass.getDissipation());
    glass.setRollingFrictionCoefficient(rollingFrictionGlass);
    glass.setRollingStiffness(0.2 * glass.getLoadingStiffness());
    glass.setRollingDissipation(0.2 * glass.getDissipation());
    pGlass = this->speciesHandler.copyAndAddObject(glass);

    LinearPlasticViscoelasticFrictionSpecies MCC;
    MCC.setHandler(&this->speciesHandler);
    MCC.setDensity(rhoMCC);
    Mdouble mMCC = MCC.getMassFromRadius(dMCC/2.0);
    MCC.setCollisionTimeAndRestitutionCoefficient(tc,rMCC,mMCC);
    MCC.setUnloadingStiffnessMax(5.0*MCC.getLoadingStiffness());
    MCC.setPenetrationDepthMax(0.05*dMCC);
    MCC.setSlidingFrictionCoefficient(slidingFrictionMCC);
    MCC.setSlidingStiffness(2.0 / 7.0 * MCC.getLoadingStiffness());
    MCC.setSlidingDissipation(2.0 / 7.0 * MCC.getDissipation());
    MCC.setRollingFrictionCoefficient(rollingFrictionMCC);
    MCC.setRollingStiffness(0.2 * MCC.getLoadingStiffness());
    MCC.setRollingDissipation(0.2 * MCC.getDissipation());
    pMCC = this->speciesHandler.copyAndAddObject(MCC);


    SphericalParticle p0;
    p0.setSpecies(pSteel);
    p0.setRadius(dSteel/2.);

    double percXMin = 0.0;
    double percXMax = 1.0;
    double percYMin = 0.0;
    double percYMax = 1.0;
    double percZMin = 0.65;
    double percZMax = 0.95;

    const unsigned int numToInsert = 1;
    for (unsigned int i=0; i<numToInsert; i++)
    {
        Vec3D getRandomPos;
        Vec3D getRandomVel;
        int failcounter = 0;

        do
        {
            getRandomPos.X = DPMBase::random.getRandomNumber(DPMBase::getXMin() + p0.getRadius() + percXMin*(DPMBase::getXMax()-DPMBase::getXMin()), percXMax*DPMBase::getXMax() - p0.getRadius());
            getRandomPos.Y = DPMBase::random.getRandomNumber(DPMBase::getYMin() + p0.getRadius() + percYMin*(DPMBase::getYMax()-DPMBase::getYMin()), percYMax*DPMBase::getYMax() - p0.getRadius());
            getRandomPos.Z = DPMBase::random.getRandomNumber(DPMBase::getZMin() + p0.getRadius() + percZMin*(DPMBase::getZMax()-DPMBase::getZMin()), percZMax*DPMBase::getZMax() - p0.getRadius());

            getRandomPos = {0.5*DPMBase::getXMax(), 0.5*DPMBase::getYMax(), 0.9*DPMBase::getZMax()};

            p0.setPosition(getRandomPos);
            p0.setVelocity(getRandomVel);
            failcounter++;
            if (failcounter == 1000) {logger(INFO,"Failcounter reached"); break;}

        } while (!MercuryBase::checkParticleForInteraction(p0));

        if (failcounter != 1000)
        {
            this->particleHandler.copyAndAddObject(p0);
        }
        else {break;}
    }

    std::cout << "Placed " << this->particleHandler.getNumberOfObjects() << " out of " << numToInsert <<  " particles" << std::endl;

    InfiniteWall wall;
    wall.setSpecies(pSteel);
    wall.set(Vec3D(-1,0,0),{DPMBase::getXMin(),0.0,0.0});
    this->wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D( 1,0,0),{DPMBase::getXMax(),0.0,0.0});
    this->wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0,-1,0),{0.0,DPMBase::getYMin(),0.0});
    this->wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0, 1,0),{0.0,DPMBase::getYMax(),0.0});
    this->wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0, 0,-1),{0.0,0.0,DPMBase::getZMin()});
    this->wallHandler.copyAndAddObject(wall);
}

template <class ELEMENT>
void DieFilling<ELEMENT>::pinBC()
{
    logger(DEBUG,"Call to pinBC()");

    // Boundaries are numbered:
    // 0 is at the bottom
    // 1 2 3 4 from the front  proceeding anticlockwise
    // 5 is at the top

    // ---------------------------------------------------------------------------------------------- //
    // FOR UNIFORM INFLOW in z-ir
    //
    //   ^   ^   ^  ^   ^
    //   |   |   |  |   |
    //
    for (unsigned iBound = 0; iBound < this->mesh_pt()->nboundary(); iBound++)
    {
        if (iBound == 0)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin u and v and w
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(1);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(2);
            }
        }
        else if (iBound == 5)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin u and v and w
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(1);
            }
        }
        else
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin all velocities
                for (unsigned i = 0; i < 3; i++)
                {
                    this->mesh_pt()->boundary_node_pt(iBound, iNode)->pin(i);
                }
            }
        }
    }
}

template <class ELEMENT>
void DieFilling<ELEMENT>::setBC()
{
    logger(DEBUG,"Call to setBC()");

    // Pin redudant pressure dofs
    //Fix 3-th pressure value in first element to 0.0.
    dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(0))->fix_pressure(1,0.0);

    // Boundaries are numbered:
    // 0 is at the bottom
    // 1 2 3 4 from the front  proceeding anticlockwise
    // 5 is at the top

    // Periodic flow
    //double zVelIn = 0.05+ 0.3*sin(9.*getTime());

    // Uniform inflow
    double zVelIn = 0.;// getInflowVel(getTime());// getInflowVel(getTime());

    for (unsigned iBound = 0; iBound < this->mesh_pt()->nboundary(); iBound++)
    {
        if (iBound == 0)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,zVelIn);
            }
        }
        else if (iBound == 5)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
            }
        }
        else if (iBound == 1 || iBound == 3)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
        else if (iBound == 2)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++) {
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
        else if (iBound == 4)
        {
            unsigned long int nNode = this->mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++) {
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                this->mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
    }
}

#endif //MERCURYDPM_DIEFILLING_H
