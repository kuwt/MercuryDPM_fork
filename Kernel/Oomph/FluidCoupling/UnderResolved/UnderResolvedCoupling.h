//
// Created by mitchel on 6/4/19.
//

#ifndef MERCURYDPM_UNDERRESOLVEDCOUPLING_H
#define MERCURYDPM_UNDERRESOLVEDCOUPLING_H

#include "MercuryDieFilling.h"
#include "OomphDieFilling.h"
#include <vector>
#include "Oomph/OomphHelpers.h"

// Equation headers
#include "Elements/AndersonJackson.h"

class UnderResolvedCoupling;
namespace getDataFromElement
{
    UnderResolvedCoupling *ptrToCoupledClass = nullptr;
    double getVoidageOfElement(const double &, const oomph::Vector<double> &);
    double getVoidageOfElement_byEl(const int &);
    double getdVoidagedxOfElement_byEl(const int &, const int &);
    double getdVoidagedtOfElement_byEl(const int &);
    oomph::Vector<double> getBodyForceByCoupling_byEl(const int &);
}


class UnderResolvedCoupling : public OomphDieFilling<oomph::RefineableAJQCrouzeixRaviartElement<3> >, public MercuryDieFilling
{
public:
    
    UnderResolvedCoupling() {};
    
    UnderResolvedCoupling(const double& xMin, const double& xMax,
                    const double& yMin, const double& yMax,
                    const double& zMin, const double& zMax,
                    const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                    oomph::TimeStepper* time_stepper_pt = &oomph::Mesh::Default_TimeStepper) :
            xMin_(xMin), xMax_(xMax), yMin_(yMin), yMax_(yMax), zMin_(zMin), zMax_(zMax), nx_(nx), ny_(ny), nz_(nz),
            adaptEveryNFluidTimesteps_(5), updateCouplingEveryNParticleTimesteps_(1)
    {
        getDataFromElement::ptrToCoupledClass = this;
        mercVoidage::dpmPointer = this;
        oomph::Problem::add_time_stepper_pt(new oomph::BDF<2>);

        setXMin(xMin_); setXMax(xMax_);
        setYMin(yMin_); setYMax(yMax_);
        setZMin(zMin_); setZMax(zMax_);
        
        oomph::Problem::mesh_pt() = new oomph::RefineableSimpleCubicMesh<oomph::RefineableAJQCrouzeixRaviartElement<3> >(nx_,ny_,nz_,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_,time_stepper_pt);
        // Set error estimator
        oomph::Z2ErrorEstimator* error_estimator_pt=new oomph::Z2ErrorEstimator;
        mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
    
        pinBC();
        setBC();
        
        // Pass pointer to Reynolds number to elements
        unsigned long nelem=mesh_pt()->nelement();
        for (unsigned e=0;e<nelem;e++)
        {
            auto currEl = dynamic_cast<oomph::RefineableAJQCrouzeixRaviartElement<3> *>(mesh_pt()->element_pt(e));
    
            currEl->re_pt()= &Re_;
            currEl->re_st_pt() = &ReSt_;
            currEl->re_invfr_pt() = &Re_InvFr_;
    
            // USED FOR PLOTTING AND SEMI-RESOLVED
            //dynamic_cast<RefineableAJQCrouzeixRaviartElement<3> *>(mesh_pt()->element_pt(e))->voidage_fct_pt() = &getDataFromElement::getVoidageOfElement;
    
            currEl->voidage_fct_pt_byEl() = &getDataFromElement::getVoidageOfElement_byEl;
            currEl->dvoidage_dx_fct_pt_byEl() = &getDataFromElement::getdVoidagedxOfElement_byEl;
            currEl->dvoidage_dt_fct_pt_byEl() = &getDataFromElement::getdVoidagedtOfElement_byEl;
            
            currEl->body_force_fct_pt_by_coupling() = &getDataFromElement::getBodyForceByCoupling_byEl;
    
            // NEEDED FOR COUPLING THROUGH ELEMENT NUMBER
            currEl->set_number(e); // Adding global number to element
        }
        std::cout <<"Number of equations: " << oomph::Problem::assign_eqn_numbers() << std::endl;

        oomph::Problem::Max_newton_iterations = 3000;
        oomph::Problem::Newton_solver_tolerance = 1e-10;
        oomph::Problem::Max_residuals = 1e3;
    }

    oomph::RefineableSimpleCubicMesh<oomph::RefineableAJQCrouzeixRaviartElement<3>>* mesh_pt()
    {
        return dynamic_cast<oomph::RefineableSimpleCubicMesh<oomph::RefineableAJQCrouzeixRaviartElement<3>>*>(oomph::Problem::mesh_pt());
    }
    
    // General functions for solving time-dependent problems
    void setTimeStepOomph(const double& dtOomph) {dtOomph_ = dtOomph;}
    double getTimeStepOomph() {return dtOomph_;}
    void setSaveCountOomph(const int& saveCountOomph) {saveCountOomph_ = saveCountOomph;}
    int getSaveCountOomph() {return saveCountOomph_;}
    void syncTimeSteps();
    void settleParticles();
    
    virtual void pinBC();
    virtual void setBC();
    void solveSystem(oomph::DocInfo &doc_info);
    void setAdaptEveryNFluidTimesteps(const int nFluidTimesteps) {adaptEveryNFluidTimesteps_ = nFluidTimesteps;}
    int getAdaptEveryNFluidTimesteps() {return adaptEveryNFluidTimesteps_;}
    void setUpdateCouplingEveryNParticleTimesteps(const int nParticleTimesteps) {updateCouplingEveryNParticleTimesteps_ = nParticleTimesteps;}
    int getUpdateCouplingEveryNParticleTimesteps() {return updateCouplingEveryNParticleTimesteps_;}
    
    // Functions regarding scaling
    void setReynoldsNumber(const double& Re) {Re_ = Re;}
    double getReynoldsNumber() {return Re_;}
    void setReynoldsStrouhalNumber(const double& ReSt) {ReSt_ = ReSt;}
    double getReynoldsStrouhalNumber() {return ReSt_;}
    void setReynolds_InverseFroudeNumber(const double& Re_InvFr) {Re_InvFr_ = Re_InvFr;}
    double getReynolds_InverseFroudeNumber() {return Re_InvFr_;}

    double getU() {return U_;}
    double getL() {return L_;}
    double getRmu() {return Rmu_;}
    double getRrho() {return Rrho_;}
    double getT() {return T_;}
    double getg() {return g_;}

    void setU(double U) {U_ = U;}
    void setL(double L) {L_ = L;}
    void setRmu(double Rmu) {Rmu_ = Rmu;}
    void setRrho(double Rrho) {Rrho_ = Rrho;}
    void setT(double T) {T_ = T;}
    void setg(double g) {g_ = g;}

    void scaleReynolds() {Re_ *= (getRmu())/(getRrho() * getU() * getL());}
    void scaleReynoldsStrouhal() {ReSt_ *= (getU()*getT())/getL();}
    void scaleReynolds_InverseFroude(){Re_InvFr_ *= (getL()*getL()*getg())/(getU()*getU());}

    double getVelocityScaling() {return (U_);}// Scaling when going from dimensionless to dimensional
    double getLengthScaling() {return (L_);}// Scaling when going from dimensionless to dimensional
    double getPressureScaling() {return (fluidDynamicViscosity_ * U_)/(L_);} // Scaling when going from dimensionless to dimensional
    double getTimeScaling() {return (T_);}// Scaling when going from dimensionless to dimensional
    double getGravityScaling() {return (g_);}// Scaling when going from dimensionless to dimensional
    double getBodyForceScaling() {return (U_/fluidDynamicViscosity_)/(L_*L_);}// Scaling when going from dimensionless to dimensional
    double getSourceScaling() {return U_/L_;}// Scaling when going from dimensionless to dimensional
    double getDensityScaling() {return Rrho_;}// Scaling when going from dimensionless to dimensional
    double getDynViscosityScaling() {return Rmu_;}// Scaling when going from dimensionless to dimensional

    // Specific functions for time-dependent boundary conditions
    double getInflowVel(const double &);
    bool myHeaviSide(const double &, const double &);

    // Actions during solving
    void actionsAfterTimeStep() override;
    void actions_before_adapt() override;
    void actions_after_adapt() override;

    // Functions regarding fluid->particle interaction
    void updateInteractionForces();
    void getS(oomph::Vector<double> &s, const unsigned &iPart, const int &pInEl);
    void updatePOnPart(const unsigned &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s);
    Vec3D getPressureGradient(const int &pIsInEl, const oomph::Vector<double> &s);
    void updateFdOnPart(const unsigned &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s);
    void updateShearOnPart(const unsigned &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s);
    void updateSaffmanLift(const unsigned &iPart, const int &pInEl, const oomph::Vector<double> &s);
    void updateMagnusLift(const unsigned &iPart, const int &pInEl, const oomph::Vector<double> &s);
    
    Vec3D getPOnPart(const unsigned int &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s);
    Vec3D getFdOnPart(const unsigned &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s);
    
    //Functions regarding particle->fluid interaction
    void addCouplingToResiudals();

    // Functions keeping track of particle-element-voidage information
    void generateLists();
    void updateLists();
    void setListsLengths();
    void updateListOfPInElem();
    void updateListOfVoidageInElem();
    void updateListOfdVoidagedxInElem();
    void updateListOfdVoidagedtInElem();
    void updateNeighbourList();

    void updateListsForSpecificParticle(const unsigned &iPart);
    void updateListOfPInElem(const unsigned &iPart);
    void updateListOfVoidageInElem(const unsigned &iEl);
    void updateListOfdVoidagedxInElem(const unsigned &iEl);
    void updateListOfdVoidagedxAroundElem(const unsigned &iEl);
    void updateHistoryValuesVoidage(const unsigned &iEl);

    void updateListOfdVoidagedtAfterAdapt(const std::vector< std::vector<double> > &histValuesBeforeAdapt, const std::vector<double> &dVoidagedtBeforeAdapt); // Used for passing history to sons upon adapt()

    // Functions regarging adaptivity
    void stickLeavesOfNeighbourRecursiveInDir(oomph::Tree* mainTree_, const int &iTree, std::vector<oomph::Tree*> &allNeighbouringLeaves, const int &dir_);
    void stickLeavesOfEqualOrSmallerSizedNeighbourRecursiveInDir(oomph::Tree* neighbourTree, std::vector<oomph::Tree*> &allNeighbouringLeaves, const int &dir_);
    void getInternalIndex(oomph::Tree* Tree_, int &internalIndex);

    void getNeighbourIndex(int indexTree_, int dir_, int &indexInternalNeighbour_);
    void getExternalNeighbour(oomph::Tree* Tree_, const int &iTree_, int dir_, int &differenceInLevel_, std::vector<int> &indexPath_, oomph::Tree* &neighbourPtr_);
    int convertIndexFromPath(std::vector<int> &indexPath_, const int &differenceInLevel_, const int &dir_);

    std::vector<oomph::Tree*> getNeighboursOfEl(const int &iEl_, const int &dir_) {return listOfRefinableNeighbours_[iEl_][dir_];}
    
    // Access functions of particle-element-voidage information
    int getPInElemFromList(const int &iP_) { return listOfPInElem_[iP_];}
    double getVoidageInElemFromList(const int &iEl_) { return listOfVoidageInElem_[iEl_];}
    std::vector<long int> getNeighbourElementNumbers(const int &iEl_) {return listOfNeighbours_[iEl_];} // unrefineable version
    double getdVoidagedxInElemFromList(const int &iEl_, const int &dir_) {return listOfdVoidagedxInElem_[iEl_][dir_];}
    double getdVoidagedtInElemFromList(const int &iEl_) {return listOfdVoidagedtInElem_[iEl_];}
    oomph::Vector<double> getBodyForceInElemByCoupling(const int &iEl_);
    
    // General functions
    void setInflowVel(double inflowVel) {inflowVel_ = inflowVel;}
    double getInflowVel() {return inflowVel_;}
    
    void setFluidDensity(double fluidDensity) {fluidDensity_ = fluidDensity;}
    double getFluidDensity() {return fluidDensity_;}
    void setFluidDynamicViscosity(double fluidDynamicViscosity) {fluidDynamicViscosity_ = fluidDynamicViscosity;}
    double getFluidDynamicViscosity() {return fluidDynamicViscosity_;}
    void setFluidKinematicViscosity(double fluidKinematicViscosity) {fluidKinematicViscosity_ = fluidKinematicViscosity;}
    double getFluidKinematicViscosity() {return fluidKinematicViscosity_;}
    
    void setBodyForceFraction(double bodyForceFraction) {bodyForceFraction_ = bodyForceFraction;}
    double getBodyForceFraction() {return bodyForceFraction_;}
    void resetBodyForceFraction() {bodyForceFraction_ = 0.0;}
    
    void outputNeighboursOfElements();
    
    // Flags
    void setSettling(bool flag) {settlingOn_ = flag;}
    void setAdaptOn(bool flag) {adaptOn_ = flag;}
    
    // Member variables
    std::vector<unsigned int> listOfPInElem_;
    std::vector<double> listOfVoidageInElem_;
    std::vector< std::vector<double> > listOfdVoidagedxInElem_;
    std::vector<double> listOfdVoidagedtInElem_;
    std::vector< std::vector<double> > historyValuesVoidage_;
    std::vector< std::vector<long int> > listOfNeighbours_;
    
    std::vector< std::vector< std::vector<oomph::Tree*> > > listOfRefinableNeighbours_;
    
    oomph::Vector<oomph::Tree*> prevTrees_;
    oomph::Vector<oomph::Tree*> prevLeaves_;
    
    
private:
    double inflowVel_ = 0.0;
    double bodyForceFraction_ = 1.0;
    
    double Re_ = 0.0;
    double ReSt_ = 0.0;
    double Re_InvFr_ = 0.0;
    
    double U_ = 1.0;
    double L_ = 1.0;
    double Rmu_ = 1.0;
    double Rrho_ = 1.0;
    double T_ = 1.0;
    double g_ = 1.0;
    
    double dtOomph_ = 0;
    int saveCountOomph_ = 0;
    double xMin_;
    double xMax_;
    double yMin_;
    double yMax_;
    double zMin_;
    double zMax_;
    unsigned int nx_;
    unsigned int ny_;
    unsigned int nz_;
    
    int adaptEveryNFluidTimesteps_;
    int updateCouplingEveryNParticleTimesteps_;

    bool settlingOn_ = true;
    bool adaptOn_ = false;
    
    double fluidDensity_;
    double fluidKinematicViscosity_; // At 25 degrees C
    double fluidDynamicViscosity_; // At 25 degrees C
};

#endif //MERCURYDPM_UNDERRESOLVEDCOUPLING_H
