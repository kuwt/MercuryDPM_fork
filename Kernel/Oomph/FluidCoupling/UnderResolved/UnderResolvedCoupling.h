//
// Created by mitchel on 6/4/19.
//

#ifndef MERCURYDPM_UNDERRESOLVEDCOUPLING_H
#define MERCURYDPM_UNDERRESOLVEDCOUPLING_H

#include "MercuryDieFilling.h"
#include "OomphDieFilling.h"
#include <vector>
#include "Oomph/OomphHelpers.h"
#include "Math/Vector.h"

// Equation headers
#include "Elements/AndersonJackson.h"

#include "../../../../oomph-lib/src/meshes/simple_cubic_mesh.h"
#include "UnderResolvedMesh.h"

template<class ELEMENT>
class UnderResolvedCoupling;

namespace getDataFromElement
{
    template<class ELEMENT>
    UnderResolvedCoupling<ELEMENT> *ptrToCoupledClass = nullptr;
    double getVoidageOfElement(const double &, const oomph::Vector<double> &);
    double getVoidageOfElement_byEl(const int &);
    double getdVoidagedxOfElement_byEl(const int &, const int &);
    double getdVoidagedtOfElement_byEl(const int &);
    void getBodyForceByCoupling_byEl(const int &, const oomph::Vector<double> &);
}

template <class ELEMENT>
class UnderResolvedCoupling : public UnderResolvedMesh<ELEMENT>, public MercuryDieFilling
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
        oomph::Problem::add_time_stepper_pt(new oomph::BDF<2>);

        setXMin(xMin_); setXMax(xMax_);
        setYMin(yMin_); setYMax(yMax_);
        setZMin(zMin_); setZMax(zMax_);
        
        oomph::Problem::mesh_pt() = new oomph::RefineableSimpleCubicMesh<ELEMENT >(nx_,ny_,nz_,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_,time_stepper_pt);

        // Set error estimator
        oomph::Z2ErrorEstimator* error_estimator_pt=new oomph::Z2ErrorEstimator;
        mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
    
        pinBC();
        setBC();

        // Pass pointer to Reynolds number to elements
        unsigned long nelem=mesh_pt()->nelement();
        for (unsigned e=0;e<nelem;e++)
        {
            auto currEl = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
    
            currEl->re_pt()= &Re_;
            currEl->re_st_pt() = &ReSt_;
            currEl->re_invfr_pt() = &Re_InvFr_;
    
            // USED FOR PLOTTING AND SEMI-RESOLVED
            //dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(e))->voidage_fct_pt() = &getDataFromElement::getVoidageOfElement;

            // NEEDED FOR COUPLING THROUGH ELEMENT NUMBER
            currEl->set_number(e); // Adding global number to element
        }
        std::cout <<"Number of equations: " << oomph::Problem::assign_eqn_numbers() << std::endl;

        oomph::Problem::Max_newton_iterations = 3000;
        oomph::Problem::Newton_solver_tolerance = 1e-10;
        oomph::Problem::Max_residuals = 1e3;
    }

    oomph::RefineableSimpleCubicMesh<ELEMENT>* mesh_pt()
    {
        return dynamic_cast<oomph::RefineableSimpleCubicMesh<ELEMENT>*>(oomph::Problem::mesh_pt());
    }

    // General functions for solving time-dependent problems
    void setTimeStepOomph(const double& dtOomph) {dtOomph_ = dtOomph;}
    double getTimeStepOomph() {return dtOomph_;}
    void setSaveCountOomph(const int& saveCountOomph) {saveCountOomph_ = saveCountOomph;}
    int getSaveCountOomph() {return saveCountOomph_;}
    void syncTimeSteps();
    void settleParticles();
    
    virtual void pinBC() {}
    virtual void setBC() {}
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

    // Functions to pass info to fluid element
    void passInfoToFluid();

    // Functions keeping track of particle-element-voidage information
    void prescribeVoidage(double voidageValue, double from, double to, int dir);

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
    void getBodyForceInElemByCoupling(const int &iEl_, oomph::Vector<double> &force_);

    // General functions
    virtual void setInflowVel(double inflowVel) {inflowVel_ = inflowVel;}
    virtual double getInflowVel() {return inflowVel_;}
    
    void setFluidDensity(double fluidDensity) {fluidDensity_ = fluidDensity;}
    double getFluidDensity() {return fluidDensity_;}
    void setFluidDynamicViscosity(double fluidDynamicViscosity) {fluidDynamicViscosity_ = fluidDynamicViscosity;}
    double getFluidDynamicViscosity() {return fluidDynamicViscosity_;}
    void setFluidKinematicViscosity(double fluidKinematicViscosity) {fluidKinematicViscosity_ = fluidKinematicViscosity;}
    double getFluidKinematicViscosity() {return fluidKinematicViscosity_;}

    void setIncrementalFluidSolve(bool x) {incrementalFluidSolve_ = x;}
    bool getIncrementalFluidSolve() {return incrementalFluidSolve_;}

    void setBodyForceFraction(double bodyForceFraction) {bodyForceFraction_ = bodyForceFraction;}
    double getBodyForceFraction() {return bodyForceFraction_;}
    void resetBodyForceFraction() {bodyForceFraction_ = 0.0;}
    
    void outputNeighboursOfElements();
    void doc_solution(oomph::DocInfo&) override;
    void doc_voidage(oomph::DocInfo&) override;
    void doc_element(oomph::DocInfo&) override;


    // Flags
    void setSettling(bool flag) {settlingOn_ = flag;}
    void setAdaptOn(bool flag) {adaptOn_ = flag;}

    // Flags for active interaction forces
    void setInteractionForceP(bool flag) {interactionForceP_ = flag;}
    bool getInteractionForceP() {return interactionForceP_;}
    void setInteractionForceFd(bool flag) {interactionForceFd_ = flag;}
    bool getInteractionForceFd() {return interactionForceFd_;}
    void setInteractionForceShear(bool flag) {interactionForceShear_ = flag;}
    bool getInteractionForceShear() {return interactionForceShear_;}
    void setInteractionForceMagnus(bool flag) {interactionForceMagnus_ = flag;}
    bool getInteractionForceMagnus() {return interactionForceMagnus_;}
    void setInteractionForceSaffman(bool flag) {interactionForceSaffman_ = flag;}
    bool getInteractionForceSaffman() {return interactionForceSaffman_;}

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
    bool incrementalFluidSolve_ = false;
    
    double fluidDensity_;
    double fluidKinematicViscosity_; // At 25 degrees C
    double fluidDynamicViscosity_; // At 25 degrees C

    // Set bools for which forces contribute to interactionForces
    bool interactionForceP_ = false;
    bool interactionForceFd_ = false;
    bool interactionForceShear_ = false;
    bool interactionForceMagnus_ = false;
    bool interactionForceSaffman_ = false;
};

/*!
 * \details General solve routine for a Particle-Fluid coupled simulation (with oomph-lib and MercuryDPM)
 * @param[in] oomph::DocInfo &doc_info, contains info of oomph outputs
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::solveSystem(oomph::DocInfo &doc_info_)
{
    logger(INFO,"Entered solveSystem");

    syncTimeSteps();

    // Creating lists of particles elements and voidages for first timestep
    logger(INFO, "Update lists");
    updateLists();

    logger(INFO, "Completed Update lists");

    if (settlingOn_) settleParticles();

    logger(INFO,"\n\n ---------------------------------------- \n\n             Solving the problem            \n\n ---------------------------------------- \n\n");

    while ( (getTime() < getTimeMax()) || (this->time_pt()->time() < getTimeMax()) )
    {
        if (static_cast<int>(getTime()/getTimeStep()) % dataFile.getSaveCount() == 0)
        {
            logger(INFO, "MercuryDPM/oomph-lib time = % / % , tmax = % ", getTime(), this->time_pt()->time(),
                   getTimeMax());
        }

        //FIXME Hacked in place but should be actionsBeforeTimestep
        passInfoToFluid();

        // If merc is at lower time than oomph, do particle solve
        if (getTime() <= this->time_pt()->time())
        {
            logger(DEBUG, "MercuryDPM/oomph-lib time = % / % , tmax = % ", getTime(), this->time_pt()->time(),
                   getTimeMax());

            logger(DEBUG, "Call to computeOneTimeStep()");
            computeOneTimeStep();
        }

        // If oomph is at lower time than merc, do fluid solve
        else if (getTime() > this->time_pt()->time())
        {
            logger(DEBUG, "MercuryDPM/oomph-lib time = % / % , tmax = % ", getTime(), this->time_pt()->time(),
                   getTimeMax());

            logger(DEBUG, "Call to unsteady_newton_solve()");
            oomph::Problem::unsteady_newton_solve(getTimeStepOomph());
            //doubly_adaptive_unsteady_newton_solve(getTimeStepOomph(),1e-4, 1, false);

            if (!getIncrementalFluidSolve()) {
                if (static_cast<int>(this->time_pt()->time() / getTimeStepOomph()) % getSaveCountOomph() == 0) {
                    doc_solution(doc_info_);
                    doc_voidage(doc_info_);
                    doc_element(doc_info_);
                    doc_info_.number()++;
                }
            }
            else if (getIncrementalFluidSolve())
            {
                logger(INFO,"Call to unsteady_newton_solve() with increasing fraction fo body force coupling");
                int n = 10;
                double incrementBodySolve = 1.0/n;
                for (int i = 0; i <= n; i++)
                {
                    logger(DEBUG,"In bodyForceScaling with i = %, and t = %",i,this->time_pt()->time());
                    setBodyForceFraction(static_cast<double>(i) * incrementBodySolve);
                    logger(INFO,"getBodyForceFraction() = %", getBodyForceFraction());
                    oomph::Problem::unsteady_newton_solve(getTimeStepOomph());

                    logger(INFO,"Call to unsteady_newton_solve() with increasing fraction fo Re");
                    int n2 = 10;
                    double incrementRe = 1.0/n2 * getReynoldsNumber();
                    double newRe = 0.0;
                    for (int i = 0; i <= n2; i++)
                    {
                        setReynoldsNumber(newRe);
                        scaleReynolds();
                        setReynoldsStrouhalNumber(newRe);
                        scaleReynoldsStrouhal();
                        setReynolds_InverseFroudeNumber(newRe);
                        scaleReynolds_InverseFroude();

                        newRe += incrementRe;
                        oomph::Problem::steady_newton_solve(getTimeStepOomph());
                    }
                }
            }

            logger(DEBUG,"Call to updateListOfdVoidagedtInElem()");
            updateListOfdVoidagedtInElem();
            // Updating boundary conditions as it is dependent on time
            actionsAfterTimeStep();

            if (adaptOn_ && (static_cast<int>(this->time_pt()->time()/getTimeStepOomph()) % getAdaptEveryNFluidTimesteps() == 0) )
            {
                logger(DEBUG,"Call to adapt()");
                oomph::Problem::adapt();
                setBC();
            }
        }
    }
};

/*!
 * \details Output for oomph-lib standard output
 * @param[in] oomph::DocInfo &doc_info, contains info of oomph outputs
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::doc_solution(oomph::DocInfo& doc_info)
{
    std::ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts;
    npts=5;

    // Output solution
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();
}

/*!
 * \details Output for oomph-lib voidage (fluid volume fraction)
 * @param[in] oomph::DocInfo &doc_info, contains info of oomph outputs
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::doc_voidage(oomph::DocInfo& doc_info)
{
    std::ofstream some_file;
    char filename[100];

    // number of plot points
    unsigned npts = 1;

    // Output solution if using get_voidage_byEl
    sprintf(filename,"%s/voidagen%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    //mesh_pt()->output_voidage_byEl(some_file,npts); //FIXME Create function
    some_file.close();
}

/*!
 * \details Output for oomph-lib element output, containing nodes and element number
 * @param[in] oomph::DocInfo &doc_info, contains info of oomph outputs
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::doc_element(oomph::DocInfo& doc_info)
{
    std::ofstream some_file;
    char filename[100];

// number of plot points
    unsigned npts = 1;

// Output solution
    sprintf(filename, "%s/elements%i.dat", doc_info.directory().c_str(), doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file, npts);
    some_file.close();
}

/*!
 * \details Syncing time-steps for the two participants based on changing Mercury timestep slightly (MercuryDPM and oomph-lib)
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::syncTimeSteps()
{
    double outputTimeStep = getTimeStepOomph()*getSaveCountOomph();
    double timestep = DPMBase::getTimeStep();
    setSaveCount(std::ceil(outputTimeStep/timestep));
    double mercuryTimeStep = outputTimeStep/dataFile.getSaveCount();
    logger(INFO,"Changing Mercury time step from % to % to output together with oomph (ratio %:1)",getTimeStep(),mercuryTimeStep,dataFile.getSaveCount());
    setTimeStep(mercuryTimeStep);
}

/*!
 * \details Updates lists containing information for coupling during simulation.
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateLists()
{
    updateListOfPInElem();
    updateListOfVoidageInElem();
    updateListOfdVoidagedxInElem();
}

/*!
 * \details Letting the particles settle due to gravity before the coupling simulation is started.
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::settleParticles()
{
    logger(INFO, "Inside settleParticles");

    unsigned it = 0;
    while (getKineticEnergy() > 1e-3 || it < 100)
    {
        it++;
        computeOneTimeStep();

        if ((it % 1000)==0)
        {
            logger(INFO,"getKineticEnergy() = % ",getKineticEnergy());
        }
    }
    logger(INFO,"Particles have settled");
    setTime(0.);
}

/*!
 * \details Actions before adapt, copying tree-structure before it is adapted to keep access to history values
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::actions_before_adapt()
{
    logger(DEBUG,"Call to actions_before_adapt");

    // Make copy before grid is adapted
    prevLeaves_.clear();
    prevTrees_.clear();
    mesh_pt()->forest_pt()->stick_leaves_into_vector(prevLeaves_);
    mesh_pt()->forest_pt()->stick_all_tree_nodes_into_vector(prevTrees_);
    logger(DEBUG,"prevLeaves.size() = %, prevTrees.size() = %",prevLeaves_.size(),prevTrees_.size());
}

/*!
 * \details Actions after adapt, renumber elements, update lists
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::actions_after_adapt()
{
    logger(DEBUG,"Call to actions_after_adapt");

    logger(DEBUG,"Renumbering elements");
    unsigned long nelem=mesh_pt()->nelement();
    for (unsigned e=0;e<nelem;e++)
    {
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(e))->set_number(e); // Adding global number to element
    }

    logger(DEBUG,"We need to store dedt before overwriting");
    auto historyValuesVoidageBeforeAdapt = historyValuesVoidage_;
    auto dVoidagedtBeforeAdapt = listOfdVoidagedtInElem_;

    logger(DEBUG,"Call to generateLists()");
    generateLists();

    logger(DEBUG,"Pass history values for dedt");
    updateListOfdVoidagedtAfterAdapt(historyValuesVoidageBeforeAdapt, dVoidagedtBeforeAdapt);
}

/*!
 * \details Actions after timestep, update drag, pressure and shear contributions to force
 */
//
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::actionsAfterTimeStep()
{
    logger(DEBUG, "Call to actions_after_timestep()");

    updateInteractionForces();

    setBC();
}

/*!
 * \details Updating interaction forces between fluid and particles
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateInteractionForces()
{
    logger(DEBUG,"Call to updateInteractionForces()");
    for (unsigned int iPart = 0; iPart < particleHandler.getNumberOfObjects(); iPart++)
    {
        int pInEl = getPInElemFromList(iPart);
        logger(DEBUG,"pInEl = %",pInEl);

        if (pInEl == -1)
        {
            continue;
        }

        double voidInEl = getVoidageInElemFromList(pInEl);

        logger(DEBUG,"voidInEl = %",voidInEl);

        oomph::Vector<double> s(3,0.0);
        // getS also updates every list if particle moved element
        logger(DEBUG,"Call to getS()");
        getS(s,iPart,pInEl);

        logger(DEBUG,"Call to updatePOnPart()");
        updatePOnPart(iPart, pInEl, voidInEl, s);

        logger(DEBUG,"Call to updateFdOnPart()");
        updateFdOnPart(iPart, pInEl, voidInEl, s);

        logger(DEBUG,"Call to updateShearOnPart()");
        updateShearOnPart(iPart, pInEl, voidInEl, s);

        logger(DEBUG,"Call to updateSaffmanLift()");
        updateSaffmanLift(iPart, pInEl, s);

        logger(DEBUG,"Call to updateMagnusLift()");
        updateMagnusLift(iPart, pInEl, s);
    }

    logger(DEBUG,"Updated forces on all particles");
}

/*!
 * \details Add pressure contribution to force computation
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updatePOnPart(const unsigned int &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s)
{
    if (!getInteractionForceP()) {return;}

    BaseParticle* particlePtr = particleHandler.getObject(iPart);
    Vec3D Fp = getPOnPart(iPart, pInEl, voidInEl, s) * getPressureScaling(); // To get dimensional quantity

    particlePtr->addForce(Fp);
}

/*!
 * \details Computes pressure force on particle based on pressure gradient around particle
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
Vec3D UnderResolvedCoupling<ELEMENT>::getPOnPart(const unsigned int &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s)
{
    BaseParticle* particlePtr = particleHandler.getObject(iPart);
    double Volumei = particlePtr->getVolume();

    Vec3D pGradient = getPressureGradient(pInEl, s);
    Vec3D pForceOnP = -Volumei*(1.-voidInEl) * pGradient; // pForceOnP is non-dimensional as pGradient is non-dimensional

    return pForceOnP;
}

/*!
 * \details Add drag force contribution to force computation
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateFdOnPart(const unsigned int &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s)
{
    if (!getInteractionForceFd()) {return;}

    BaseParticle* particlePtr = particleHandler.getObject(iPart);
    Vec3D dF = getFdOnPart(iPart,pInEl,voidInEl,s) * getBodyForceScaling(); // To get dimensional quantity

    logger(DEBUG,"dF on % = {% % %}",iPart, dF.X, dF.Y, dF.Z);

    particlePtr->addForce(dF);
}

/*!
 * \details Computes drag force on particle based on drag coefficient of particle (Reynolds nr), implemented both Stokes law and Ergun + Wen-Yu
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
Vec3D UnderResolvedCoupling<ELEMENT>::getFdOnPart(const unsigned int &iPart, const int &pInEl, const double &voidInEl, const oomph::Vector<double> &s)
{
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl));

    BaseParticle* particlePtr = particleHandler.getObject(iPart);

    // Drag force is given by: Fd = 1/2 Cd \rho_f \vec{v_rel} |\vec{v_rel}|, the fluid velocity dependence on voidage is already taken care by AJ solver
    double pRad = particlePtr->getRadius();

    oomph::Vector<double> velAtPos(3,0.0);
    elPtr->interpolated_u_nst(s,velAtPos);
    Vec3D scaledFluidVel = getVelocityScaling() * convertVecFuncs::convertToVec3D(velAtPos); // Need to scale fluid velocity to dimensional
    Vec3D relVel = particlePtr->getVelocity() - scaledFluidVel;

    double relVelMagnitude = sqrt(relVel.X*relVel.X + relVel.Y*relVel.Y + relVel.Z*relVel.Z);

    double Re_p = voidInEl * getFluidDensity() * relVelMagnitude * 2.0 * pRad / getFluidDynamicViscosity();
    if (Re_p == 0)
    {
        return {0.,0.,0.};
    }
    logger(DEBUG,"Re_p = %, relVel = (% % %), rho_f = %",Re_p, relVel.X, relVel.Y, relVel.Z, getFluidDensity());
    //https://www.sciencedirect.com/science/article/pii/S1001627917302901?via%3Dihub#bib47 (Valid to Re_p < 1e6)
    // 24./Re + (4.2798.*(0.2.*Re).^0.2374)./(1+(0.2.*Re).^0.8961) + (0.4769.*((Re/263000).^(-7.5968)))./(1+((Re./263000).^(-7.6707))) + (Re.^0.8156)./461000
    double Cd = 24.0 / Re_p + (4.2798 * pow(Re_p / 5.0, 0.2374)) / (1.0 + pow(Re_p / 5.0, 0.8961)) +
                (0.4769 * pow(Re_p / 263000, 0. - 7.5968)) / (1.0 + pow(Re_p / 263000, 0. - 7.6707)) +
                (pow(Re_p, 0.8156)) / (461000);

    // ------------------------------------------------- Stokes law ------------------------------------------------- //
    // compute drag force on particle
    Vec3D dF_stokes = { -(relVel.X) ,
                        -(relVel.Y) ,
                        -(relVel.Z) };
    dF_stokes *= 1./2. * Cd * relVelMagnitude * getFluidDensity() * constants::pi*pRad*pRad;

    dF_stokes /= getBodyForceScaling();
    if (iPart == 0)
    {
        logger(DEBUG, "Re_p = %", Re_p);
        logger(DEBUG, "relVelMagnitude = %, Cd = %", relVelMagnitude, Cd);
        logger(DEBUG, "Using Stokes, Fd(iPart = %) = (%, %, %)",iPart,dF_stokes.X,dF_stokes.Y,dF_stokes.Z);
    }
    //return dF_stokes;
    // ------------------------------------------------- Stokes law ------------------------------------------------- //

    // ---------------------------------------------- Ergun + Wen-Yu ------------------------------------------------ //
    double beta;
    if (Re_p < 1e3)
    {
        Cd = 24./Re_p * (1.0+0.15*pow(Re_p,0.687));
    }
    else
    {
        Cd = 0.44;
    }

    if (voidInEl <= 0.8) // Ergun
    {
        beta = 150.0 * (1.0-voidInEl)*(1.0-voidInEl)/voidInEl * getFluidDynamicViscosity()/(4.0*pRad*pRad) +
               1.75*(1-voidInEl)*getFluidDensity()/(2.0*pRad)*relVelMagnitude;
    }
    else //Wen-Yu
    {
        beta = 0.75*Cd*(voidInEl * (1.0-voidInEl))/(2.0*pRad)*getFluidDensity()*relVelMagnitude*pow(voidInEl,-2.65);
    }

    Vec3D dF = -beta/(1.0-voidInEl) * relVel * (4./3. * constants::pi * pRad * pRad * pRad);
    dF /= getBodyForceScaling();

    if (iPart == 0)
    {
        logger(DEBUG, "Re_p = %", Re_p);
        logger(DEBUG, "beta = %, voidInEl = %", beta, voidInEl);
        logger(DEBUG, "relVelMagnitude = %, Cd = %", relVelMagnitude, Cd);
        logger(DEBUG, "Using Ergun and Wen-Yu, Fd(iPart = %) = (%, %, %)",iPart,dF.X,dF.Y,dF.Z);
    }
    // ---------------------------------------------- Ergun + Wen-Yu ------------------------------------------------ //

    return dF;
}

/*!
 * \details Add shear force contribution to force computation
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
//FIXME Code options open
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateShearOnPart(const unsigned int &iPart_, const int &pInEl_, const double &voidInEl_, const oomph::Vector<double> &s_)
{
    if (!getInteractionForceShear()) {return;}

    BaseParticle* particlePtr = particleHandler.getObject(iPart_);

    oomph::Vector<double> pPos = convertVecFuncs::convertToOomphVec(particlePtr->getPosition());
    double pRad = particlePtr->getRadius();
    double pVolume = particlePtr->getVolume();
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl_));

    oomph::Vector<double> shearForceOnP(3,0.0);
    oomph::Vector<double> torques(3,0.0);

    /// Torque computation
    // Get current vorticity at particle position
    oomph::Vector<double> omega(3,0.0);
    elPtr->get_vorticity(s_, omega);

    Vec3D scaledOmega = convertVecFuncs::convertToVec3D(omega) * getVelocityScaling()/getLengthScaling(); // Need to scale velocity gradient
    scaledOmega *= 0.5;

    /// According to Jonas Einarsson (for creeping flow, incompressible) Why is there an 8, and not a 4?
    Vec3D TorqueMethod1 = 8.0 * constants::pi * getFluidDynamicViscosity() * pRad * pRad * pRad * (scaledOmega - particlePtr->getAngularVelocity() );
    logger(DEBUG,"omega[iPart = %] = (%, %, %)",iPart_, scaledOmega.X, scaledOmega.Y, scaledOmega.Z);
    logger(DEBUG,"particleHandler.getObject(iPart_)->getAngularVelocity()[iPart = %] = (%, %, %)",iPart_, particlePtr->getAngularVelocity().X, particlePtr->getAngularVelocity().Y, particlePtr->getAngularVelocity().Z);
    logger(DEBUG,"TorqueMethod1[iPart = %] = (%, %, %)",iPart_, TorqueMethod1.X, TorqueMethod1.Y, TorqueMethod1.Z);

    particlePtr->addTorque(TorqueMethod1);
    //prevTorques = TorqueMethod1; // No need to set up torques as torque can always be computed unless particle is exactly on element boundary


    /// Shear force contribution
    oomph::Vector<double> dummyVec(3,0.0);
    oomph::GeomObject* geom_obj_pt = nullptr;
    //FIXME We could make setDistFromCOM to relative distance in s, if we find out how to get the absolute distance between the points (see dx_tauxx .../(2*setDistFromCOM))

    oomph::Vector<double> nX(3,0.0); nX[0] = 1.0;
    oomph::Vector<double> nY(3,0.0); nY[1] = 1.0;
    oomph::Vector<double> nZ(3,0.0); nZ[2] = 1.0;
    oomph::Vector<double> nX_neg(3,0.0); nX_neg[0] = -1.0;
    oomph::Vector<double> nY_neg(3,0.0); nY_neg[1] = -1.0;
    oomph::Vector<double> nZ_neg(3,0.0); nZ_neg[2] = -1.0;

//    double setDistFromCOM = pRad;
//    double dx = 2.*setDistFromCOM;
//    double dy = 2.*setDistFromCOM;
//    double dz = 2.*setDistFromCOM;
//    oomph::Vector<double> LocXMin = pPos;   LocXMin[0] -= setDistFromCOM;
//    oomph::Vector<double> LocXMax = pPos;   LocXMax[0] += setDistFromCOM;
//    oomph::Vector<double> LocYMin = pPos;   LocYMin[1] -= setDistFromCOM;
//    oomph::Vector<double> LocYMax = pPos;   LocYMax[1] += setDistFromCOM;
//    oomph::Vector<double> LocZMin = pPos;   LocZMin[2] -= setDistFromCOM;
//    oomph::Vector<double> LocZMax = pPos;   LocZMax[2] += setDistFromCOM;
//
//    oomph::Vector<double> sXMin(3,0.0);
//    oomph::Vector<double> sXMax(3,0.0);
//    oomph::Vector<double> sYMin(3,0.0);
//    oomph::Vector<double> sYMax(3,0.0);
//    oomph::Vector<double> sZMin(3,0.0);
//    oomph::Vector<double> sZMax(3,0.0);
//
//    elPtr->locate_zeta(LocXMin, geom_obj_pt, sXMin); // Takes 1.5e-5 seconds per locate_zeta()
//    elPtr->locate_zeta(LocXMax, geom_obj_pt, sXMax);
//    elPtr->locate_zeta(LocYMin, geom_obj_pt, sYMin);
//    elPtr->locate_zeta(LocYMax, geom_obj_pt, sYMax);
//    elPtr->locate_zeta(LocZMin, geom_obj_pt, sZMin);
//    elPtr->locate_zeta(LocZMax, geom_obj_pt, sZMax);

    // Using estimate for s instead of locate_zeta (speed up of factor ~14 without locate_zeta()), see above code
    double setDistFromCOM_in_s = 0.01;
    oomph::Vector<double> sXMin = s_;
    oomph::Vector<double> sXMax = s_;
    oomph::Vector<double> sYMin = s_;
    oomph::Vector<double> sYMax = s_;
    oomph::Vector<double> sZMin = s_;
    oomph::Vector<double> sZMax = s_;
    sXMin[0] -= setDistFromCOM_in_s;
    sXMax[0] += setDistFromCOM_in_s;
    sYMin[1] -= setDistFromCOM_in_s;
    sYMax[1] += setDistFromCOM_in_s;
    sZMin[2] -= setDistFromCOM_in_s;
    sZMax[2] += setDistFromCOM_in_s;
    oomph::Vector<double> LocXMin(3,0.0);
    oomph::Vector<double> LocXMax(3,0.0);
    oomph::Vector<double> LocYMin(3,0.0);
    oomph::Vector<double> LocYMax(3,0.0);
    oomph::Vector<double> LocZMin(3,0.0);
    oomph::Vector<double> LocZMax(3,0.0);
    elPtr->interpolated_x(sXMax,LocXMax);
    elPtr->interpolated_x(sXMin,LocXMin);
    elPtr->interpolated_x(sYMax,LocYMax);
    elPtr->interpolated_x(sYMin,LocYMin);
    elPtr->interpolated_x(sZMax,LocZMax);
    elPtr->interpolated_x(sZMin,LocZMin);
    double dx = LocXMax[0] - LocXMin[0];
    double dy = LocYMax[1] - LocYMin[1];
    double dz = LocZMax[2] - LocZMin[2];
    // Using estimate for s instead of locate_zeta (speed up of factor ~14 without locate_zeta())

    // Derivatives w.r.t x-dir
    oomph::Vector<double> tnXMin(3,0.0);
    oomph::Vector<double> ttXMin(3,0.0);
    oomph::Vector<double> tnXMax(3,0.0);
    oomph::Vector<double> ttXMax(3,0.0);

    elPtr->get_traction(sXMin,nX_neg,dummyVec,tnXMin,ttXMin);
    elPtr->get_traction(sXMax,nX,dummyVec,tnXMax,ttXMax);
    double tau_xx_XMin = tnXMin[0];
    double tau_xy_XMin = ttXMin[1];
    double tau_xz_XMin = ttXMin[2];
    double tau_xx_XMax = tnXMax[0];
    double tau_xy_XMax = ttXMax[1];
    double tau_xz_XMax = ttXMax[2];
    double dx_tauxx = (tau_xx_XMax + tau_xx_XMin)/(2.*dx); // If nX
    double dx_tauxy = (tau_xy_XMax + tau_xy_XMin)/(2.*dx); // If nX
    double dx_tauxz = (tau_xz_XMax + tau_xz_XMin)/(2.*dx); // If nX


    // Derivatives w.r.t y-dir
    oomph::Vector<double> tnYMin(3,0.0);
    oomph::Vector<double> ttYMin(3,0.0);
    oomph::Vector<double> tnYMax(3,0.0);
    oomph::Vector<double> ttYMax(3,0.0);

    elPtr->get_traction(sYMin,nY_neg,dummyVec,tnYMin,ttYMin);
    elPtr->get_traction(sYMax,nY,dummyVec,tnYMax,ttYMax);
    double tau_yx_YMin = ttYMin[0];
    double tau_yy_YMin = tnYMin[1];
    double tau_yz_YMin = ttYMin[2];
    double tau_yx_YMax = ttYMax[0];
    double tau_yy_YMax = tnYMax[1];
    double tau_yz_YMax = ttYMax[2];
    double dy_tauyx = (tau_yx_YMax + tau_yx_YMin)/(2.*dy); // If nY
    double dy_tauyy = (tau_yy_YMax + tau_yy_YMin)/(2.*dy); // If nY
    double dy_tauyz = (tau_yz_YMax + tau_yz_YMin)/(2.*dy); // If nY


    // Derivatives w.r.t z-dir
    oomph::Vector<double> tnZMin(3,0.0);
    oomph::Vector<double> ttZMin(3,0.0);
    oomph::Vector<double> tnZMax(3,0.0);
    oomph::Vector<double> ttZMax(3,0.0);

    elPtr->get_traction(sZMin,nZ_neg,dummyVec,tnZMin,ttZMin);
    elPtr->get_traction(sZMax,nZ,dummyVec,tnZMax,ttZMax);
    double tau_zx_ZMin = ttZMin[0];
    double tau_zy_ZMin = ttZMin[1];
    double tau_zz_ZMin = tnZMin[2];
    double tau_zx_ZMax = ttZMax[0];
    double tau_zy_ZMax = ttZMax[1];
    double tau_zz_ZMax = tnZMax[2];
    double dz_tauzx = (tau_zx_ZMax + tau_zx_ZMin)/(2.*dz); // If nZ
    double dz_tauzy = (tau_zy_ZMax + tau_zy_ZMin)/(2.*dz); // If nZ
    double dz_tauzz = (tau_zz_ZMax + tau_zz_ZMin)/(2.*dz); // If nZ


    // Shear force contribution is equal to (1-eps)*V_p * nabla dot tau
    shearForceOnP[0] = (1.0-voidInEl_)*(pVolume)* (dx_tauxx + dy_tauyx + dz_tauzx);
    shearForceOnP[1] = (1.0-voidInEl_)*(pVolume)* (dx_tauxy + dy_tauyy + dz_tauzy);
    shearForceOnP[2] = (1.0-voidInEl_)*(pVolume)* (dx_tauxz + dy_tauyz + dz_tauzz);

    Vec3D scaledShearForce = convertVecFuncs::convertToVec3D(shearForceOnP); // Need to scale shear force
    scaledShearForce *= getLengthScaling() * getLengthScaling() / getVelocityScaling();

    particlePtr->addForce(scaledShearForce);
    logger(DEBUG,"shearForceOnP[iPart = %] = (%, %, %)",iPart_, shearForceOnP[0], shearForceOnP[1], shearForceOnP[2]);

    //allPrevTorques += convertVecFuncs::convertToVec3D(torques);
    //allPrevForces += convertVecFuncs::convertToVec3D(shearForceOnP);
}

/*!
 * \details Add Saffman lift contribution to force computation
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const double &voidInEl, local voidage
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateSaffmanLift(const unsigned &iPart_, const int &pInEl_, const oomph::Vector<double> &s_)
{
    if (!getInteractionForceSaffman()) {return;}

    BaseParticle* particlePtr = particleHandler.getObject(iPart_);

    oomph::Vector<double> velAtPos(3,0.0);
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl_));

    elPtr->interpolated_u_nst(s_,velAtPos);
    Vec3D scaledFluidVel = getVelocityScaling() * convertVecFuncs::convertToVec3D(velAtPos); // Need to scale fluid velocity to dimensional
    Vec3D relVel = scaledFluidVel - particlePtr->getVelocity();

    if (relVel.X == 0 && relVel.Y == 0 && relVel.Z == 0)
    {
        return;
    }

    double pRad = particlePtr->getRadius();

    // interpolated_dudx_nst(s_,u[i],x[j])
    double dudy = elPtr->interpolated_dudx_nst(s_,0,1);
    double dudz = elPtr->interpolated_dudx_nst(s_,0,2);
    double dvdx = elPtr->interpolated_dudx_nst(s_,1,0);
    double dvdz = elPtr->interpolated_dudx_nst(s_,1,2);
    double dwdx = elPtr->interpolated_dudx_nst(s_,2,0);
    double dwdy = elPtr->interpolated_dudx_nst(s_,2,1);
    double sign_dudy = mathsFunc::sign(dudy);
    double sign_dudz = mathsFunc::sign(dudz);
    double sign_dvdx = mathsFunc::sign(dvdx);
    double sign_dvdz = mathsFunc::sign(dvdz);
    double sign_dwdx = mathsFunc::sign(dwdx);
    double sign_dwdy = mathsFunc::sign(dwdy);

    Vec3D Fsaffman = {relVel.Y*sign_dvdx*std::pow(std::abs(dvdx),0.5) + relVel.Z*sign_dwdx*std::pow(std::abs(dwdx),0.5),
                      relVel.X*sign_dudy*std::pow(std::abs(dudy),0.5) + relVel.Z*sign_dwdy*std::pow(std::abs(dwdy),0.5),
                      relVel.X*sign_dudz*std::pow(std::abs(dudz),0.5) + relVel.Y*sign_dvdz*std::pow(std::abs(dvdz),0.5)};
    Fsaffman *= getVelocityScaling()/getLengthScaling(); // Need to scale velocity gradient to dimensional

    Fsaffman *= 1.615*getFluidDensity()*std::pow(getFluidDynamicViscosity(),0.5)*4*pRad*pRad;


    particlePtr->addForce(Fsaffman);
    logger(DEBUG,"relVel[iPart = %] = (%, %, %)",iPart_, relVel.X, relVel.Y, relVel.Z);
    logger(DEBUG,"Fsaffman[iPart = %] = (%, %, %)",iPart_, Fsaffman.X, Fsaffman.Y, Fsaffman.Z);

    //allPrevForces += Fsaffman;
}

/*!
 * \details Add Magnus lift contribution to force computation
 * @param[in] const unsigned int &iPart, particle index
 * @param[in] const int &pInEl, element index containing iPart
 * @param[in] const oomph::Vector<double> &s, parameterised location of iPart in pInEl
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateMagnusLift(const unsigned &iPart_, const int &pInEl_, const oomph::Vector<double> &s_)
{
    if (!getInteractionForceMagnus()) {return;}

    BaseParticle* particlePtr = particleHandler.getObject(iPart_);

    https://www.sciencedirect.com/science/article/pii/S1001627917302901, Issues in Eulerian–Lagrangian modeling of sediment transport under saltation regime
    // Fm = 1/2 rho_f Cl A |u_rel| (u_rel) cross (omega_rel / |omega_rel|)
    // Cl = min(0.5, 0.5pRad |w_rel|/|u_rel|)
    Vec3D FMagnus ={0.0, 0.0, 0.0};
    double pRad =particlePtr->getRadius();
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl_));

    oomph::Vector<double> velAtPos(3,0.0);
    elPtr->interpolated_u_nst(s_,velAtPos);
    Vec3D scaledFluidVel = getVelocityScaling() * convertVecFuncs::convertToVec3D(velAtPos); // Need to scale fluid velocity to dimensional
    Vec3D relVel = scaledFluidVel - particlePtr->getVelocity();

    Vec3D absRelVel;
    absRelVel.X = abs(relVel.X);
    absRelVel.Y = abs(relVel.Y);
    absRelVel.Z = abs(relVel.Z);
    double magnitudeAbsRelVel = sqrt(pow(absRelVel.X,2) + pow(absRelVel.Y,2) + pow(absRelVel.Z,2));

    oomph::Vector<double> omega(3,0.0);
    elPtr->get_vorticity(s_, omega);

    Vec3D scaledOmega = convertVecFuncs::convertToVec3D(omega) * getVelocityScaling()/getLengthScaling(); // Need to scale velocity to dimensional
    scaledOmega *= 0.5;

    Vec3D relRot = scaledOmega - particlePtr->getAngularVelocity();
    Vec3D absRelRot;
    absRelRot.X = abs(relRot.X);
    absRelRot.Y = abs(relRot.Y);
    absRelRot.Z = abs(relRot.Z);
    double magnitudeAbsRelRot = sqrt(pow(absRelRot.X,2) + pow(absRelRot.Y,2) + pow(absRelRot.Z,2));

    if (magnitudeAbsRelRot != 0.0)
    {
        double Cl = std::min(0.5, 0.5 * pRad * (magnitudeAbsRelRot / magnitudeAbsRelVel));

        FMagnus = Vec3D::cross(relVel, relRot);

        if (absRelRot.X != 0.0)
        {
            FMagnus.X *= 0.5 * getFluidDensity() * constants::pi * Cl * pRad * pRad * absRelVel.X / absRelRot.X;
        }
        else
        {
            FMagnus.X = 0;
        }

        if (absRelRot.Y != 0.0)
        {
            FMagnus.Y *= 0.5 * getFluidDensity() * constants::pi * Cl * pRad * pRad * absRelVel.Y / absRelRot.Y;
        }
        else
        {
            FMagnus.Y = 0;
        }

        if (absRelRot.Z != 0.0)
        {
            FMagnus.Z *= 0.5 * getFluidDensity() * constants::pi * Cl * pRad * pRad * absRelVel.Z / absRelRot.Z;
        }
        else
        {
            FMagnus.Z = 0;
        }
    }
    particlePtr->addForce(FMagnus);
    logger(DEBUG,"FMagnus[iPart = %] = (%, %, %)",iPart_, FMagnus.X, FMagnus.Y, FMagnus.Z);

    //allPrevForces += FMagnus;
}

template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::passInfoToFluid()
{
    for (int i = 0; i<mesh_pt()->nelement(); i++)
    {
        auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(i));

        if (!elPtr->voidageIsFixed())
        {
            double localVoidage = getVoidageInElemFromList(i);
            oomph::Vector<double> dummyF = {0.,0.,0.};

            elPtr->setVoidage(localVoidage);
            getBodyForceInElemByCoupling(i, dummyF);
        }
    }
}

template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::getBodyForceInElemByCoupling(const int &iEl_, oomph::Vector<double> &force_)
{
    Vec3D Ftotal = {0.,0.,0.};
    for (int iP = 0; iP < particleHandler.getNumberOfObjects(); iP++)
    {
        if (getPInElemFromList(iP) == iEl_)
        {
            int pInEl = getPInElemFromList(iP);
            double voidInEl = getVoidageInElemFromList(pInEl);
            oomph::Vector<double> s(3,0.0);
            getS(s,iP,pInEl);

            Ftotal += getPOnPart(iP, pInEl, voidInEl, s);
            Ftotal += getFdOnPart(iP, pInEl, voidInEl, s);
            //FIXME ADD other 3 forces? they are rotations though, not normal forces so as vorticity?

            logger(DEBUG,"Ftotal on % = {% % %}",iP, Ftotal.X, Ftotal.Y, Ftotal.Z);
        }
    }

    oomph::Vector<double> posNode0(3, 0.0);
    oomph::Vector<double> posNodeEnd(3, 0.0);
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(iEl_));
    const unsigned int iNode = elPtr->nnode();
    // Note that index of end node is iNode - 1, else segmentation fault as node_pt(iNode) does not exist
    elPtr->node_pt(0)->position(posNode0);
    elPtr->node_pt(iNode - 1)->position(posNodeEnd);
    double VolumeElement = abs(posNode0[0] - posNodeEnd[0]) * abs(posNode0[1] - posNodeEnd[1]) * abs(posNode0[2] - posNodeEnd[2]);
    logger(DEBUG,"Ftotal on el % = {% % %}", iEl_, Ftotal.X, Ftotal.Y, Ftotal.Z);

    Ftotal /= VolumeElement;
    Ftotal *= getBodyForceFraction();
    //Ftotal *= (getVelocityScaling()*getVelocityScaling())/(getLengthScaling()*getFluidDynamicViscosity()); // To get right scaling with oomph-lib as they divide by L^2/U mu_ref

    Ftotal *= -1.0; // As action = -reaction

    if (Ftotal.Z > 1e-15)
    {
        logger(INFO, "Ftotal on el % = {% % %}\n", iEl_, Ftotal.X, Ftotal.Y, Ftotal.Z);
    }

    oomph::Vector<double> F = convertVecFuncs::convertToOomphVec(Ftotal);

    force_ = F; //FIXME where is this force added to the residuals?

    elPtr->setCouplingForce(F);
}

template <class ELEMENT>
Vec3D UnderResolvedCoupling<ELEMENT>::getPressureGradient(const int & pIsInEl_, const oomph::Vector<double> &s_)
{
    oomph::Vector<double> pGrad(3, 0.0);
    auto elPtr = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pIsInEl_));

    pGrad[0] = elPtr->interpolated_dpdx_nst(s_,0);
    pGrad[1] = elPtr->interpolated_dpdx_nst(s_,1);
    pGrad[2] = elPtr->interpolated_dpdx_nst(s_,2);

    return convertVecFuncs::convertToVec3D(pGrad);
}


/*!
 * \details Prescribe voidage values in certain elements that have volumetric overlap with input values in a ceratin direction
 *
 * \param[in] double voidageValue prescribed value for the voidage
 * \param[in] double from Start of domain that has prescribed voidage
 * \param[in] double to End of domain that has prescribed voidage
 * \param[in] double dir Direction of the domain
 *
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::prescribeVoidage(double voidageValue, double from, double to, int dir)
{
    for (unsigned int iEl = 0; iEl < mesh_pt()->nelement(); iEl++)
    {
        auto elPtr = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(iEl));

        oomph::Vector<double> posNode0(3, 0.0);
        oomph::Vector<double> posNodeEnd(3, 0.0);
        const unsigned int iNode = elPtr->nnode();
        elPtr->node_pt(0)->position(posNode0);
        elPtr->node_pt(iNode - 1)->position(posNodeEnd);

        if ((posNode0[dir]>from && posNode0[dir]<to) || (posNodeEnd[dir]>from && posNodeEnd[dir]<to))
        {
            elPtr->setVoidage(voidageValue);
            elPtr->setVoidageFixed(true);
        }
    }
}

/*!
 * \details Create lists as initial condition.
 *
 * \param[in]
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::generateLists()
{
    setListsLengths();
    updateListOfPInElem();
    updateListOfVoidageInElem();

    updateNeighbourList();
    updateListOfdVoidagedxInElem();
}

/*!
 * \details Clear and update 2 lists with default values.
 * 1) listOfPInElem, Contains a list with particleID and elementID it is in.
 * 2) listOfVoidageInElem, Contains a list of elementID with its voidage.
 *
 * \param[in]
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::setListsLengths()
{
    logger(DEBUG,"Call to setListsLengths()");

    // Resizing the vectors
    listOfPInElem_.resize(particleHandler.getNumberOfObjects());
    listOfVoidageInElem_.resize(mesh_pt()->nelement());
    listOfdVoidagedxInElem_.resize(mesh_pt()->nelement());
    listOfdVoidagedtInElem_.resize(mesh_pt()->nelement());
    historyValuesVoidage_.resize(mesh_pt()->nelement());

    // Setting default values
    for (unsigned int iP = 0; iP < particleHandler.getNumberOfObjects(); iP++)
    {
        listOfPInElem_[iP] = -1;
    }

    for (unsigned int iEl = 0; iEl < mesh_pt()->nelement(); iEl++)
    {

        listOfVoidageInElem_[iEl] = 1.0;
        listOfdVoidagedtInElem_[iEl] = 0.0;

        listOfdVoidagedxInElem_[iEl].resize(3);
        for (unsigned int iDir = 0; iDir < 3; iDir++)
        {
            listOfdVoidagedxInElem_[iEl][iDir] = 0.0;
        }

        historyValuesVoidage_[iEl].resize(5);
        for (unsigned int iHist = 0; iHist < 5; iHist++)
        {
            historyValuesVoidage_[iEl][iHist] = 1.0;
        }
    }
}

/*!
 * \details update lists for specific particle upon moving elements. It determines elements to update based on
 * listOfPInElem.
 *
 * \param[in] const unsgined &iPart, particle id
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListsForSpecificParticle(const unsigned &iPart)
{
    // Make backup of previous element number
    int prevEl = getPInElemFromList(iPart);

    // Update the listOfPInElem[iPart]
    updateListOfPInElem(iPart);

    // Get current element number
    int currEl = getPInElemFromList(iPart);
    logger(DEBUG,"prevEl = %, currEl = %", prevEl, currEl);

    // Update both fluid element voidages as the particle moved boundary
    updateListOfVoidageInElem(prevEl);
    updateListOfVoidageInElem(currEl);

    logger(DEBUG,"updateListOfVoidageInElem complete");

    // Update list for dvoidagedx for specific elements
    updateListOfdVoidagedxAroundElem(prevEl);
    updateListOfdVoidagedxAroundElem(currEl);

    //Note that dVoidagedt is updated only after a fluid timestep
    logger(DEBUG,"Updated elements % & % (and its neighbours for depsdx)",prevEl, currEl);
}

/*!
 * \details update listOfPInElem, Contains a list with particleID and elementID it is in.
 *
 * \param[in]
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfPInElem()
{
    for (unsigned int iP = 0; iP < particleHandler.getNumberOfObjects(); iP++)
    {
        updateListOfPInElem(iP);
    }
}


/*!
 * \details update listOfPInElem for specific particle.
 *
 * \param[in] const unsigned iPart, particle id
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfPInElem(const unsigned &iPart)
{
    logger(INFO, "Inside updateListOfPInElem(const unsigned &iPart)");

    Vec3D pPos = particleHandler.getObject(iPart)->getPosition();
    const unsigned long nEl = mesh_pt()->nelement();
    // Loop over fluid elements
    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        oomph::Vector<double> s(3, 0.0);
        oomph::GeomObject* geom_obj_pt = nullptr;
        dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(iEl))->locate_zeta(convertVecFuncs::convertToOomphVec(pPos), geom_obj_pt, s);
        if (geom_obj_pt != nullptr)
        {
            listOfPInElem_[iPart] = iEl;
            break; // Break out of fluid-elements loop when particle was found in cell
        }
    }
}

/*!
 * \details update listOfVoidageInElem, Contains a list of elementID with its voidage.
 * This function uses the information stored in listOfPInElem
 *
 * \param[in]
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfVoidageInElem()
{
    logger(DEBUG,"Call to updateListOfVoidageInElem()");

    double nEl = mesh_pt()->nelement();

    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        updateListOfVoidageInElem(iEl);
    }
}

/*!
 * \details update listOfVoidageInElem, for specific fluid element
 *
 * \param[in] const unsigned iEl, fluid element id
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfVoidageInElem(const unsigned &iEl)
{
    auto elPtr = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(iEl));

    if (elPtr->voidageIsFixed()) {return;}

    double VolumeElement;
    double VolumeParticles = 0.0;
    oomph::Vector<double> posNode0(3, 0.0);
    oomph::Vector<double> posNodeEnd(3, 0.0);

    const unsigned int iNode = elPtr->nnode();

    // Note that index of end node is iNode - 1, else segmentation fault as node_pt(iNode) does not exist
    elPtr->node_pt(0)->position(posNode0);
    elPtr->node_pt(iNode - 1)->position(posNodeEnd);

    VolumeElement = abs(posNode0[0] - posNodeEnd[0]) * abs(posNode0[1] - posNodeEnd[1]) * abs(posNode0[2] - posNodeEnd[2]);

    for (unsigned int iP = 0; iP < listOfPInElem_.size(); iP++)
    {
        if (listOfPInElem_[iP] == iEl)
        {
            VolumeParticles += particleHandler.getObject(iP)->getVolume();
        }
    }
    listOfVoidageInElem_[iEl] =  (1. - VolumeParticles / VolumeElement);
}

/*!
 * \details update listOfNeighbours
 *
 * \param[in]
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateNeighbourList()
{
    logger(INFO,"Call to updateNeighbourList()");
    int startInd = 20; // Used to access right index in enum

    //FIXME Make output code new function with doc_info as input var
    /*
        //output File for voidage containing information on what el is where
        std::ofstream some_file;
        char filename[100];

        // number of plot points
        unsigned npts = 1;

        // Output solution if using get_voidage_byEl
        sprintf(filename,"%s/updateNeighbourListVoidagen%i.dat",doc_info.directory().c_str(),doc_info.number());
        some_file.open(filename);
        //mesh_pt()->output_voidage_byEl(some_file,npts);
        some_file.close();
        doc_info.number()++;
    */

    unsigned int nTree = mesh_pt()->forest_pt()->ntree();
    logger(DEBUG,"ntree = %",nTree);

    unsigned int minLevel, maxLevel;
    mesh_pt()->get_refinement_levels(minLevel,maxLevel);
    logger(DEBUG,"minRefLevel %, maxRefLevel %", minLevel, maxLevel);

    unsigned long nEl = mesh_pt()->nelement();

    //Empty existing listOfRefinableNeighbours
    listOfRefinableNeighbours_.clear();

    // Reset sizing of listOfRefinableNeighbours
    listOfRefinableNeighbours_.resize(nEl);
    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        listOfRefinableNeighbours_[iEl].resize(6);
        for (int i = 0; i<6; i++)
        {
            listOfRefinableNeighbours_[iEl][i].push_back(nullptr);
        }
    }

    logger(DEBUG,"Going into loop, nTree = %", nTree);
    for (unsigned iTree = 0; iTree < nTree; iTree++)
    {
        oomph::Tree* treeRootPtr = mesh_pt()->forest_pt()->tree_pt(iTree);

        oomph::Vector<oomph::Tree*> allLeavesOfTree;
        treeRootPtr->stick_leaves_into_vector(allLeavesOfTree);

        logger(DEBUG,"allLeavesOfTree.size() = %",allLeavesOfTree.size());

        if (allLeavesOfTree.size() == 1) // allLeavesOfTree only contains the treeRootPtr, iow there are no other leaves except the treeRootPtr
        {
            for (int dir = 0; dir < 6; dir++)
            {
                oomph::Tree* neighbourTreeRootPtr = mesh_pt()->forest_pt()->tree_pt(iTree)->neighbour_pt(startInd+dir);
                if (neighbourTreeRootPtr == nullptr)
                {
                    logger(DEBUG,"There is no neighbour in  dir %",dir);
                    // No need to pushback anything as there is no neighbour in this direction
                }
                else
                {
                    logger(DEBUG,"There is a neighbour in  dir %",dir);

                    oomph::Vector<oomph::Tree*> neighbourTreeRootPtrsLeaves;

                    stickLeavesOfEqualOrSmallerSizedNeighbourRecursiveInDir(neighbourTreeRootPtr, neighbourTreeRootPtrsLeaves, dir);
                    for (auto iLeaf : neighbourTreeRootPtrsLeaves)
                    {
                        listOfRefinableNeighbours_[treeRootPtr->object_pt()->number()][dir].push_back(iLeaf); //Pushback the other treeRootPtr
                    }
                }
            }
        }
        else
        {
            for (auto iLeaf : allLeavesOfTree)
            {
                for (int dir = 0; dir < 6; dir++)
                {
                    oomph::Vector<oomph::Tree*> neighbourTreePtrs;
                    logger(DEBUG, "We get here");
                    stickLeavesOfNeighbourRecursiveInDir(iLeaf, iTree, neighbourTreePtrs, dir);
                    for (auto iLeaf2 : neighbourTreePtrs)
                    {
                        logger(DEBUG,"Pushed back %", iLeaf2->object_pt()->number());
                        listOfRefinableNeighbours_[iLeaf->object_pt()->number()][dir].push_back(iLeaf2); //Pushback the other treeRootPtr
                    }
                }
            }
        }
    }

    outputNeighboursOfElements();
}

/*!
 * \details creates a vector of Tree* that are neighbouring to input Tree* in a certain direction.
 * It needs iTree to find the neighbours on TreeRoot level (see getExternalNeighbour)
 *
 * @param[in] Tree* Tree_, pointer to Tree of which neighbours need to be determined
 * @param[in] const int &iTree, index of the Tree that Tree_ is part of
 * @param[in] std::vector<Tree*> & allNeighbouringLeaves_, vector in which all neighbouring leaves will be stored
 * @param[in] const int &dir_, direction in which the elements neighbour needs to be determined following the oomph-lib
 *              numbering convention in direction (-x = 0, +x = 1, ... , +z = 5)
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::stickLeavesOfNeighbourRecursiveInDir(oomph::Tree* Tree_, const int &iTree, std::vector<oomph::Tree*> &allNeighbouringLeaves_, const int &dir_)
{
    logger(DEBUG,"call to stickLeavesOfNeighbourRecursiveInDir for iEl % in dir %",Tree_->object_pt()->number(),dir_);
    logger(DEBUG,"Tree_->level() == %",Tree_->level());

    int indexTree_ = -1; // = Tree_->son_type();

    getInternalIndex(Tree_, indexTree_);
    logger(DEBUG,"indexTree_ = %", indexTree_);

    int indexInternalNeighbour = 0;
    int differenceInLevel = 0;
    std::vector<int> indexPath;

    getNeighbourIndex(indexTree_, dir_, indexInternalNeighbour);
    oomph::Tree* neighbourPtr = nullptr;

    logger(DEBUG,"indexInternalNeighbour = %", indexInternalNeighbour);

    if (indexInternalNeighbour == -1) // Neighbour can not be found internally
    {
        logger(DEBUG,"Try to find neighbour externally");
        getExternalNeighbour(Tree_, iTree, dir_ , differenceInLevel, indexPath, neighbourPtr); // Ignore for now that the neighbour can have sons
    }
    else
    {
        neighbourPtr = Tree_->father_pt()->son_pt(indexInternalNeighbour);
    }

    logger(DEBUG,"neighbourPtr found");
    if (neighbourPtr == nullptr)
    {
        logger(DEBUG,"No neighbour in this direction %",dir_);
        return;
        // No need to push anything pack
    }

    // Make sure that only the right sons are considered
    oomph::Tree* neighbourTreeSonAtSameLevel;
    neighbourTreeSonAtSameLevel = neighbourPtr;

    logger(DEBUG,"neighbourTreeSonAtSameLevel->number = %", neighbourTreeSonAtSameLevel->object_pt()->number());
    logger(DEBUG,"neighbourTreeSonAtSameLevel->nsons() = %", neighbourTreeSonAtSameLevel->nsons());
    logger(DEBUG,"differenceInLevel = %", differenceInLevel);

    while (neighbourTreeSonAtSameLevel->nsons() > 0 && differenceInLevel > 0)
    {
        int indexInNeighbourFromPath = convertIndexFromPath(indexPath, differenceInLevel, dir_);
        neighbourTreeSonAtSameLevel = neighbourTreeSonAtSameLevel->son_pt(indexInNeighbourFromPath);

        differenceInLevel--;
    }
    if (differenceInLevel >= 0) // Neighbour is equal sized, but it may still have sons
    {
        stickLeavesOfEqualOrSmallerSizedNeighbourRecursiveInDir(neighbourTreeSonAtSameLevel, allNeighbouringLeaves_, dir_);
    }
    else // Neighbour is bigger, so it is the only adjacent element
    {
        allNeighbouringLeaves_.push_back(neighbourTreeSonAtSameLevel);
    }
}

/*!
 * \details Find the internal index of Tree_ in its father, following oomph-lib number convention.
 * If the element has no father (e.g. Tree_ is a TreeRoot*), internalIndex is assigned -1
 *
 * @param[in] Tree* Tree_, pointer to Tree of which neighbours need to be determined
 * @param[in] &internalIndex
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::getInternalIndex(oomph::Tree* Tree_, int &internalIndex)
{
    logger(DEBUG,"call to getInternalIndex");

    oomph::Tree* fathPtr = Tree_->father_pt();
    if (fathPtr == nullptr) //
    {
        logger(DEBUG,"father Pointer = nullptr");
    }

    internalIndex = -1;
    std::vector<int> elNumbers;
    int elNumber;

    for (unsigned int i = 0; i<fathPtr->nsons(); i++)
    {
        elNumber = fathPtr->son_pt(i)->object_pt()->number();
        if (elNumber == Tree_->object_pt()->number())
        {
            internalIndex = i;
            return;
        }
    }

    if (internalIndex == -1)
    {
        logger(INFO,"internalIndex could not be found");
    }
}

/*!
 * \details Assigns index of internal neighbour to indexInternalNeighbour. If no internal neighbour is found, assign -1.
 *
 * @param[in] const int indexTree_, index of a Tree in its father
 * @param[in] const int dir_, the direction to check
 * @param[in] int &indexInternalNeighbour_, contains index of neighbour.
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::getNeighbourIndex(const int indexTree_, const int dir_, int &indexInternalNeighbour_)
{
    logger(DEBUG,"call to getNeighbourIndex");
    logger(DEBUG,"indexTree_ = %, dir_ = %, indexInternalNeighbour_ = %", indexTree_, dir_, indexInternalNeighbour_);

    std::vector<int> internalNeighbour; // Neighbour, index corresponds to dir
    if (indexTree_ == 0)
    {
        internalNeighbour = {-1, 1, -1, 2, -1, 4};
    }
    else if (indexTree_ == 1)
    {
        internalNeighbour = {0, -1, -1, 3, -1, 5};
    }
    else if (indexTree_ == 2)
    {
        internalNeighbour = {-1, 3, 0, -1, -1, 6};
    }
    else if (indexTree_ == 3)
    {
        internalNeighbour = {2, -1, 1, -1, -1, 7};
    }
    else if (indexTree_ == 4)
    {
        internalNeighbour = {-1, 5, -1, 6, 0, -1};
    }
    else if (indexTree_ == 5)
    {
        internalNeighbour = {4, -1, -1, 7, 1, -1};
    }
    else if (indexTree_ == 6)
    {
        internalNeighbour = {-1, 7, 4, -1, 2, -1};
    }
    else if (indexTree_ == 7)
    {
        internalNeighbour = {6, -1, 5, -1, 3, -1};
    }

    indexInternalNeighbour_ = internalNeighbour[dir_];
    logger(DEBUG,"indexTree_ = %, dir_ = %, indexInternalNeighbour_ = %", indexTree_, dir_, indexInternalNeighbour_);

}

/*!
 * \details Find the neighbour of Tree_, use when no internal neighbour could be found. neighbourPtr points to an equal
 * or larger sized element compared to Tree_ (i.o.w. neighbourPtr->Level() <= Tree_->Level())
 *
 * @param[in] Tree* Tree_, pointer to Tree of which neighbours need to be determined
 * @param[in] const int &iTree, index of the Tree that Tree_ is part of
 * @param[in] const int &dir_, direction in which the elements neighbour needs to be determined following the oomph-lib
 *              numbering convention in direction (-x = 0, +x = 1, ... , +z = 5)
 * @param[in] int &differenceInLevel, difference in level between Tree_ and neigbourPtr_
 * @param[in] std::vector<int> &indexPath, contain indices of the path taken to find the neighbour of Tree_
 * @param[in] Tree* &neighbourPtr_, pointer to neighbouring Tree
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::getExternalNeighbour(oomph::Tree* Tree_, const int &iTree, const int dir_, int &differenceInLevel, std::vector<int> &indexPath, oomph::Tree* &neighbourPtr_)
{
    //std::cout<<std::endl;
    logger(DEBUG,"call to getExternalNeighbour for el % in idir %, iTree %,", Tree_->object_pt()->number(),dir_, iTree);
    int startInd = 20;
    differenceInLevel++;

    int internalIndex;
    getInternalIndex(Tree_, internalIndex);
    indexPath.push_back(internalIndex);
    logger(DEBUG,"internalIndex = %", internalIndex);


    int indexFather;
    if (Tree_->father_pt() == Tree_->root_pt())
    {
        logger(DEBUG,"Tree_->father_pt() == treeRootPtr, hence make neighbourPtr other tree");
        oomph::Tree* externalNeighbour = mesh_pt()->forest_pt()->tree_pt(iTree)->neighbour_pt(startInd+dir_);
        neighbourPtr_ = externalNeighbour;
        return;
    }
    else
    {
        logger(DEBUG,"(Tree_->father_pt() != nullptr)");
    }

    getInternalIndex(Tree_->father_pt(),indexFather);
    logger(DEBUG,"indexFather = %", indexFather);

    int indexFatherNeighbour;


    getNeighbourIndex(indexFather, dir_, indexFatherNeighbour);

    logger(DEBUG,"indexFather = %, dir_ = %, indexFatherNeighbour = %", indexFather, dir_, indexFatherNeighbour);

    if (indexFatherNeighbour == -1)
    {
        getExternalNeighbour(Tree_->father_pt(), iTree, dir_, differenceInLevel, indexPath, neighbourPtr_);
    }
    else
    {
        neighbourPtr_ = Tree_->father_pt()->father_pt()->son_pt(indexFatherNeighbour);
    }
}

/*!
 * \details Converts the index found in indexPath[differenceInLevel - 1] to account for refined neighbour Tree_.
 * In the end this results in only finding the neighbouring leaves of Tree_. Return the converted index, either the
 * index of the neighbouring sonPtr or -1, in case of no neighbour in that direction
 *
 * @param[in] std::vector<int> &indexPath, contain indices of the path taken to find the neighbour of Tree_
 * @param[in] int &differenceInLevel, difference in level between Tree_ and neigbourPtr_
 * @param[in] const int &dir_, direction in which the elements neighbour needs to be determined following the oomph-lib
 *              numbering convention in direction (-x = 0, +x = 1, ... , +z = 5)
 * @param[out] int convertedIndex, converted index from indexPath.
 */
template <class ELEMENT>
int UnderResolvedCoupling<ELEMENT>::convertIndexFromPath(std::vector<int>& indexPath, const int& differenceInLevel_, const int& dir_)
{
    int convertedIndex;
    int indexToConvert = indexPath[differenceInLevel_-1]; // -1 to get indices right

    std::vector<int> possibleIndices; // indices
    if (indexToConvert == 0)
    {
        possibleIndices = {1, -1, 2, -1, 4, -1};
    }
    else if (indexToConvert == 1)
    {
        possibleIndices = {-1, 0, 2, -1, 5, -1};
    }
    else if (indexToConvert == 2)
    {
        possibleIndices = {3, -1, -1, 0, 6, -1};
    }
    else if (indexToConvert == 3)
    {
        possibleIndices = {-1, 2, -1, 1, 7, -1};
    }
    else if (indexToConvert == 4)
    {
        possibleIndices = {5, -1, 6, -1, -1, 0};
    }
    else if (indexToConvert == 5)
    {
        possibleIndices = {-1, 4, 7, -1, -1, 1};
    }
    else if (indexToConvert == 6)
    {
        possibleIndices = {7, -1, -1, 4, -1, 2};
    }
    else if (indexToConvert == 7)
    {
        possibleIndices = {-1, 6, -1, 5, -1, 3};
    }

    convertedIndex = possibleIndices[dir_];

    return convertedIndex;
}

/*!
 * \details Stick all leaves of the neighbouring tree that is in contact with Tree_ in a certain direction into a vector.
 * Recursively so, in case of the neighbouring tree having sons.
 *
 * @param[in] Tree* &neighbourTree_, pointer to neighbouring Tree
 * @param[in] std::vector<Tree*> & allNeighbouringLeaves_, vector in which all neighbouring leaves will be stored
 * @param[in] const int &dir_, direction in which the elements neighbour needs to be determined following the oomph-lib
 *              numbering convention in direction (-x = 0, +x = 1, ... , +z = 5)
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::stickLeavesOfEqualOrSmallerSizedNeighbourRecursiveInDir(oomph::Tree* neighbourTree_, std::vector<oomph::Tree*> &allNeighbouringLeaves_, const int &dir_)
{
    std::vector<unsigned int> iSonsForMinX = {1, 3, 5, 7};
    std::vector<unsigned int> iSonsForMaxX = {0, 2, 4, 6};
    std::vector<unsigned int> iSonsForMinY = {2, 3, 6, 7};
    std::vector<unsigned int> iSonsForMaxY = {0, 1, 4, 5};
    std::vector<unsigned int> iSonsForMinZ = {4, 5, 6, 7};
    std::vector<unsigned int> iSonsForMaxZ = {0, 1, 2, 3};

    std::vector<unsigned int> doThisVec;
    if (dir_ == 0)
    {
        doThisVec = iSonsForMinX;
    }
    else if (dir_ == 1)
    {
        doThisVec = iSonsForMaxX;
    }
    else if (dir_ == 2)
    {
        doThisVec = iSonsForMinY;
    }
    else if (dir_ == 3)
    {
        doThisVec = iSonsForMaxY;
    }
    else if (dir_ == 4)
    {
        doThisVec = iSonsForMinZ;
    }
    else if (dir_ == 5)
    {
        doThisVec = iSonsForMaxZ;
    }

    unsigned int numsons = neighbourTree_->nsons();
    if (numsons > 0)
    {
        for(auto i : doThisVec)
        {
            stickLeavesOfEqualOrSmallerSizedNeighbourRecursiveInDir(neighbourTree_->son_pt(i), allNeighbouringLeaves_, dir_);
        }
    }
    else
    {
        allNeighbouringLeaves_.push_back(neighbourTree_);
    }
}

/*!
 * \details update listOfdVoidagedxInElem
 *
 * \param[in]
 * \param[out]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfdVoidagedxInElem()
{
    logger(INFO,"Call to updateListOfdVoidagedxInElem()");
    double nEl = mesh_pt()->nelement();
    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        updateListOfdVoidagedxInElem(iEl);
    }
}


/*!
 * \details update listOfdVoidagedxInElem for elements surrounding iEl
 *
 * \param[in]
 * \param[out]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfdVoidagedxAroundElem(const unsigned &iEl)
{
/*    // non-refinable version
    std::vector<long int> itsNeighboursTmp = getNeighbourElementNumbers(iEl);
    // You do not need to update it for the element itself, as this function is called twice. Once for the element,
    // and once for the element that the particle moved to. However, if the particle moves diagonally I [Mitchel],
    // am not sure how the code picks this up. Therefore iEl is pushed back just to be sure
    std::vector<long int> updateTheseElements = itsNeighboursTmp;
    updateTheseElements.push_back(iEl);

    for (unsigned int ithUpdate = 0; ithUpdate < updateTheseElements.size(); ithUpdate++)
    {
        long int updateThisEl = updateTheseElements[ithUpdate];
        if (updateThisEl != -999)
        {
            updateListOfdVoidagedxInElem(updateThisEl);
        }
    }    */

    //refineable version
    for (unsigned iDir =0; iDir < 6; iDir++)
    {
        std::vector<oomph::Tree*> neighbours = getNeighboursOfEl(iEl, iDir);
        for (auto iNeighbour = 0; iNeighbour < neighbours.size(); iNeighbour++)
        {
            if (neighbours[iNeighbour] != nullptr)
            {
                updateListOfdVoidagedxInElem(neighbours[iNeighbour]->object_pt()->number());
            }
        }
    }
}

/*!
 * \details update listOfdVoidagedxInElem for specific iEl
 *
 * \param[in] const unsigned iEl, fluid element id
 * \param[out]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfdVoidagedxInElem(const unsigned &iEl)
{
    /*// Working for unrefined elements by getVoidageInElemFromList() using the fact that unstored neighbours are assigned value -999
    std::vector<long int> itsNeighbours = getNeighbourElementNumbers(iEl);
    for (unsigned int iDir = 0; iDir < 3; iDir++)
    {
        long int neighbour0 = itsNeighbours[2*iDir];
        long int neighbour2 = itsNeighbours[2*iDir+1];
        logger(DEBUG,"Neighbours: % %",neighbour0,neighbour2);
        double eps0, eps1, eps2;

        // Separate cases for what the neighbour numbers are. -999 means that it has no neighbour on that side.
        // Overwriting neighbour value to make code less repeatable
        // Neighbour only in negative direction
        if (neighbour0 != -999 && neighbour2 == -999)
        {
            neighbour2 = iEl;
        }
            // Neighbour only in positive direction
        else if (neighbour0 == -999 && neighbour2 != -999)
        {
            neighbour0 = iEl;
        }
            // No neighbours on either side, default dvoidage_dx to zero
        else
        {
            listOfdVoidagedxInElem[iEl][iDir] = 0.0;
            continue; //Skip computation below
        }

        eps0 = getVoidageInElemFromList(neighbour0);
        eps1 = getVoidageInElemFromList(iEl);
        eps2 = getVoidageInElemFromList(neighbour2);

        oomph::Vector<double> s0(3,0.0);
        oomph::Vector<double> x0(3,0.0);
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(neighbour0))->interpolated_x(s0,x0);
        oomph::Vector<double> s1(3,0.0);
        oomph::Vector<double> x1(3,0.0);
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(iEl))->interpolated_x(s1,x1);
        oomph::Vector<double> s2(3,0.0);
        oomph::Vector<double> x2(3,0.0);
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(neighbour2))->interpolated_x(s2,x2);

        logger(DEBUG,"Voidage: % % %", eps0, eps1,eps2);
        logger(DEBUG,"x: % % %",x0[iDir],x1[iDir],x2[iDir]);

        // See Oct 14 Notes for derivation: 2*f'(x1) = (f(x2)-f(x1))/h2 - (f(x0)-f(x1))/h1
        // Note that is does not take depsdx into account for adapted grids, it takes centroid of element:
        // It is correct for un-adapted grids as then dh = 0;
        // __________________
        // |   |   |         |
        // | x0| x1|         |
        // |___|___|    x2   |
        // |   |   |         |
        // |   |   |         |
        // -------------------

        if ((neighbour0 != iEl) && (neighbour2 != iEl)) // Using both neighbours
        {
            listOfdVoidagedxInElem[iEl][iDir] = 0.5* ((eps2 - eps1)/(x2[iDir]-x1[iDir]) - (eps0 - eps1)/(x1[iDir] - x0[iDir]));
        }
        else if ((neighbour0 != iEl) && (neighbour2 == iEl)) // Using only left neighbours
        {
            listOfdVoidagedxInElem[iEl][iDir] = (eps1 - eps0)/(x1[iDir] - x0[iDir]);
        }
        else if ((neighbour0 == iEl) && (neighbour2 != iEl)) // Using only right neighbours
        {
            listOfdVoidagedxInElem[iEl][iDir] = (eps2 - eps1)/(x2[iDir] - x1[iDir]);
        }
        else
        {
            listOfdVoidagedxInElem[iEl][iDir] = 0.0;
        }

        logger(DEBUG,"1) listOfdVoidagedxInElem[iEl %][iDir %] = %", iEl, iDir, listOfdVoidagedxInElem[iEl][iDir]);
    }*/


    logger(DEBUG,"Element %",iEl);
    // Make it working for refineablemeshes
    for (unsigned int iDir = 0; iDir < 3; iDir++)
    {
        auto neighboursMinDir = listOfRefinableNeighbours_[iEl][2*iDir];
        auto neighboursMaxDir = listOfRefinableNeighbours_[iEl][2*iDir+1];

        int levelIEl = dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(iEl))->tree_pt()->level(); // neighboursMaxDir[0]->level();

        double level;
        std::vector<double> voidages;
        double avgVoidageLeft = 0.0;
        double voidageMiddle = getVoidageInElemFromList(iEl);
        double avgVoidageRight = 0.0;

        double avgPosLeft = 0.0;
        oomph::Vector<double> s1(3,0.0);
        oomph::Vector<double> x1(3,0.0);
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(iEl))->interpolated_x(s1,x1);
        double avgPosRight = 0.;

        int sumNonNullMin = 0;
        int sumNonNullMax = 0;

        //FIXME make better approx for depsdx

        for (unsigned iNeighbour = 0; iNeighbour < neighboursMinDir.size(); iNeighbour++)
        {
            if (neighboursMinDir[iNeighbour] != nullptr)
            {
                level = neighboursMinDir[iNeighbour]->level();
                sumNonNullMin++;

                double iWeight = std::min(1.0, 1.0 / (pow(2, level - levelIEl)));
                avgVoidageLeft += iWeight * getVoidageInElemFromList(neighboursMinDir[iNeighbour]->object_pt()->number());
                logger(DEBUG,"elLeft = %", neighboursMinDir[iNeighbour]->object_pt()->number());

                oomph::Vector<double> s0(3,0.0);
                oomph::Vector<double> x0(3,0.0);
                dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(neighboursMinDir[iNeighbour]->object_pt()->number()))->interpolated_x(s0,x0);

                avgPosLeft += iWeight * x0[iDir];
            }
        }
        for (unsigned iNeighbour = 0; iNeighbour < neighboursMaxDir.size(); iNeighbour++)
        {
            if (neighboursMaxDir[iNeighbour] != nullptr)
            {
                level = neighboursMaxDir[iNeighbour]->level();
                sumNonNullMax++;

                double iWeight = std::min(1.0, 1.0 / (pow(2, level - levelIEl)));
                avgVoidageRight += iWeight * getVoidageInElemFromList(neighboursMaxDir[iNeighbour]->object_pt()->number());
                logger(DEBUG,"elRight = %", neighboursMaxDir[iNeighbour]->object_pt()->number());

                oomph::Vector<double> s2(3,0.0);
                oomph::Vector<double> x2(3,0.0);
                dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(neighboursMaxDir[iNeighbour]->object_pt()->number()))->interpolated_x(s2,x2);

                avgPosRight += iWeight * x2[iDir];
            }
        }

        //logger(INFO,"%, %, %,\t %, %, %",
        //        avgVoidageLeft, voidageMiddle, avgVoidageRight, avgPosLeft, x1[iDir], avgPosRight);
        logger(DEBUG, "sumNonNullMin = %, sumNonNullMax = %", sumNonNullMin, sumNonNullMax);
        logger(DEBUG,"avgVoidageLeft = %, voidageMiddle = %, avgVoidageRight = %, avgPosLeft = %, posMiddle = %, avgPosRight = %",
               avgVoidageLeft, voidageMiddle, avgVoidageRight, avgPosLeft, x1[iDir], avgPosRight);
        if (sumNonNullMin > 0 && sumNonNullMax > 0) // Neighbours on both sides
        {
            listOfdVoidagedxInElem_[iEl][iDir] = 0.5* ((avgVoidageRight - voidageMiddle)/(avgPosRight-x1[iDir]) - (avgVoidageLeft - voidageMiddle)/(x1[iDir] - avgPosLeft));
        }
        else if (sumNonNullMin == 0 && sumNonNullMax > 0) // Using only right neighbours
        {
            listOfdVoidagedxInElem_[iEl][iDir] = (avgVoidageRight - voidageMiddle)/(avgPosRight-x1[iDir]);
        }
        else if (sumNonNullMin > 0 && sumNonNullMax == 0) // Using only left neighbours
        {
            listOfdVoidagedxInElem_[iEl][iDir] = (voidageMiddle - avgVoidageLeft)/(x1[iDir] - avgPosLeft);
        }
        else if (sumNonNullMin == 0 && sumNonNullMax == 0)
        {
            listOfdVoidagedxInElem_[iEl][iDir] = 0.0;
        }

        logger(DEBUG,"listOfdVoidagedxInElem[iEl %][iDir %] = %", iEl, iDir, listOfdVoidagedxInElem_[iEl][iDir]);
    }

    logger(DEBUG,"dvoidagedx = %, %, % \n", listOfdVoidagedxInElem_[iEl][0], listOfdVoidagedxInElem_[iEl][1], listOfdVoidagedxInElem_[iEl][2]);

}

/*!
 * \details update listOfdVoidagedtInElem. Hard-coded 5 history values of voidage, using one-sided 5pt
 * stencil to determine dvoidage/dt. Update historyValuesVoidage as well.
 *
 * \param[in]
 * \param[out]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfdVoidagedtInElem()
{
    logger(DEBUG,"Call to updateListOfdVoidagedtInElem()");

    double nEl = mesh_pt()->nelement();
    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        updateHistoryValuesVoidage(iEl);
        std::vector<double> histValuesVoidage = historyValuesVoidage_[iEl];

        logger(DEBUG,"histValuesVoidage = % % % % %", histValuesVoidage[0], histValuesVoidage[1], histValuesVoidage[2], histValuesVoidage[3], histValuesVoidage[4]);

        listOfdVoidagedtInElem_[iEl] = 1.0 / (12.0 * getTimeStepOomph()) * (3.0 * histValuesVoidage[0]
                                                                            - 16.0 * histValuesVoidage[1]
                                                                            + 36.0 * histValuesVoidage[2]
                                                                            - 48.0 * histValuesVoidage[3]
                                                                            + 25.0 * histValuesVoidage[4]);
        logger(DEBUG,"dVoidagedt[iEl = %]: % ", iEl, listOfdVoidagedtInElem_[iEl]);
    }
}

/*!
 * \details update listOfdVoidagedtInElem upon adapting mesh. It uses the Trees and Leaves of the previous timestep to
 * find out if it was a tree/leave or non-existent in the previous timestep (prevLeaves/prevTrees). Based on this
 * information is transferred to the adapted mesh. (hard coded for OcTree as it uses 8 sons in the loop)
 *
 * \param[in] const std::vector< std::vector<double> > &histValuesBeforeAdapt, contain history values of all elements before adapt
 * \param[out]const std::vector<double> &dVoidagedtBeforeAdapt, contains dvoidagedt of all elements before adapt
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateListOfdVoidagedtAfterAdapt(const std::vector< std::vector<double> > &histValuesBeforeAdapt, const std::vector<double> &dVoidagedtBeforeAdapt)
{

    int nEl = mesh_pt()->nelement();
    logger(DEBUG,"nEl = %",nEl);
    int atPrevEl = 0;
    int cntForRef = 0;

    oomph::Vector<oomph::Tree*> currLeaves;
    oomph::Vector<oomph::Tree*> currTrees;
    mesh_pt()->forest_pt()->stick_leaves_into_vector(currLeaves);
    mesh_pt()->forest_pt()->stick_all_tree_nodes_into_vector(currTrees);

    logger(DEBUG,"prevLeaves.size() = %, prevTrees.size() = %", prevLeaves_.size(),prevTrees_.size());
    logger(DEBUG,"currLeaves.size() = %, currTrees.size() = %", currLeaves.size(),currTrees.size());

    for (unsigned iEl = 0; iEl < nEl; iEl++)
    {
        logger(DEBUG,"Element %, atPrevEl = %", iEl, atPrevEl);

        oomph::Tree* currLeaf = currLeaves[iEl];
        bool wasTree = false;
        bool wasLeaf = false;

        for (auto iPrevLeaf : prevLeaves_)
        {
            if (iPrevLeaf == currLeaf)
            {
                wasLeaf = true;
                break;
            }
        }
        for (auto iPrevTree : prevTrees_)
        {
            if (iPrevTree == currLeaf)
            {
                wasTree = true;
                break;
            }
        }

        if (wasLeaf && wasTree)
        {
            logger(DEBUG,"Element % wasLeaf and wasTree",iEl);

            historyValuesVoidage_[iEl] = histValuesBeforeAdapt[atPrevEl];
            listOfdVoidagedtInElem_[iEl] = dVoidagedtBeforeAdapt[atPrevEl];
            atPrevEl+=1; // Increment by 1 as not (un)refined
        }
        else if (!wasLeaf && !wasTree)
        {
            logger(DEBUG, "Element % is created due to refinement", iEl);

            historyValuesVoidage_[iEl] = histValuesBeforeAdapt[atPrevEl];
            listOfdVoidagedtInElem_[iEl] = dVoidagedtBeforeAdapt[atPrevEl];

            cntForRef++; // It will find its 7 other sons immediatley after iEl;
            if (cntForRef % 8 == 0)
            {
                cntForRef = 0;
                atPrevEl += 1;
            }
        }
        else if (!wasLeaf && wasTree)
        {
            logger(DEBUG,"Element % is combined with 8 sons of prevTime",iEl);
            std::vector<double> histOfSon;
            int sizeOfHist = histValuesBeforeAdapt[atPrevEl].size();
            histOfSon.resize(sizeOfHist);
            double dEpsdt = 0;

            for (int iSon = 0; iSon < 8; iSon++)
            {
                int sonNr = atPrevEl+iSon;
                logger(DEBUG,"sonNr %", sonNr);

                for (int i =0; i < sizeOfHist; i++)
                {
                    histOfSon[i] += (histValuesBeforeAdapt[sonNr][i] / 8);
                }
                dEpsdt += dVoidagedtBeforeAdapt[sonNr] / 8;
            }

            historyValuesVoidage_[iEl] = histOfSon;
            listOfdVoidagedtInElem_[iEl] = dEpsdt;

            atPrevEl += 8;
        }
        else
        {
            logger(WARN,"You should never see this, see updateListOfdVoidagedtAfterAdapt as there is an error in your refinement");
        }
    }
}

/*!
 * \details Update history values for specific element
 *
 * \param[in] const unsigned int IEl_, element id
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::updateHistoryValuesVoidage(const unsigned int &iEl_)
{
    std::vector<double> histValuesVoidage = historyValuesVoidage_[iEl_];

    // Overwrite old data
    for (unsigned int i = 0; i < histValuesVoidage.size() - 1; i++)
    {
        histValuesVoidage[i] = histValuesVoidage[i+1];
    }

    // Add current voidage as last entry
    histValuesVoidage[histValuesVoidage.size()-1] = getVoidageInElemFromList(iEl_);
}

// OUTPUT ALL ELEMENTS WITH ITS NEIGHBOURING ELEMENTS
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::outputNeighboursOfElements()
{
    unsigned long int nEl = mesh_pt()->nelement();
    for (unsigned int iEl = 0; iEl < nEl; iEl++)
    {
        logger(DEBUG, "iEL: %", iEl);
        for (unsigned int iDir = 0; iDir < 6; iDir++)
        {
            logger(DEBUG, "\tiDir: %", iDir);
            logger(DEBUG, "\t\tNeighbours:");

            for (unsigned int iNeighbour = 0; iNeighbour < listOfRefinableNeighbours_[iEl][iDir].size(); iNeighbour++)
            {
                if (listOfRefinableNeighbours_[iEl][iDir][iNeighbour] != nullptr)
                {
                    logger(DEBUG, "\t\t\t %", listOfRefinableNeighbours_[iEl][iDir][iNeighbour]->object_pt()->number());
                }
            }
        }
    }
}

/*!
 * \details get s for specific particle. Also checks if particle has moved element.
 *
 * \param[in] oomph::Vector<double> s_, const unsigned iPart_, const int pInEl_
 * \param[in]
 */
template <class ELEMENT>
void UnderResolvedCoupling<ELEMENT>::getS(oomph::Vector<double> &s_, const unsigned &iPart_, const int &pInEl_)
{
    BaseParticle* particlePtr = particleHandler.getObject(iPart_);

    oomph::Vector<double> pPos = convertVecFuncs::convertToOomphVec(particlePtr->getPosition());
    oomph::GeomObject* geom_obj_pt = nullptr;

    logger(DEBUG,"Check if pPos is in pIsInEl");
    dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl_))->locate_zeta(pPos, geom_obj_pt, s_);
    if (geom_obj_pt == nullptr)
    {
        logger(DEBUG,"particle % moved through element boundary",iPart_);
        logger(DEBUG,"s = % % %",s_[0], s_[1], s_[2]);
        updateListsForSpecificParticle(iPart_);

        int pInEl_redo = getPInElemFromList(iPart_);
        oomph::Vector<double> sNew(3,0.0);
        dynamic_cast<ELEMENT *>(mesh_pt()->element_pt(pInEl_redo))->locate_zeta(pPos, geom_obj_pt, sNew);

        if (geom_obj_pt == nullptr)
        {
            // logger(DEBUG,"pPos[iPart = %] = (%, %, %)",iPart_, pPos[0], pPos[1], pPos[2]);
            // logger(DEBUG,"s[iPart = %] = (%, %, %)",iPart_, s[0], s[1], s[2]);
            // updateListsForSpecificParticle(iPart_);

            logger(WARN,"Problem, after update still particle is outside of all elements");
            return; // FIXME For now just simply adding no forces to the particle
        }
        else
        {
            s_ = sNew;
        }
    }
}

/*!
 * Function that returns the amplification factor of a parabolic function based on time and stroke
 *
 * @param[in] time_
 * @return double ampFactor
 */
/*template <class ELEMENT>
double UnderResolvedCoupling<ELEMENT>::getInflowVel(const double &time_)
{
    double stroke = 0.14;
    double inflowVel = 0.;
    double timeSplit = 0.15;
    double timePeriod = 0.5;

    double getVal = timeSplit* myHeaviSide(fmod(time_, timePeriod),timeSplit) +
                    (timePeriod-timeSplit)* (1-myHeaviSide(fmod(time_, 0.5),timeSplit));
    if (getVal == timeSplit)
    {
        inflowVel = constants::pi/(2.* getVal) * stroke * sin(1./getVal * constants::pi * fmod(time_,timePeriod));
    }
    else if (getVal == (timePeriod-timeSplit))
    {
        inflowVel = - constants::pi/(2.* getVal) * stroke * sin(1./getVal * constants::pi * (fmod(time_,timePeriod) - timeSplit));
    }

    return inflowVel;
}*/


/*!
 * Heaviside function that returns true or false based on entry values a_ and b_.
 * Returns 0 if a_ > b_, 1 if b_ >= a_
 *
 * @param[in] a_
 * @param[in] b_
 * @return bool
 */
template <class ELEMENT>
bool UnderResolvedCoupling<ELEMENT>::myHeaviSide(const double &a_, const double &b_)
{
    return b_ >= a_;
}

namespace getDataFromElement
{
//    double getVoidageOfElement_byEl(const int& elNr_)
//    {
//        return ptrToCoupledClass->getVoidageInElemFromList(elNr_);
//    }
//    double getdVoidagedxOfElement_byEl(const int& elNr_, const int& dir_)
//    {
//        return ptrToCoupledClass->getdVoidagedxInElemFromList(elNr_,dir_);
//    }
//    double getdVoidagedtOfElement_byEl(const int& elNr_)
//    {
//        return 0.0;
//        //return ptrToCoupledClass->getdVoidagedtInElemFromList(elNr_);
//    }
//    void getBodyForceByCoupling_byEl(const int& elNr_, oomph::Vector<double> &force_)
//    {
//        ptrToCoupledClass->getBodyForceInElemByCoupling(elNr_,force_);
//    }

    std::clock_t totaltime = 0.0;
    int long cnt = 0;
    void printTotalTime()
    {
        std::cout << "TotalTime in getVoidageOfElement = " << 1000.*totaltime / CLOCKS_PER_SEC << std::endl;
        std::cout << "count getVoidageOfElement calls = " << cnt << std::endl;
    }
    double getVoidageOfElement(const double& time, const oomph::Vector<double>& x)
    {
        std::clock_t start = clock();
        cnt++;

//        double nEl = ptrToCoupledClass->mesh_pt()->nelement();
//
//        for (unsigned iEl = 0; iEl < nEl; iEl++)
//        {
//            oomph::Vector<double> s(3, 0.0);
//            oomph::GeomObject* geom_obj_pt = nullptr;
//            dynamic_cast<ELEMENT *>(ptrToCoupledClass->mesh_pt()->element_pt(
//                    iEl))->locate_zeta(x, geom_obj_pt, s);
//            if (geom_obj_pt != nullptr)
//            {
//                std::clock_t end = clock();
//                totaltime += (end-start);
//
//                return ptrToCoupledClass->listOfVoidageInElem_[iEl];
//            }
//        }
        return 1.0;
    }
}


#endif //MERCURYDPM_UNDERRESOLVEDCOUPLING_H
