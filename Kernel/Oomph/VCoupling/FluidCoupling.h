//
// Created by mitchel on 11/7/22.
//

#ifndef MERCURYDPM_FLUIDCOUPLING_H
#define MERCURYDPM_FLUIDCOUPLING_H

#include "VCoupling.h"
#include "FluidCoupledElement.h"

using namespace oomph;

template<class M, class O>
class FluidCoupling : public VolumeCoupling<M, O>
{
public:
    
    typedef typename O::ELEMENT ELEMENT;
    
    void solveParticleFluidCoupling()
    {
        // compute nStep
        unsigned nStep = O::getOomphTimeStep()/M::getTimeStep();
        if (nStep==0) {
            // if oomph has a smaller time step than Mercury decrease Mercury timestep
            logger(INFO, "Set nStep %, change mercuryTimeStep from % to %", nStep, M::getTimeStep(), O::getOomphTimeStep());
            nStep = 1;
            M::setTimeStep(O::getOomphTimeStep());
        } else {
            logger(INFO, "Set nStep %, change oomphTimeStep from % to %", nStep, O::getOomphTimeStep(), nStep * M::getTimeStep());
            O::setOomphTimeStep(nStep * M::getTimeStep());
        }
    
        // call solve routine
        solveParticleFluidCoupling(nStep);
    }
    
    // solve volume coupled OomphMercuryProblem
    void solveParticleFluidCoupling(unsigned nStep)
    {
        // check whether time steps are set
        logger.assert_always(O::getOomphTimeStep()>0,"Oomph time step not initialised");
        logger.assert_always(M::getTimeStep()>0,"Mercury time step not initialised");
    
        unsigned nDone = 0; //< last written file number
        while (M::getTime() < M::getTimeMax())
        {
            this->actionsBeforeOomphTimeStep();
            // solve the coupled problem for one time step
            logger(INFO,"computeOneTimeStepForFluidCoupling");
            computeOneTimeStepForFluidCoupling(nStep);
            /*
            // write outputs of the oomphProb; this is slaved to the vtk output of Mercury, i.e. an oomph-lib output get written everytime a Mercury vtk file gets written
            if (M::getParticlesWriteVTK() && M::getVtkWriter()->getFileCounter() > nDone)
            {
                O::doc_paraview(doc_info);
                //O::writeToVTK(); //FIXME No vtkwriter for fluid yet
                nDone = M::getVtkWriter()->getFileCounter();
            }
            */
        }
    }
    
    void computeOneTimeStepForFluidCoupling(unsigned nStepsMercury)
    {
        logger(INFO,"solveMercury");
        BaseCoupling<M,O>::solveMercury(nStepsMercury);
        logger(INFO,"solveOomph");
        BaseCoupling<M,O>::solveOomph();
    }
    
    // Set and get functions for fractions of interaction coupling forces
    void setFractionDragForce(double df_) {fractionDragForce = df_;}
    void setFractionShearForce(double sf_) {fractionShearForce = sf_;}
    void setFractionSaffmanLift(double sl_) {fractionSaffmanLift = sl_;}
    void setFractionMagnusLift(double ml_) {fractionMagnusLift = ml_;}
    void setFractionPressureGradient(double pg_) {fractionPressureGradient = pg_;}
    void activateAllInteractionForces()
    {
        fractionDragForce = 1.0;
        fractionShearForce = 1.0;
        fractionSaffmanLift = 1.0;
        fractionMagnusLift = 1.0;
        fractionPressureGradient = 1.0;
    }
    
    double getFractionDragForce() {return fractionDragForce;}
    double getFractionShearForce() {return fractionShearForce;}
    double getFractionSaffmanLift() {return fractionSaffmanLift;}
    double getFractionMagnusLift() {return fractionMagnusLift;}
    double getFractionPressureGradient() {return fractionPressureGradient;}
    

    
    void updateFluidVolumeFractions()
    {
        for (unsigned i = 0; i < this->fluid_mesh_pt()->nelement(); i++)
        {
            ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->fluid_mesh_pt()->element_pt(i));
            updateFluidVolumeFractionForElement(el_pt);
        }
    }
    
    void updateFluidVolumeFractionForElement(ELEMENT* el_pt)
    {
        double VolumeElement;
        double VolumeParticles = 0.0;
        oomph::Vector<double> posNode0(3, 0.0);
        oomph::Vector<double> posNodeEnd(3, 0.0);
        
        const unsigned int iNode = el_pt->nnode();
        
        // Note that index of end node is iNode - 1, else segmentation fault as node_pt(iNode) does not exist
        el_pt->node_pt(0)->position(posNode0);
        el_pt->node_pt(iNode - 1)->position(posNodeEnd);
        
        VolumeElement = abs(posNode0[0] - posNodeEnd[0]) * abs(posNode0[1] - posNodeEnd[1]) * abs(posNode0[2] - posNodeEnd[2]);
        
        for (unsigned int iP = 0; iP < this->particleHandler.getNumberOfObjects(); iP++)
        {
            Vec3D pPos = this->particleHandler.getObject(iP)->getPosition();
            
            oomph::Vector<double> s(3, 0.0);
            GeomObject* geom_obj_pt = nullptr;
            el_pt->locate_zeta(convertVecFuncs::convertToOomphVec(pPos), geom_obj_pt, s);
            
            if (geom_obj_pt != nullptr)
            {
                VolumeParticles += this->particleHandler.getObject(iP)->getVolume();
            }
        }
        
        el_pt->setFluidVolumeFraction(1. - VolumeParticles / VolumeElement);
    }
    
private:
    // Fractions of interaction coupling forces, by default all fractions are 0.0;
    double fractionDragForce = 0.0;
    double fractionShearForce = 0.0;
    double fractionSaffmanLift = 0.0;
    double fractionMagnusLift = 0.0;
    double fractionPressureGradient = 0.0;
};

#endif //MERCURYDPM_FLUIDCOUPLING_H
