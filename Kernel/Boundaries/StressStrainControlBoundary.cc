//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provid->d that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provid->d with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "StressStrainControlBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"
#include "MpiDataClass.h"
#include "MpiContainer.h"
#include "DPMBase.h"
#include "PeriodicBoundary.h"
#include "LeesEdwardsBoundary.h"

/*!
 * \details constructor, set all the parameters to zero.
 */
StressStrainControlBoundary::StressStrainControlBoundary()
        : BaseBoundary()
{
    stressGoal_.setZero();
    strainRate_.setZero();
    gainFactor_.setZero();
    isStrainRateControlled_ = true;
    stressTotal_.setZero();
    dstrainRate_.setZero();
    lengthBox_.setZero();
    centerBox_.setZero();
    relativeToCenter_.setZero();
    integratedShift_ = 0.0;
    //
    logger(DEBUG, "StressStrainControlBoundary::StressStrainControlBoundary() finished");
}

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 */
StressStrainControlBoundary* StressStrainControlBoundary::copy() const
{
    return new StressStrainControlBoundary(*this);
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void StressStrainControlBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " stressGoal " << stressGoal_;
    os << " strainRate " << strainRate_;
    os << " gainFactor " << gainFactor_;
    os << " isStrainRateControlled " << isStrainRateControlled_;
    os << " integratedShift " << integratedShift_;
//    os << leesEdwardsBoundaries_.size() << ' ';
//    for (const LeesEdwardsBoundary& b : leesEdwardsBoundaries_)
//    {
//        os << b << ' ';
//    }
//
//    os << periodicBoundaries_.size() << ' ';
//    for (const PeriodicBoundary& b : periodicBoundaries_)
//    {
//        os << b << ' ';
//    }
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in] is      the istream
 */
void StressStrainControlBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> stressGoal_;
    is >> dummy >> strainRate_;
    is >> dummy >> gainFactor_;
    is >> dummy >> isStrainRateControlled_;
    set(stressGoal_, strainRate_, gainFactor_, isStrainRateControlled_);
    is >> dummy >> integratedShift_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object class' name
 */
std::string StressStrainControlBoundary::getName() const
{
    return "StressStrainControlBoundary";
}

/*!
 * \details This is where the stress-strain control is implemented and
 * the boundary will be checked each timestep to make sure the target
 * stress/strain are achieved.
 */
void StressStrainControlBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& particleHandler)
{
    checkPeriodicLeesEdwardsBoundariesAfterParticlesMove(particleHandler);
    
    //Real Stress and StrainRate Control every time step
    logger.assert(getHandler() != nullptr,
                  "you need to set the handler of this boundary before the parameters can be set");
    
    determineLengthAndCentre();
    
    //this activate only the stress control
    if (stressGoal_.XX != 0 || stressGoal_.YY != 0 || stressGoal_.ZZ != 0 ||
        stressGoal_.XY != 0)
    {
        computeStrainRate();
    }
    activateStrainRateControl(particleHandler);
}

/*!
 * \details This function checks the boundaries after the particles are moved,
 * to serve the stress-control for the following steps.
 */
void StressStrainControlBoundary::checkPeriodicLeesEdwardsBoundariesAfterParticlesMove(ParticleHandler& particleHandler)
{
    // Call Boundaries
    for (PeriodicBoundary& b : periodicBoundaries_)
    {
        b.checkBoundaryAfterParticlesMove(particleHandler);
    }
    for (LeesEdwardsBoundary& b : leesEdwardsBoundaries_)
    {
        b.checkBoundaryAfterParticlesMove(particleHandler);
    }
}

/*!
 * \details This function determines both the length in x,y,z direction
 * and center location of the domain.
 */
void StressStrainControlBoundary::determineLengthAndCentre()
{
    const DPMBase* const dpm = getHandler()->getDPMBase();
    lengthBox_.X = (dpm->getXMax() - dpm->getXMin());
    lengthBox_.Y = (dpm->getYMax() - dpm->getYMin());
    lengthBox_.Z = (dpm->getZMax() - dpm->getZMin());
    // Box length Lx, Lyand Lz, and center point C of the box
    centerBox_.X = (dpm->getXMax() + dpm->getXMin()) / 2.0;
    centerBox_.Y = (dpm->getYMax() + dpm->getYMin()) / 2.0;
    centerBox_.Z = (dpm->getZMax() + dpm->getZMin()) / 2.0;
    
}

/*!
 * \details This function is used to compute the new strainrate tensor based on
 * the stress differences between the actual and user set values.
 */
void StressStrainControlBoundary::computeStrainRate()
{
    //Set strainrate change tensor to zero
    dstrainRate_.setZero();
    
    //Get the timestep dt
    const Mdouble timeStep = getHandler()->getDPMBase()->getTimeStep();
    
    Matrix3D dstrainRate_;
    // calculate the stress total and average over the volume
    Matrix3D stressTotal_ = getHandler()->getDPMBase()->getTotalStress();
    
    // amount by which the strainrate tensor has to be changed
    if (stressGoal_.XX != 0)
    {
        dstrainRate_.XX = dstrainRate_.XX + gainFactor_.XX * timeStep * (stressTotal_.XX - stressGoal_.XX);
        logger(VERBOSE, "StressXX = %",stressTotal_.XX);
    }
    if (stressGoal_.YY != 0)
    {
        dstrainRate_.YY = dstrainRate_.YY + gainFactor_.YY * timeStep * (stressTotal_.YY - stressGoal_.YY);
        logger(VERBOSE, "StressYY = %",stressTotal_.YY);
    }
    if (stressGoal_.ZZ != 0)
    {
        dstrainRate_.ZZ = dstrainRate_.ZZ + gainFactor_.ZZ * timeStep * (stressTotal_.ZZ - stressGoal_.ZZ);
        logger(VERBOSE, "StressZZ = %",stressTotal_.ZZ);
    }
    if (stressGoal_.XY != 0)
    {
        dstrainRate_.XY = dstrainRate_.XY + gainFactor_.XY * timeStep * (stressTotal_.XY - stressGoal_.XY);
        logger(VERBOSE, "dstrainRate.XY = %",dstrainRate_.XY);
        logger(VERBOSE, "StressXY = %",stressTotal_.XY);
        logger(VERBOSE, "StressYX = %",stressTotal_.YX);
    }
    
    //  Update the strainrate tensor
    strainRate_ = strainRate_ + dstrainRate_;
}

/*!
 * \details This function activate the strainrate control based on user inputs, it takes only
 * the strainRate_ tensor and move particles based on this tensor, also move boudnaries based
 * on the new domain size.
 */
void StressStrainControlBoundary::activateStrainRateControl(const ParticleHandler& particleHandler)
{
    const DPMBase* const dpm = getHandler()->getDPMBase();
    // update the domain size based on the new strainrate tensor
    updateDomainSize();
    
    //  Move the boundaries to next time step
    if (strainRate_.XY != 0 || stressGoal_.XY != 0)
    {
        determineStressControlledShearBoundaries();
    }
    else
    {
        //  Move the boundary in x direction to next time step
        periodicBoundaries_[0].set(Vec3D(1.0, 0.0, 0.0), dpm->getXMin(), dpm->getXMax());
        //  Move the boundary in y direction to next time step
        periodicBoundaries_[1].set(Vec3D(0.0, 1.0, 0.0), dpm->getYMin(), dpm->getYMax());
        //  Move the boundary in z direction to next time step
        periodicBoundaries_[2].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
    }
    
    //  Give the strain-rate for all particles and move them to next timestep before integration
    if (isStrainRateControlled_)
    {
        for (auto& p : particleHandler)
        {
            relativeToCenter_.X = p->getPosition().X - centerBox_.X;
            relativeToCenter_.Y = p->getPosition().Y - centerBox_.Y;
            relativeToCenter_.Z = p->getPosition().Z - centerBox_.Z;
            p->move(Vec3D(strainRate_.XX * dpm->getTimeStep() * relativeToCenter_.X +
                          strainRate_.XY * dpm->getTimeStep() * relativeToCenter_.Y,
                          strainRate_.YY * dpm->getTimeStep() * relativeToCenter_.Y,
                          strainRate_.ZZ * dpm->getTimeStep() * relativeToCenter_.Z));
        }
    }
}

/*!
 * \details This function is used to upate the domain size based on the strainRate tensor,
 * note that the system is symmetric and thereofore we have to update boundary in both Min and Max.
 */
void StressStrainControlBoundary::updateDomainSize()
{
    DPMBase* const dpm = getHandler()->getDPMBase();
    //  Change the system size according to next time step
    dpm->setXMax(dpm->getXMax() + 0.5 * lengthBox_.X * strainRate_.XX * dpm->getTimeStep());
    dpm->setXMin(dpm->getXMin() - 0.5 * lengthBox_.X * strainRate_.XX * dpm->getTimeStep());
    dpm->setYMax(dpm->getYMax() + 0.5 * lengthBox_.Y * strainRate_.YY * dpm->getTimeStep());
    dpm->setYMin(dpm->getYMin() - 0.5 * lengthBox_.Y * strainRate_.YY * dpm->getTimeStep());
    dpm->setZMax(dpm->getZMax() + 0.5 * lengthBox_.Z * strainRate_.ZZ * dpm->getTimeStep());
    dpm->setZMin(dpm->getZMin() - 0.5 * lengthBox_.Z * strainRate_.ZZ * dpm->getTimeStep());
    
    // Box length Lx, Lyand Lz for next time step
    lengthBox_.X = (dpm->getXMax() - dpm->getXMin());
    lengthBox_.Y = (dpm->getYMax() - dpm->getYMin());
    lengthBox_.Z = (dpm->getZMax() - dpm->getZMin());
}

/*!
 * \details This function determines the boundary for stress-controlled shear situation.
 * In this case, the Lees-Edwards boundary in x-y directions and Normal periodic boundary
 * in z direction are combined together as a cuboid shear box.
 */
void StressStrainControlBoundary::determineStressControlledShearBoundaries()
{
    const DPMBase* const dpm = getHandler()->getDPMBase();
    //Determine the current shear velocity and shift of the Lees-Edwards boundary
    const Mdouble velocity_xy = strainRate_.XY * lengthBox_.Y;
    integratedShift_ += velocity_xy * dpm->getTimeStep();
    
    // Determine how to move boundaries based on strainrate control or boundary control
    if (isStrainRateControlled_)
    {
        //  Move the Lees-Edwards boundary in z direction to next time step
        leesEdwardsBoundaries_[0].set(
                [this](Mdouble time)
                { return integratedShift_; },
                [](Mdouble time UNUSED)
                { return 0; },
                dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
    }
    else
    {
        //  Update the velocity of Lees-Edwards boundary in x-y direction to next time step
        leesEdwardsBoundaries_[0].set(
                [this](Mdouble time)
                { return integratedShift_; },
                [velocity_xy](Mdouble time UNUSED)
                { return velocity_xy; },
                dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
    }
    //  Move the boundary in z direction to next time step
    periodicBoundaries_[0].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
}

/*!
 * \details This function sets the inputs for the whole StressStrainControlBoundary
 * based on the user inputs.
 */
void
StressStrainControlBoundary::set(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
                                 bool isStrainRateControlled)
{
    periodicBoundaries_.clear();
    leesEdwardsBoundaries_.clear();

    isStrainRateControlled_ = isStrainRateControlled;
    stressGoal_ = stressGoal;
    strainRate_ = strainRate;
    gainFactor_ = gainFactor;
    
    logger.assert_always(getHandler() != nullptr,
                         "you need to set the handler of this boundary before the parameters can be set");
    const DPMBase* dpm = getHandler()->getDPMBase();
    
    // add sensible checks if parameters are valid
    
    if (stressGoal.XY == 0 && strainRate.XY == 0)
    {
        logger(INFO, "Shear rate is zero, setting up three periodic boundaries");
        // Set up new box periodic boundaries, no simple shear activated
        PeriodicBoundary boundary;
        boundary.setHandler(getHandler());
        boundary.set(Vec3D(1.0, 0.0, 0.0), dpm->getXMin(), dpm->getXMax());
        periodicBoundaries_.push_back(boundary);
        boundary.set(Vec3D(0.0, 1.0, 0.0), dpm->getYMin(), dpm->getYMax());
        periodicBoundaries_.push_back(boundary);
        boundary.set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        periodicBoundaries_.push_back(boundary);
    }
    else if (stressGoal.XY != 0 || strainRate.XY != 0)
    {
        logger(INFO, "Shear rate is not zero, setting up Lees-Edwards in xy and a periodic boundary in z");
        
        PeriodicBoundary boundary;
        boundary.setHandler(getHandler());
        boundary.set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        periodicBoundaries_.push_back(boundary);
        // Lees Edwards bc in y direction & periodic boundary in x direction, simple shear boundary activated
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.setHandler(getHandler());
        leesEdwardsBoundary.set(
                [this](Mdouble time UNUSED)
                { return integratedShift_; },
                [](Mdouble time UNUSED)
                { return 0; },
                dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
        leesEdwardsBoundaries_.push_back(leesEdwardsBoundary);
    }
    else
    {
        logger(ERROR, "This should not happen; please check implementation of StressStrainControlBoundary::set");
    }
    
}

/*!
 * \details This function is used to create ghost particles such that
 * the stress are calculated correctly after restart from another configuration.
 */
void StressStrainControlBoundary::createPeriodicParticles(ParticleHandler& particleHandler)
{
    // Call Boundaries
    for (PeriodicBoundary& b : periodicBoundaries_)
    {
        b.createPeriodicParticles(particleHandler);
    }
    for (LeesEdwardsBoundary& b : leesEdwardsBoundaries_)
    {
        b.createPeriodicParticles(particleHandler);
    }
}