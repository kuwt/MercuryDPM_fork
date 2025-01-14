                                                                                                                        //Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
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
#include "MpiContainer.h"
#include "DPMBase.h"
#include "PeriodicBoundary.h"
#include "LeesEdwardsBoundary.h"


/*!
 * \details constructor, set all the parameters that boundary needs as inputs to zero.
 */
StressStrainControlBoundary::StressStrainControlBoundary()
        : BaseBoundary()
{
    stressGoal_.setZero();
    strainRate_.setZero();
    pGainFactor_.setZero();
    iGainFactor_.setZero();
    isStrainRateControlled_ = true;
    integratedShift_ = 0.0;
    //
    logger(DEBUG, "StressStrainControlBoundary::StressStrainControlBoundary() finished");
}
Mdouble StressStrainControlBoundary::computeStressError(){
    Matrix3D stress = getHandler()->getDPMBase()->getTotalStress();
    return (stress.XX- stressGoal_.XX);
}
/*!
 * \details Copy method; creates a copy on the boundary and returns its pointer.
 */
StressStrainControlBoundary* StressStrainControlBoundary::copy() const
{
    return new StressStrainControlBoundary(*this);
}

/*!
 * \details Writes boundary's properties to an ostream.
 * \param[in] os    the ostream.
 */
void StressStrainControlBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " stressGoal " << stressGoal_;
    os << " strainRate " << strainRate_;
    os << " pGainFactor " << pGainFactor_;
    os << " iGainFactor " << iGainFactor_;
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
 * \details Reads the boundary properties from an istream.
 * \param[in] is      the istream.
 */
void StressStrainControlBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> stressGoal_;
    is >> dummy >> strainRate_;
    is >> dummy >> pGainFactor_;
    is >> dummy >> iGainFactor_;
    is >> dummy >> isStrainRateControlled_;
    set(stressGoal_, strainRate_, pGainFactor_, isStrainRateControlled_, iGainFactor_);
    is >> dummy >> integratedShift_;
}

/*!
 * \details Returns the name of the object class.
 * \return      the object class' name.
 */
std::string StressStrainControlBoundary::getName() const
{
    return "StressStrainControlBoundary";
}

/*!
 * \details This is where the stress-strain control is implemented and
 * the boundary will be checked at each time step to make sure the target
 * stress/strain rate are achieved.
 */
void StressStrainControlBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& particleHandler)
{
    checkPeriodicLeesEdwardsBoundariesAfterParticlesMove(particleHandler);
    
    //Real Stress and StrainRate Control every time step
    logger.assert_debug(getHandler() != nullptr,
                  "you need to set the handler of this boundary before the parameters can be set");
    
    determineLengthAndCentre();
    
    //this checks if stressGoal matrix is non-zero and then activate only the stress control
    if (stressGoal_.XX != 0 || stressGoal_.YY != 0 || stressGoal_.ZZ != 0 ||
        stressGoal_.XY != 0)
    {
        //the strain rate in the corresponding direction will be calculated
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
//    const DPMBase* const dpm = getHandler()->getDPMBase();
}

/*!
 * \details This function is used to compute the new strain rate tensor based on
 * the stress differences between the actual and user set target values.
 * This is based on a PI-controller,
 *    strainRate = pGain * (stress-stressGoal) + iGain * integral(stress-stressGoal)
 */
void StressStrainControlBoundary::computeStrainRate()
{
    // calculate the stress total and average over the volume
    Matrix3D stress = getHandler()->getDPMBase()->getTotalStress();

    // controller, sets the amount by which the strain rate tensor has to be changed
    double timeStep = getHandler()->getDPMBase()->getTimeStep();
    // amount by which the strain rate tensor has to be changed
    if (stressGoal_.XX != 0) {
        strainRate_.XX = xx.apply(stress.XX - stressGoal_.XX,timeStep);
    }
    if (stressGoal_.YY != 0) {
        strainRate_.YY = yy.apply(stress.YY - stressGoal_.YY,timeStep);
    }
    if (stressGoal_.ZZ != 0) {
        strainRate_.ZZ = zz.apply(stress.ZZ - stressGoal_.ZZ,timeStep);
    }
    if (stressGoal_.XY != 0) {
        strainRate_.XY = xy.apply(stress.XY - stressGoal_.XY,timeStep);
    }
}

/*!
 * \details This function activate the strain rate control based on user inputs, it takes only
 * the strainRate_ tensor and move particles based on this tensor, also move boundaries based
 * on the new domain size.
 */
void StressStrainControlBoundary::activateStrainRateControl(const ParticleHandler& particleHandler)
{
    const DPMBase* const dpm = getHandler()->getDPMBase();
    // update the domain size based on the new strain rate tensor
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
    
    //  Give the strain-rate for all particles and move them to next time step before integration
    if (isStrainRateControlled_)
    {
        for (auto& p : particleHandler)
        {
            // Box length Lx, Ly and Lz, and center point C of the box
            Vec3D centerBox;
            centerBox.X = (dpm->getXMax() + dpm->getXMin()) / 2.0;
            centerBox.Y = (dpm->getYMax() + dpm->getYMin()) / 2.0;
            centerBox.Z = (dpm->getZMax() + dpm->getZMin()) / 2.0;
            Vec3D relativeToCenter;
            relativeToCenter.X = p->getPosition().X - centerBox.X;
            relativeToCenter.Y = p->getPosition().Y - centerBox.Y;
            relativeToCenter.Z = p->getPosition().Z - centerBox.Z;
            p->move(Vec3D(strainRate_.XX * dpm->getTimeStep() * relativeToCenter.X +
                          strainRate_.XY * dpm->getTimeStep() * relativeToCenter.Y,
                          strainRate_.YY * dpm->getTimeStep() * relativeToCenter.Y,
                          strainRate_.ZZ * dpm->getTimeStep() * relativeToCenter.Z));
        }
    }
}

/*!
 * \details This function is used to update the domain size based on the strainRate tensor,
 * note that the system is symmetric and therefore we have to update boundary in both Min and Max.
 */
void StressStrainControlBoundary::updateDomainSize()
{
    DPMBase* const dpm = getHandler()->getDPMBase();
    // Box length Lx, Ly and Lz for next time step
    Vec3D lengthBox;
    lengthBox.X = (dpm->getXMax() - dpm->getXMin());
    lengthBox.Y = (dpm->getYMax() - dpm->getYMin());
    lengthBox.Z = (dpm->getZMax() - dpm->getZMin());
    //  Change the system size according to next time step
    dpm->setXMax(dpm->getXMax() + 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
    dpm->setXMin(dpm->getXMin() - 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
    dpm->setYMax(dpm->getYMax() + 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
    dpm->setYMin(dpm->getYMin() - 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
    dpm->setZMax(dpm->getZMax() + 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
    dpm->setZMin(dpm->getZMin() - 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
}

/*!
 * \details This function determines the boundary for stress-controlled shear situation.
 * In this case, the Lees-Edwards boundary in x-y directions and normal periodic boundary
 * in z direction are combined together as a cuboid shear box.
 */
void StressStrainControlBoundary::determineStressControlledShearBoundaries()
{
    const DPMBase* const dpm = getHandler()->getDPMBase();
    //Determine the current shear velocity and shift of the Lees-Edwards boundary
    const Mdouble velocity_xy = strainRate_.XY * (dpm->getXMax() - dpm->getXMin());
    integratedShift_ += velocity_xy * dpm->getTimeStep();
    
    // Determine how to move boundaries based on strain rate control or boundary control
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
 * \param[in] stressGoal        The target stress tensor that needs to achieve at the end of the deformation.
 * \param[in] strainRate        The target strain rate tensor, cannot be set non-zero with target Stress,
 *                              i.e. stressGoal.XX != 0, then strainRate.XX = 0.
 * \param[in] pGain        The incremental factor for the stress control, usually a very small value ~ 0.0001.
 * \param[in] isStrainRateControlled        The boolean key to determine whether particles are moved by
 *                                          strain rate affine movement (true) or dragged by boundary itself (false).
 * \param[in] iGain        The incremental factor for the integral stress control
 */
void StressStrainControlBoundary::set(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& pGain, bool isStrainRateControlled, const Matrix3D& iGain)
{
    periodicBoundaries_.clear();
    leesEdwardsBoundaries_.clear();

    isStrainRateControlled_ = isStrainRateControlled;
    stressGoal_ = stressGoal;
    strainRate_ = strainRate;
    pGainFactor_ = pGain;
    iGainFactor_ = iGain;

    logger.assert_always(stressGoal.XZ == 0, "Shear stress in XZ cannot be controlled; use shear stress in XY instead");
    logger.assert_always(stressGoal.ZX == 0, "Shear stress in ZX cannot be controlled; use shear stress in XY instead");
    logger.assert_always(stressGoal.YZ == 0, "Shear stress in YZ cannot be controlled; use shear stress in XY instead");
    logger.assert_always(stressGoal.ZY == 0, "Shear stress in ZY cannot be controlled; use shear stress in XY instead");
    logger.assert_always(strainRate.XZ == 0, "Strain rate in XZ cannot be controlled; use strain rate in XY instead");
    logger.assert_always(strainRate.ZX == 0, "Strain rate in ZX cannot be controlled; use strain rate in XY instead");
    logger.assert_always(strainRate.YZ == 0, "Strain rate in YZ cannot be controlled; use strain rate in XY instead");
    logger.assert_always(strainRate.ZY == 0, "Strain rate in ZY cannot be controlled; use strain rate in XY instead");
    logger.assert_always(stressGoal.XZ != 0 ? (pGain.XZ != 0 || iGain.XZ != 0) : true,
                         "You need to set a gain factor in XZ in order to control stress");
    logger.assert_always(stressGoal.XX != 0 ? (pGain.XX != 0 || iGain.XX != 0) : true,
                         "You need to set a gain factor in XX in order to control stress");
    logger.assert_always(stressGoal.YY != 0 ? (pGain.YY != 0 || iGain.YY != 0) : true,
                         "You need to set a gain factor in YY in order to control stress");
    logger.assert_always(stressGoal.ZZ != 0 ? (pGain.ZZ != 0 || iGain.ZZ != 0) : true,
                         "You need to set a gain factor in ZZ in order to control stress");

    logger.assert_always(getHandler() != nullptr, "You need to set the handler of this boundary");
    const DPMBase *dpm = getHandler()->getDPMBase();

    // set gains for the strain-rate controllers
    xx.set(pGainFactor_.XX, iGainFactor_.XX);
    xy.set(pGainFactor_.XY, iGainFactor_.XY);
    yy.set(pGainFactor_.YY, iGainFactor_.YY);
    zz.set(pGainFactor_.ZZ, iGainFactor_.ZZ);
    // reset the integrated error to zero
    xx.reset();
    xy.reset();
    yy.reset();
    zz.reset();

    // if there is no shear (compression only)
    if (stressGoal_.XY == 0 && strainRate_.XY == 0)
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
    else
    {
        logger(INFO, "Shear rate is not zero, setting up Lees-Edwards in xy and a periodic boundary in z");

        PeriodicBoundary boundary;
        boundary.setHandler(getHandler());
        boundary.set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        periodicBoundaries_.push_back(boundary);

        // Lees Edwards bc in y direction & periodic boundary in x direction, simple shear boundary activated
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.setHandler(getHandler());
        double integratedShift = integratedShift_;
        leesEdwardsBoundary.set(
                [&integratedShift](Mdouble time UNUSED) { return integratedShift; },
                [](Mdouble time UNUSED) { return 0; },
                dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
        leesEdwardsBoundaries_.push_back(leesEdwardsBoundary);
    }
}

void StressStrainControlBoundary::setStrainRate(const Matrix3D& strainRate)
{
    strainRate_ = strainRate;
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
