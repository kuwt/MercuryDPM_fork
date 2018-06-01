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
 * \details constructor
 */
StressStrainControlBoundary::StressStrainControlBoundary()
        : BaseBoundary()
{
    stressGoal_.setZero();
    strainRate_.setZero();
    gainFactor_.setZero();
    isStrainRateControlled_ = true;
    stressTotal.setZero();
    stressKinetic.setZero();
    stressStatic.setZero();
    dstrainRate.setZero();
    lengthBox.setZero();
    centerBox.setZero();
    relativeToCenter.setZero();
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
    for (const LeesEdwardsBoundary& b : leesEdwardsBoundary_)
    {
        os << b << "\n";
    }
    for (const PeriodicBoundary& b : periodicBoundary_)
    {
        os << b << "\n";
    }
    BaseBoundary::write(os);
    
    
    os << "\n";
    os << " stressGoal " << stressGoal_;
    os << " strainRate " << strainRate_;
    os << " gainFactor " << gainFactor_;
    os << " isStrainRateControlled " << isStrainRateControlled_;
    os << " integratedShift " << integratedShift_;
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in] is        the istream
 */
void StressStrainControlBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> stressGoal_;
    is >> dummy >> strainRate_;
    is >> dummy >> gainFactor_;
    is >> dummy >> isStrainRateControlled_;
    is >> dummy >> integratedShift_;
    set(stressGoal_, strainRate_, gainFactor_, isStrainRateControlled_);
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
 */
void StressStrainControlBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    // Call Boundaries
    for (PeriodicBoundary& b : periodicBoundary_)
    {
        b.checkBoundaryAfterParticlesMove(pH);
    }
    for (LeesEdwardsBoundary& b : leesEdwardsBoundary_)
    {
        b.checkBoundaryAfterParticlesMove(pH);
    }
    //Real Stress and StrainRate Control every time step
    logger.assert_always(getHandler() != nullptr,
                         "you need to set the handler of this boundary before the parameters can be set");
    DPMBase* dpm = getHandler()->getDPMBase();
    
    
    lengthBox.X = (dpm->getXMax() - dpm->getXMin());
    lengthBox.Y = (dpm->getYMax() - dpm->getYMin());
    lengthBox.Z = (dpm->getZMax() - dpm->getZMin());
    centerBox.X = (dpm->getXMax() + dpm->getXMin()) / 2.0;
    centerBox.Y = (dpm->getYMax() + dpm->getYMin()) / 2.0;
    centerBox.Z = (dpm->getZMax() + dpm->getZMin()) / 2.0;// Box length Lx, Lyand Lz, and center point C of the box
    //static double integratedShift = 0.0; //!Used for updating the shift of Lees-Edwards boundary
    
    
    if (stressGoal_.XX == 0 && stressGoal_.YY == 0 && stressGoal_.ZZ == 0 &&
        stressGoal_.XY == 0) //this activate only the strainrate control in x-y-z
    {
        //  Change the system size according to next time step
        dpm->setXMax(dpm->getXMax() + 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
        dpm->setXMin(dpm->getXMin() - 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
        dpm->setYMax(dpm->getYMax() + 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
        dpm->setYMin(dpm->getYMin() - 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
        dpm->setZMax(dpm->getZMax() + 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
        dpm->setZMin(dpm->getZMin() - 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
        
        // Box length Lx, Lyand Lz for next time step
        lengthBox.X = (dpm->getXMax() - dpm->getXMin());
        lengthBox.Y = (dpm->getYMax() - dpm->getYMin());
        lengthBox.Z = (dpm->getZMax() - dpm->getZMin());
        
        //  Move the boundaries to next time step
        if (strainRate_.XY != 0 || stressGoal_.XY != 0)
        {
            //Determine the current shear velocity and shift of the Lees-Edwards boundary
            double velocity_xy = strainRate_.XY * lengthBox.Y;
            integratedShift_ += velocity_xy * dpm->getTimeStep();
            
            // Determine how to move boundaries based on strainrate control or boundary control
            if (isStrainRateControlled_)
            {
                //  Move the Lees-Edwards boundary in z direction to next time step
                leesEdwardsBoundary_[0].set(
                        [this](double time)
                        { return integratedShift_; },
                        [](double time UNUSED)
                        { return 0; },
                        dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
            }
            else
            {
                //  Update the velocity of Lees-Edwards boundary in x-y direction to next time step
                leesEdwardsBoundary_[0].set(
                        [this](double time)
                        { return integratedShift_; },
                        [velocity_xy](double time UNUSED)
                        { return velocity_xy; },
                        dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
            }
            std::cout << "leesEdwardsShift" << leesEdwardsBoundary_[0].getCurrentShift();
            //  Move the boundary in z direction to next time step
            periodicBoundary_[0].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
            
        }
        else
        {
            //  Move the boundary in x direction to next time step
            periodicBoundary_[0].set(Vec3D(1.0, 0.0, 0.0), dpm->getXMin(), dpm->getXMax());
            //  Move the boundary in y direction to next time step
            periodicBoundary_[1].set(Vec3D(0.0, 1.0, 0.0), dpm->getYMin(), dpm->getYMax());
            //  Move the boundary in z direction to next time step
            periodicBoundary_[2].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        }
        
        //  Give the strain-rate for all particles and move them to next timestep before integration
        if (isStrainRateControlled_)
        {
            for (auto& p : pH)
            {
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
    else
    {
        dstrainRate.setZero();
        
        
        //calculate stress for kinetic part
        stressKinetic = dpm->getKineticStress();
        
        //calculate the static stress tensor
        stressStatic = dpm->getStaticStress();
        
        // calculate the stress total and average over the volume
        stressTotal = dpm->getTotalStress();
        
        // amount by which the strainrate tensor has to be changed
        if (stressGoal_.XX != 0)
        {
            dstrainRate.XX = dstrainRate.XX + gainFactor_.XX * dpm->getTimeStep() * (stressTotal.XX - stressGoal_.XX);
            std::cout << "StressXX = " << stressTotal.XX << std::endl;
        }
        if (stressGoal_.YY != 0)
        {
            dstrainRate.YY = dstrainRate.YY + gainFactor_.YY * dpm->getTimeStep() * (stressTotal.YY - stressGoal_.YY);
            std::cout << "StressYY = " << stressTotal.YY << std::endl;
        }
        if (stressGoal_.ZZ != 0)
        {
            dstrainRate.ZZ = dstrainRate.ZZ + gainFactor_.ZZ * dpm->getTimeStep() * (stressTotal.ZZ - stressGoal_.ZZ);
            std::cout << "StressZZ = " << stressTotal.ZZ << std::endl;
        }
        if (stressGoal_.XY != 0)
        {
            dstrainRate.XY = dstrainRate.XY + gainFactor_.XY * dpm->getTimeStep() * (stressTotal.XY - stressGoal_.XY);
            std::cout << "dstrainRate.XY = " << dstrainRate.XY << std::endl;
            std::cout << "StressXY = " << stressTotal.XY << std::endl;
            std::cout << "StressYX = " << stressTotal.YX << std::endl;
        }
        
        
        //  Update the strainrate tensor
        strainRate_ = strainRate_ + dstrainRate;
        
        //  Change the system size according to next time step
        dpm->setXMax(dpm->getXMax() + 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
        dpm->setXMin(dpm->getXMin() - 0.5 * lengthBox.X * strainRate_.XX * dpm->getTimeStep());
        dpm->setYMax(dpm->getYMax() + 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
        dpm->setYMin(dpm->getYMin() - 0.5 * lengthBox.Y * strainRate_.YY * dpm->getTimeStep());
        dpm->setZMax(dpm->getZMax() + 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
        dpm->setZMin(dpm->getZMin() - 0.5 * lengthBox.Z * strainRate_.ZZ * dpm->getTimeStep());
        
        // Box length Lx, Lyand Lz for next time step
        lengthBox.X = (dpm->getXMax() - dpm->getXMin());
        lengthBox.Y = (dpm->getYMax() - dpm->getYMin());
        lengthBox.Z = (dpm->getZMax() - dpm->getZMin());
        
        //  Move the boundaries to next time step
        if (strainRate_.XY != 0 || stressGoal_.XY != 0)
        {
            //Determine the current shear velocity and shift of the Lees-Edwards boundary
            double velocity_xy = strainRate_.XY * lengthBox.Y;
            integratedShift_ += velocity_xy * dpm->getTimeStep();
            
            // Determine how to move boundaries based on strainrate control or boundary control
            if (isStrainRateControlled_)
            {
                //  Move the Lees-Edwards boundary in z direction to next time step
                leesEdwardsBoundary_[0].set(
                        [this](double time)
                        { return integratedShift_; },
                        [](double time UNUSED)
                        { return 0; },
                        dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
            }
            else
            {
                //  Update the velocity of Lees-Edwards boundary in x-y direction to next time step
                leesEdwardsBoundary_[0].set(
                        [this](double time)
                        { return integratedShift_; },
                        [velocity_xy](double time UNUSED)
                        { return velocity_xy; },
                        dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
            }
            //  Move the boundary in z direction to next time step
            periodicBoundary_[0].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        }
        else
        {
            //  Move the boundary in x direction to next time step
            periodicBoundary_[0].set(Vec3D(1.0, 0.0, 0.0), dpm->getXMin(), dpm->getXMax());
            //  Move the boundary in y direction to next time step
            periodicBoundary_[1].set(Vec3D(0.0, 1.0, 0.0), dpm->getYMin(), dpm->getYMax());
            //  Move the boundary in z direction to next time step
            periodicBoundary_[2].set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        }
        //  Give the strain-rate for all particles and move them to next timestep before integration
        if (isStrainRateControlled_)
        {
            for (auto& p : pH)
            {
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
}


void
StressStrainControlBoundary::set(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
                                 bool isStrainRateControlled)
{
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
        periodicBoundary_.push_back(boundary);
        boundary.set(Vec3D(0.0, 1.0, 0.0), dpm->getYMin(), dpm->getYMax());
        periodicBoundary_.push_back(boundary);
        boundary.set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        periodicBoundary_.push_back(boundary);
    }
    else if (stressGoal.XY != 0 || strainRate.XY != 0)
    {
        logger(INFO, "Shear rate is not zero, setting up Lees-Edwards in xy and a periodic boundary in z");
        
        PeriodicBoundary boundary;
        boundary.setHandler(getHandler());
        boundary.set(Vec3D(0.0, 0.0, 1.0), dpm->getZMin(), dpm->getZMax());
        periodicBoundary_.push_back(boundary);
        // Lees Edwards bc in y direction & periodic boundary in x direction, simple shear boundary activated
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.setHandler(getHandler());
        leesEdwardsBoundary.set(
                [this](double time UNUSED)
                { return integratedShift_; },
                [](double time UNUSED)
                { return 0; },
                dpm->getXMin(), dpm->getXMax(), dpm->getYMin(), dpm->getYMax());
        leesEdwardsBoundary_.push_back(leesEdwardsBoundary);
    }
    else
    {
        logger(ERROR, "This should not happen; please check implementation of StressStrainControlBoundary::set");
    }
    //std::cout<< "leesEdwardsShift" << leesEdwardsBoundary_[0].getCurrentShift();
    
}

void StressStrainControlBoundary::createPeriodicParticles(ParticleHandler& pH)
{
    // Call Boundaries
    for (PeriodicBoundary& b : periodicBoundary_)
    {
        b.createPeriodicParticles(pH);
    }
    for (LeesEdwardsBoundary& b : leesEdwardsBoundary_)
    {
        b.createPeriodicParticles(pH);
    }
}
















