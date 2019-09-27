//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef StressStrainControlBoundary_H
#define StressStrainControlBoundary_H

#include "BaseBoundary.h"
#include "Math/ExtendedMath.h"
#include "Side.h"
#include "LeesEdwardsBoundary.h"
#include "PeriodicBoundary.h"

class PeriodicBoundaryHandler;

/*!
 * \class StressStrainControlBoundary
 * \brief Perodic boundary that can be strain/stress controlled and achieve different
 * deformation modes.
 * \details Inherits from BaseObject, this boundary take matricies inputs from user
 * and determine which boundaries should be combined together and therefore achieve
 * different deformation mode: e.g. triaxial, bi-aixial, uni-axial and simple shear.
 */
class StressStrainControlBoundary : public BaseBoundary
{
public:
    /// \brief default constructor
    StressStrainControlBoundary();
    
    ///\brief copy constructor
    StressStrainControlBoundary(const StressStrainControlBoundary& b) = default;
    
    /// \brief destructor
    virtual ~StressStrainControlBoundary() = default;
    
    /// \brief Reads the object's id_ from given istream
    void read(std::istream& is) override;
    
    /// \brief Adds object's id_ to given ostream
    void write(std::ostream& os) const override;
    
    /// \brief Sets the name of the boundary
    std::string getName() const override;
    
    /// \brief Used to create a copy of the object
    StressStrainControlBoundary* copy() const override;
    
    /// \brief Virtual function that does things to particles, each timestep after particles have moved.
    void checkBoundaryAfterParticlesMove(ParticleHandler& particleHandler) override;
    
    /// \brief Sets the inputs from user to initialize the boundaries
    void set(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
             bool isStrainRateControlled);

    //void actionsBeforeTimeLoop() override;

    /// \brief Create the periodic particles after read in from a restart file to attain right information
    void createPeriodicParticles(ParticleHandler& particleHandler) override;
    
    //helper functions
    
    /// \brief Call the boundary and update them based on the new domain size after the stress-control movement
    void checkPeriodicLeesEdwardsBoundariesAfterParticlesMove(ParticleHandler& particleHandler);
    
    /// \brief Determine the length in x,y,z and center location of domain
    void determineLengthAndCentre();
    
    /// \brief Activate the strainrate control based on user's boolean input
    void activateStrainRateControl(const ParticleHandler& particleHandler);
    
    /// \brief compute the change of strainrate tensor based on the stress difference and then update the tensor
    void computeStrainRate();
    
    /// \brief Determines stress-controlled shear Lees-Edwards boundary in x-y direction and normal periodic in z direction
    void determineStressControlledShearBoundaries();

    /// \brief update the domain to new sizes
    void updateDomainSize();

    /// \brief accesses the strain rate tensor
    Matrix3D getStrainRate() const {return strainRate_;}

    /// \brief accesses the stressGoal_
    Matrix3D getStressGoal() const {return stressGoal_;}

    /// \brief accesses the gainFactor
    Matrix3D getGainFactor() const {return gainFactor_;}

private:

    //Set by the user

    /*!
     * \brief Stores the stress value the boundary should attain.
     * \details Unused if the stressGoal values are set to zero.
     */
    Matrix3D stressGoal_, strainRate_, gainFactor_;
    /// The boolean input, true means switch on the strain rate control
    bool isStrainRateControlled_;

    //Set each time step in checkBoundariesAfterParticlesMove

    /// Box length in x-y-z
    Vec3D lengthBox_;
    /// Center position of the domain
    Vec3D centerBox_;
    /// Particle position relative to the center of domain
    Vec3D relativeToCenter_;

    // Set in the constructor, incremented in checkBoundariesAfterParticlesMove

    /// Shift integrated for all the time when using Lees-Edwards Boundary
    Mdouble integratedShift_;

    //Defined in the set function

    /// Store boundaries into a vector for the pushback
    /// Note, there is always either no LeesEdwardsBoundary and 3 PeriodicBoundary,
    /// or 1 LeesEdwardsBoundary and 1 PeriodicBoundary
    std::vector<LeesEdwardsBoundary> leesEdwardsBoundaries_;
    std::vector<PeriodicBoundary> periodicBoundaries_;
};

#endif
