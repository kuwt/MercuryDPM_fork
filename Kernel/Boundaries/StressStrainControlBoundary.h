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
 * \brief 
 * \details Inherits from BaseObject
 */
class StressStrainControlBoundary : public BaseBoundary
{
public:
    /*!
     * \brief default constructor.
     */
    StressStrainControlBoundary();
    
    /*!
     * \brief copy constructor
     */
    StressStrainControlBoundary(const StressStrainControlBoundary& b) = default;
    
    /*!
     * \brief destructor
     */
    virtual ~StressStrainControlBoundary() = default;
    
    /*!
     * \brief Reads the object's id_ from given istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Adds object's id_ to given ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Sets the name of the boundary
     */
    std::string getName() const;
    
    /*!
     * \brief Used to create a copy of the object
     */
    virtual StressStrainControlBoundary* copy() const;
    
    /*!
     * \brief Virtual function that does things to particles, each timestep after particles have moved.
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    void set(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
             bool isStrainRateControlled);
    
    void createPeriodicParticles(ParticleHandler& pH) override;


private:
    
    /**
     * \brief Stores the stress value the boundary should attain.
     * \details Unused if values are set to nan. Default value is ...
     */
    Matrix3D stressGoal_, strainRate_, gainFactor_;
    
    Matrix3D stressTotal;  //!stress components calculation variables
    Matrix3D stressKinetic;  //!stress components calculation variables
    Matrix3D stressStatic;  //!stress components calculation variables
    Matrix3D dstrainRate; //!Stress control strainrate tensor change per timestep
    Vec3D lengthBox; //!Box length in x-y-z
    Vec3D centerBox; //!center position of the box
    Vec3D relativeToCenter; //!particle position relative to the center of box
    
    std::vector<LeesEdwardsBoundary> leesEdwardsBoundary_;
    
    std::vector<PeriodicBoundary> periodicBoundary_;
    
    bool isStrainRateControlled_;
    
    double integratedShift_;
};

#endif
