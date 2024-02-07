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

#include "DisplacementField.h"
#include "Particles/BaseParticle.h"
#include <cmath>

namespace CGFields
{
    DisplacementField::DisplacementField()
    {
        setZero();
        logger(DEBUG, "DisplacementField::DisplacementField() finished");
    }
    
    void DisplacementField::writeNames(std::ostream& os, unsigned countVariables)
    {
        os << countVariables + 1 << "-" << countVariables + 3 << ":displacementMomentum " //orientation
           << countVariables + 4 << "-" << countVariables + 9 << ":displacementMomentumFlux ";
    }
    
    void DisplacementField::write(std::ostream& os) const
    {
        os << displacementMomentum_ << " "
           << displacementMomentumFlux_ << " ";
    }
    
    void DisplacementField::output(std::ostream& os) const
    {
        os << "displacementMomentum " << displacementMomentum_ << " ";
        os << "displacementMomentumFlux " << displacementMomentumFlux_ << " ";
    }
    
    void DisplacementField::setZero()
    {
        displacementMomentum_.setZero();
        displacementMomentumFlux_.setZero();
    }
    
    DisplacementField DisplacementField::getSquared() const
    {
        DisplacementField DisplacementField;

        DisplacementField.displacementMomentum_ = Vec3D::square(displacementMomentum_);
        DisplacementField.displacementMomentumFlux_ = MatrixSymmetric3D::square(displacementMomentumFlux_);
        return DisplacementField;
    }
    
    DisplacementField& DisplacementField::operator=(const DisplacementField& P)= default;
    
    DisplacementField& DisplacementField::operator+=(const DisplacementField& P)
    {
        displacementMomentum_ += P.displacementMomentum_;
        displacementMomentumFlux_ += P.displacementMomentumFlux_;
        return *this;
    }
    
    DisplacementField& DisplacementField::operator-=(const DisplacementField& P)
    {
        displacementMomentum_ -= P.displacementMomentum_;
        displacementMomentumFlux_ -= P.displacementMomentumFlux_;
        return *this;
    }
    
    DisplacementField& DisplacementField::operator/=(Mdouble a)
    {
        displacementMomentum_ /= a;
        displacementMomentumFlux_ /= a;
        return *this;
    }
    
    DisplacementField DisplacementField::operator*(Mdouble a) const
    {
        DisplacementField p;
        p.displacementMomentum_ = displacementMomentum_ * a;
        p.displacementMomentumFlux_ = displacementMomentumFlux_ * a;
        return p;
    }
    
    void DisplacementField::addParticleStatistics(Mdouble phi, const DisplacementField& currentInteraction)
    {
        displacementMomentum_ += currentInteraction.getDisplacementMomentum() * phi;
        displacementMomentumFlux_ += currentInteraction.getDisplacementMomentumFlux() * phi;
    }
    
    void DisplacementField::setFields(const BaseParticle& p)
    {
        CGHandler* cgHandler = getCG()->getHandler();
        DPMBase* dpmBase = getCG()->getHandler()->getDPMBase();  
        
        Vec3D displacement = p.getDisplacement2(dpmBase->getXMin(), dpmBase->getXMax(), 
                                                dpmBase->getYMin(), dpmBase->getYMax(), 
                                                dpmBase->getZMin(), dpmBase->getZMax(),
                                                dpmBase->getTime() - cgHandler->getPreviousEvaluationTime());
        
        displacementMomentum_ = displacement * p.getMass();
        displacementMomentumFlux_ = MatrixSymmetric3D::selfDyadic(displacement) * p.getMass();
        
    }
    
    
    void DisplacementField::setCylindricalFields(const BaseParticle& p)
    {
        setFields(p);
    }
    
}
