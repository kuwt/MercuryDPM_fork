//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include "OrientationField.h"
#include "Particles/BaseParticle.h"
#include <cmath>

namespace CGFields
{
    OrientationField::OrientationField()
    {
        setZero();
        logger(DEBUG, "OrientationField::OrientationField() finished");
    }
    
    void OrientationField::writeNames(std::ostream& os, unsigned countVariables)
    {
        os << countVariables + 1 << ":Orientation "; //orientation
    }
    
    void OrientationField::write(std::ostream& os) const
    {
        os << orientation_;
    }
    
    void OrientationField::output(std::ostream& os) const
    {
        os << "Orientation " << orientation_;
    }
    
    void OrientationField::setZero()
    {
        orientation_.setZero();
    }
    
    OrientationField OrientationField::getSquared() const
    {
        OrientationField orientationField;
        orientationField.orientation_ = MatrixSymmetric3D::square(orientation_);
        return orientationField;
    }
    
    OrientationField& OrientationField::operator=(const OrientationField& P)= default;
    
    OrientationField& OrientationField::operator+=(const OrientationField& P)
    {
        orientation_ += P.orientation_;
        return *this;
    }
    
    OrientationField& OrientationField::operator-=(const OrientationField& P)
    {
        orientation_ -= P.orientation_;
        return *this;
    }
    
    OrientationField& OrientationField::operator/=(Mdouble a)
    {
        orientation_ /= a;
        return *this;
    }
    
    OrientationField OrientationField::operator*(Mdouble a) const
    {
        OrientationField p;
        p.orientation_ = orientation_ * a;
        return p;
    }
    
    void OrientationField::addParticleStatistics(Mdouble phi, const OrientationField& currentInteraction)
    {
        orientation_ += currentInteraction.getOrientation() * phi;
    }
    
    void OrientationField::setFields(const BaseParticle& p)
    {
        Vec3D orientation = p.getOrientation().getAxis();
        orientation_.XX = orientation.X * orientation.X;
        orientation_.XY = orientation.X * orientation.Y;
        orientation_.YY = orientation.Y * orientation.Y;
        orientation_.XZ = orientation.X * orientation.Z;
        orientation_.YZ = orientation.Y * orientation.Z;
        orientation_.ZZ = orientation.Z * orientation.Z;
        logger(DEBUG, "orientation: %", orientation_);
    }
    
    
    void OrientationField::setCylindricalFields(const BaseParticle& p)
    {
        setFields(p);
    }
    
}
