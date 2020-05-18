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

#ifndef CUBEDELETIONBOUNDARY_H
#define CUBEDELETIONBOUNDARY_H

#include "DeletionBoundary.h"
#include "BaseBoundary.h"
#include "Math/Vector.h"

class ParticleHandler;

class BaseParticle;

class CubeDeletionBoundary : public DeletionBoundary
{
public:
    /*!
     * \brief default constructor
     */
    CubeDeletionBoundary();
    
    /*!
     * \brief destructor
     */
    ~CubeDeletionBoundary() override;
    
    /*!
     * \brief Copy method; creates copy on the heap and returns a pointer to it.
     */
    CubeDeletionBoundary* copy() const override;
    
    /*!
     * \brief Sets boundary position based on two opposite corners.
     */
    void set(Vec3D posMin, Vec3D posMax);
    
    /*!
     * \brief Returns a negative value if and only if the particle is inside the boundary
     */
    Mdouble getDistance(const Vec3D& position) const override;
    
    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

private:
    
    /*!
     * \brief Minimal and maximal positions defining the boundary's boundaries.
     */
    Vec3D posMin_, posMax_;
};

#endif
