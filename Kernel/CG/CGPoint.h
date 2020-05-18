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
#ifndef CGPOINT_H
#define    CGPOINT_H

#include "CG/Fields/StandardFields.h"
#include "CG/Fields/OrientationField.h"

/*!
 * \brief Combines the position of the CGPoint (e.g. x, y, z), the parameters 
 * of the coarse-graining function (e.g. width and cutoff) and the fields to be 
 * evaluated (e.g., density, momentum, stress).
 * \details The class is combines the properties of a StandardFields and a CGFunctions
 * class.
 * It  contains two functions that depend on both the fields and the cg function,
 * evaluateParticle and evaluateContact.
 * 
 * See StandardFields and CGFunctions for more details.
 */
template<class Coordinates, class Fields=CGFields::StandardFields>
class CGPoint : public Fields
{
public:
    
    typedef Coordinates CoordinatesType;
    
    CGPoint() = default;
    
    CGPoint(const CGPoint& orig) = default;
    
    virtual ~CGPoint() = default;
    
    /*!
     * \brief Combines the write functions of the two base classes Coordinates and Fields.
     */
    void write(std::ostream& os) const
    {
        coordinates.write(os);
        Fields::write(os);
    }

public:
    
    Coordinates coordinates;
    
};

#endif

