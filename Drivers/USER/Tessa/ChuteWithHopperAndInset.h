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

#ifndef CHUTEWITHHOPPERANDINSET_H
#define CHUTEWITHHOPPERANDINSET_H
#include "ChuteWithHopper.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

class ChuteWithHopperAndInset : public ChuteWithHopper
{
    
///protected variables
protected:
    
/// private variables
private:
    /// The height of the inset
    Mdouble insetHeight;
    ///The width of the inset
    Mdouble insetWidth;
    /// The angle of the inset (input in degrees, usage in radians)
    Mdouble insetAngle;

/// public variables
public:
    
    ///The default constructor  
    ChuteWithHopperAndInset();
    
    ///The actually constructor
    void constructor();
    ///todo check wether the inset and hopper are colliding 
    ///todo check wether the resulting opening between the hopper and 
    ///the inset is 
    
    /// set function for insetHeight, insetWidth, insetAngle
    void set_Inset(double height, double width, double angle);
    
    /// get function for insetHeight, insetWidth, insetAngle
    double get_InsetHeight();
    double get_InsetWidth();
    double get_InsetAngle();
    
    virtual void setupInitialConditions();
    
    void add_Inset();
    
    virtual void write(std::ostream& os, bool print_all);
    /// Give a warning somewhere when the hopper is not raised (variable: ChuteWithHopper::lift) enough and the finite wall (inset) and hopper are colliding

    LinearViscoelasticFrictionSpecies* species;
};
#endif
