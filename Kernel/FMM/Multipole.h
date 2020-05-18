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
#ifndef MULTIPOLE_H_
#define MULTIPOLE_H_

#include "Math/NumericalVector.h"
#include "Math/Vector.h"
#include <complex>
#include <vector>

class Multipole
{
public:
    Multipole(int p, NumericalVector<>* squaredFactorials, Vec3D location);
    
    virtual ~Multipole();
    
    //Multipole manipulations
    virtual void computeMultipoleExpansion();
    
    NumericalVector<std::complex<Mdouble>> TranslateMultipoleExpansionTo(Vec3D location);
    
    NumericalVector<std::complex<Mdouble>> convertMultipoleToLocal(Vec3D location);
    
    //add a multipole to this existing multipole
    void addMultipoleCoefficients(NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients);
    
    NumericalVector<std::complex<Mdouble>> getExpansionCoefficients()
    {
        return multipoleExpansionCoefficients_;
    }
    
    void setExpansionCoefficients(NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients)
    {
        multipoleExpansionCoefficients_ = multipoleExpansionCoefficients;
    }
    
    NumericalVector<>* getSquaredFactorials()
    {
        return squaredFactorials_;
    }
    
    int getP()
    {
        return p_;
    }

protected:
    int p_; // order of truncation
    NumericalVector<>* squaredFactorials_; // required for multipole manipulations
    Vec3D location_; // location of the multipole
    NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients_; // coefficients of the multipole
};

#endif /* MULTIPOLE_H_ */
