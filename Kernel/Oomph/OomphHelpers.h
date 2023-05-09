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

#ifndef OOMPHHELPERS_H
#define OOMPHHELPERS_H
#include "Vector.h"

namespace oomph {
    /*!
     * \details Adds all elements of the vector to an output stream.
     * NB: this is a global function and a friend of the Vec3D class!
     * \param[in] os    the output stream,
     * \param[in] a     The vector of interest
     * \return          the output stream with vector elements added
     */
    std::ostream &operator<<(std::ostream &os, const Vector<double> &a) {
        if (not a.empty()) {
            os << a[0];
        }
        for (int i = 1; i < a.size(); i++) {
            os << ' ' << a[i];
        }
        return os;
    }
}

namespace convertVecFuncs
{
    Vec3D convertToVec3D(oomph::Vector<double>& Vec)
    {
        Vec3D returnVec;
        returnVec.X = Vec[0];
        returnVec.Y = Vec[1];
        returnVec.Z = Vec[2];
        return returnVec;
    }

    oomph::Vector<double> convertToOomphVec(const Vec3D& Vec)
    {
        oomph::Vector<double> returnVec(3,0.0);
        returnVec[0] = Vec.X;
        returnVec[1] = Vec.Y;
        returnVec[2] = Vec.Z;
        return returnVec;
    }
}

#endif //OOMPHHELPERS_H
