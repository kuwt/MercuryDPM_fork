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

#ifndef STLTRIANGLE_H
#define STLTRIANGLE_H

#include "Math/Vector.h"

class STLTriangle
{
public:
    STLTriangle()
    {};
    
    STLTriangle(const Vec3D newNormal, const Vec3D newVertex1, const Vec3D newVertex2, const Vec3D newVertex3)
    {
        normal = newNormal;
        vertex1 = newVertex1;
        vertex2 = newVertex2;
        vertex3 = newVertex3;
    }
    
    bool isEqualTo(const STLTriangle answer, double toll)
    {
        if (!normal.isEqualTo(answer.normal, toll)) return false;
        if (!vertex1.isEqualTo(answer.vertex1, toll)) return false;
        if (!vertex2.isEqualTo(answer.vertex2, toll)) return false;
        if (!vertex3.isEqualTo(answer.vertex3, toll)) return false;
        return true;
    }
    
    Vec3D normal;
    Vec3D vertex1;
    Vec3D vertex2;
    Vec3D vertex3;
    
};

#endif //STLTRIANGLE_H
