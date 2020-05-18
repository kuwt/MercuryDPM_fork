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

#include <Math/Vector.h>
#include <Math/Helpers.h>
#include "Logger.h"

int main()
{
    logger(INFO,"Checks several Vec3D unit operations");
    Vec3D a = {1,2,3};
    Vec3D b = {4,5,6};
    Vec3D c;
    
    c = b - a; //same as c = b.operator-(a);
    helpers::check(c,{3,3,3},0,"Vector subtraction");
    
    c = - a; //same as c = operator-(a);
    helpers::check(c,{-1,-2,-3},0,"Negation of a vector");
    
    c -= a; //same as c.operator-=(a);
    helpers::check(c,{-2,-4,-6},0,"Vector self-subtraction");
    
    c = b + a; //same as c = b.operator-(a);
    helpers::check(c,{5,7,9},0,"Vector addition");
    
    c = a; //same as c = operator+(a); //missing operator
    //helpers::check(c,{1,2,3},0,"Plus-operation of a vector");
    
    c += a; //same as c.operator-=(a);
    helpers::check(c,{2,4,6},0,"Vector self-addition");
    
    helpers::check(Vec3D::dot(a,b),32,0,"Dot product");
    
    helpers::check(Vec3D::cross(a,b),{-3,6,-3},0,"Cross product");
    
    return 0;
}



