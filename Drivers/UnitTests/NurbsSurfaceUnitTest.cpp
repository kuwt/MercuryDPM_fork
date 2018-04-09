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

#include "Nurbs/NurbsSurface.h"
#include "Math/Helpers.h"
#include <sstream>

int main()
{
    //define quarter circle
    std::vector<double> knotsU = {0,0,0,1,1,1};
    std::vector<double> knotsV = {0,0,1,1};
    std::vector<std::vector<Vec3D>> controlPoints = {{{1,0,0},{1,0,2}},{{1,1,0},{1,1,2}},{{0,1,0},{0,1,2}}} ;
    std::vector<std::vector<Mdouble>> weights = {{1,1},{1,1},{2,2}};
    NurbsSurface nurbs(knotsU,knotsV,controlPoints,weights);
    Vec3D val = nurbs.evaluate(0.5,0.5);
    helpers::check(val,{.6,.8,1},std::numeric_limits<double>::epsilon(),"Nurbs evaluation failed");
    std::cout << val << std::endl;

//    unsigned n=100;
//    std::stringstream ss;
//    for (unsigned i=0; i<=n; ++i) {
//        for (unsigned j = 0; j <= n; ++j) {
//            ss << nurbs.evaluate((double)i/n, (double)j/n) << '\n';
//        }
//    }
//    helpers::writeToFile("NurbsSurfaceUnitTest.txt",ss.str());
    //std::cout << ss.str();
    return 0;
}
