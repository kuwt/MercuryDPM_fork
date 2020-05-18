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

#include "Math/Helpers.h"
#include "Math/ExtendedMath.h"
#include <sstream>

int main()
{
    std::stringstream ss;
    for (Mdouble x=-6; x<6; x+=0.01)
    {
        Mdouble e = mathsFunc::exp(x);
        Mdouble l = mathsFunc::log(e);
        Mdouble c = mathsFunc::cos(x);
        Mdouble s = mathsFunc::sin(x);
        ss << x << " "
           << e << " "
           << l << " "
           << c << " "
           << s << "\n";
    }
    std::cout << "Execute 'gnuplot MathCheck.gnu' to view output" << std::endl;
    helpers::writeToFile("MathCheck.data",ss.str());
    helpers::writeToFile("MathCheck.gnu",
                         "set xlabel 'x'\n"
                         "p 'MathCheck.data' u 1:(log($2)) t 'log(Exp(x))',"
                         "  'MathCheck.data' u 1:3 t 'Log(Exp(x))', x,"
                         "  'MathCheck.data' u 1:5 t 'Sin(x)', sin(x),"
                         "  'MathCheck.data' u 1:4 t 'Cos(x)', cos(x)x"
                         "\n");
    //WARNING: sin and cos are really bad outside (-pi,pi)
    Mdouble d = 1e-15;
    std::cout << "Accuracy of sin: " << mathsFunc::sin(constants::pi-d)-mathsFunc::sin(constants::pi+d) << std::endl;
    std::cout << "Accuracy of cos: " << mathsFunc::cos(constants::pi)+1 << std::endl;
}
