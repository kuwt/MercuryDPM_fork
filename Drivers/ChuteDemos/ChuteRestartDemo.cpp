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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string.h>

#include "Chute.h"

class ChuteRestartDemo : public Chute
{
    void actionsOnRestart() override
    {
        setName("ChuteRestartDemo");
        setTimeMax(0.5);
        setSaveCount(100);
        fStatFile.setFileType(FileType::NO_FILE);
    }
};

/*!
 * This code tests if chute problems restart properly. As can be seen above, it is very similar to restarting any other
 * Mercury driver.
 */
int main(int argc, char* argv[])
{
    logger.assert_always(argc > 1, "Please provide a restart file to this code. Usage: ./ChuteRestartDemo -r "
            "\"ChuteDemo.restart\"");
    bool isRestart = false;
    for (unsigned int i = 0; i < argc; ++i)
    {
        if (!strcmp(argv[i], "-restart") || !strcmp(argv[i], "-r"))
        {
            isRestart = true;
        }
    }
    logger.assert_always(isRestart, "Please provide a restart file to this code. Usage: ./ChuteRestartDemo -r "
            "ChuteDemo.restart");
    ChuteRestartDemo problem;
    problem.solve(argc, argv);
}
