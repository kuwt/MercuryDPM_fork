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

#include <Species/LinearViscoelasticSpecies.h>
#include "DPMBase.h"

/**
 * Test whether savecount can be correctly reset.
 */
class SaveCountUnitTest : public DPMBase
{
public:
    void setupInitialConditions() override
    {
        setTimeStep(1);
        setTimeMax(40);
        setSaveCount(10);
        setDomain({0,0,0},{1,1,1});
        speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

    void actionsBeforeTimeStep() override
    {
        if (getNumberOfTimeSteps()==18) {
            logger(INFO,"changing savecount to 5 at time %",getTime());
            setSaveCount(5);
        } else if (getNumberOfTimeSteps()==33) {
            logger(INFO,"changing savecount to 2 at time %",getTime());
            setSaveCount(2);
        }
    }

    SaveCountUnitTest() : DPMBase()
    {
        setName("SaveCountUnitTest");
        solve();

        //check unit test
        logger.assert_always(dataFile.getCounter()==9,
                             "Number of time steps written is %, should be 9", dataFile.getCounter());
        logger.assert_always(dataFile.getLastSavedTimeStep(),
                             "Last time steps written is %, should be 40", dataFile.getLastSavedTimeStep());
        //if nothing fails, write success message
        logger(INFO,"Unit test % finished successfully",getName());
    }
};

int main(int argc, char* argv[])
{
    SaveCountUnitTest problem;
    return 0;
}
