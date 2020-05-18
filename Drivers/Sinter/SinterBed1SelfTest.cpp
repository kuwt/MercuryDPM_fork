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

#include "Mercury3D.h"
#include <Species/SinterFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <iomanip>
#include <cassert>

using helpers::to_string;
using helpers::writeToFile;
using helpers::readFromFile;

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class Sintering : public Mercury3D
{
public:

    //set default values
    Sintering()
    {
        setName("SinterBed0");
        readRestartFile();
        setRestarted(false);
        setName("SinterBed1");
        species = dynamic_cast<SinterFrictionSpecies *>(speciesHandler.getObject(0));
        assert(species);
    }

    /** creates custom console output */
    void printTime() const override
    {
        static int counter = 0;
        if (++counter == 1)
            writeEneHeader(std::cout);
        writeEneTimeStep(std::cout);
    }

    /** creates custom ene header */
    void writeEneHeader(std::ostream &os) const override
    {
        os << "relTime\trelMeanOverlap\n";
    }

    Mdouble getMeanPlasticOverlap() const
    {
        Mdouble sum = 0;
        for (BaseInteraction* const p : interactionHandler)
        {
            auto q = dynamic_cast<SinterInteraction*>(p);
            logger.assert_always(q!= nullptr,"no SinterInteraction type");
            sum += q->getPlasticOverlap();
        }
        return sum / interactionHandler.getNumberOfObjects();
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanDiameter = 2.0*particleHandler.getMeanRadius();
        //logger(INFO,"d % r %",meanPlasticOverlap,meanRadius);
        return sqrt(meanOverlap/meanDiameter);
    }

    /** creates custom ene output */
    void writeEneTimeStep(std::ostream &os) const override
    {
        os << std::setw(12) << getTime()
           << "\t" << std::setw(12) << getMeanRelativeContactRadius()
           << std::endl;
    }

//    void actionsAfterTimeStep() override
//    {
//        if (getTime()>0.1*getTimeMax())
//        {
//            setSaveCount(200);
//            restartFile.setSaveCount(getTimeMax()/getTimeStep()/5.0+1.0);
//        }
//    }

    SinterFrictionSpecies *species; //pointer to the particle properties
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Sintering s;
    s.setSaveCount(100);
    s.fStatFile.setSaveCount(300);

    Mdouble sinterRate = readFromFile("in","sinterRate",0.25*0.25*mathsFunc::square(0.3));
    s.species->setSinterRate(sinterRate);

    Mdouble sinterTime = readFromFile("in","sinterTime",24);
    s.setTimeMax(sinterTime);
    s.restartFile.setSaveCount(s.getTimeMax()/s.getTimeStep()/5.0+1.0);
    s.fStatFile.setFileType(FileType::ONE_FILE);
    s.restartFile.setFileType(FileType::MULTIPLE_FILES); // so we can restart from different times

    writeToFile("SinterBed1.gnu",
                          "set xlabel 't'\n"
                          "set ylabel 'x/d'\n"
                          "p '"+s.getName()+".ene' u 1:2 w lp, 0.075*x**0.5\n"
    );
    writeToFile("SinterBed1b.gnu",
                          "set xlabel 't'\n"
                          "set ylabel 'x/d'\n"
                          "p '"+s.getName()+".fstat' u 1:(sqrt($7/"+to_string(2.0*s.particleHandler.getMeanRadius())+")), 0.3*x**0.5\n"
    );

    s.species->setSinterType(SINTERTYPE::CONSTANT_RATE);
    s.solve();

    std::cout << "Run with ./SinterBed1 >out; the execute 'gnuplot SinterBed1.gnu' to view output" << std::endl;

    //if I add friction (10,1), it creates ruptures
    return 0;
}
