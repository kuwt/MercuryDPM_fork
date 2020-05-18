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
#include <Species/LinearPlasticViscoelasticSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <iomanip>

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class Compression : public Mercury3D{
public:

    //set default values
    Compression()
    {
        setName("Sintering");
        readRestartFile();
        setRestarted(false);
        //fStatFile.setFileType(FileType::MULTIPLE_FILES);
        setName("Compression");
        lid = dynamic_cast<InfiniteWall*>(wallHandler.getLastObject());
    }

    void setupInitialConditions() override
	{
    }

    Mdouble getPressure() const
    {
        Mdouble pressureGradient = 2.0* (maxPressure-initialPressure)/getTimeMax();
        Mdouble increasingPressure = initialPressure + pressureGradient*getTime();
        Mdouble decreasingPressure = initialPressure + pressureGradient*(getTimeMax() - getTime());
        return std::min(increasingPressure,decreasingPressure);
    }

    void actionsAfterTimeStep() override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;

            static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
            static Mdouble particleArea = constants::pi * mathsFunc::square(particleHandler.getObject(0)->getRadius());
            static Mdouble stiffness = dynamic_cast<const LinearPlasticViscoelasticSpecies*>(speciesHandler.getObject(0))->getLoadingStiffness();
            // amount by which the pressure has to be increased
            Mdouble dPressure = lid->getForce().Z/lidArea - getPressure();
            // amount by which position should be changed to achieve the right pressure
            Mdouble dZ = dPressure * particleArea / stiffness;
            lid->setVelocity(Vec3D(0.0,0.0,dZ * 100.0));
        }
    }


    void printTime() const override
    {
        static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
        std::cout << "t=" << std::left << std::setw(4) << getTime()
            << " Pressure " << std::left << std::setw(6) << getPressure()
            << " Lid " << std::left << std::setw(6) << lid->getPosition().Z
            << std::endl;
    }

    InfiniteWall* lid;
    Mdouble initialPressure;
    Mdouble maxPressure;

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Compression s;
    s.initialPressure = 200;
    s.maxPressure = 2000.0;
    s.setTimeMax(10.0);
    s.solve();
}
