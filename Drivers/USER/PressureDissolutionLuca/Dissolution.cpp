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

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class Dissolution : public Mercury3D{
public:

    //set default values
    Dissolution()
    {
        setName("Precompression");
        readRestartFile();
        setRestarted(false);
        setName("Dissolution");
        lid = dynamic_cast<InfiniteWall*>(wallHandler.getLastObject());
        if (lid==nullptr) {
            std::cerr << "Error: lid should be a infinite wall" << std::endl;
        }
    }

    //set slave variables (i.e. compute stiffness, dissipation, time step, wall and particle positions)
    void setupInitialConditions() override
	{
    }

    void actionsBeforeTimeStep() override
    {
        for (auto p : particleHandler)
        {
            if (p->getRadius()<limitRadius) {
                //\todo
                //particleHandler.removeObject(p->getIndex());
            } else if (p->getIndSpecies()==0) {
                p->setRadius(p->getRadius()*(1.0-capRockDissolutionRate*getTimeStep()));
            } else if (p->getIndSpecies()==1) {
                p->setRadius(p->getRadius()*(1.0-pureSilicaDissolutionRate*getTimeStep()));
            } else if (p->getIndSpecies()==2) {
                p->setRadius(p->getRadius()*(1.0-pureCalciumCarbonateDissolutionRate*getTimeStep()));
            }
        }
    }

    //override continueSolve function such that the code stops 
    //when the packing is relaxed (Ekin<1e-5*Eela) and 
    //the pressure on the lid is correct (|p/p_lid-1|<1e-3)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;

            static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
            static Mdouble particleArea = constants::pi * mathsFunc::square(particleHandler.getObject(0)->getRadius());
            static Mdouble stiffness = dynamic_cast<const LinearPlasticViscoelasticSpecies*>(speciesHandler.getObject(0))->getLoadingStiffness();
            // amount by which the pressure has to be increased
            Mdouble dPressure = lid->getForce().Z/lidArea - pressure;
            // amount by which position should be changed to achieve the right pressure
            Mdouble dV = dPressure * particleArea / stiffness /getTimeStep();
            //std::cout << "dP/P" << dPressure/pressure << " Z" << dZ << std::endl;
            lid->setVelocity(Vec3D(0.0,0.0,dV/50.0));

            if (std::abs(dPressure)<1e-3*pressure && getKineticEnergy()<1e-5*getElasticEnergy())
            {
                printTime();
                return false;
            }
        }
        return true;
    }

    void printTime() const override
    {
        static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
        std::cout << "t=" << getTime()
            << " Ene " << getKineticEnergy()/getElasticEnergy()
            << " |dP/P| " << 1.0-(lid->getForce().Z/lidArea)/pressure
            << " lidZ " << lid->getPosition().Z
            << std::endl;
    }

    InfiniteWall* lid;

    Mdouble pressure;

    Mdouble capRockDissolutionRate; //relative dissolution rate
    Mdouble pureSilicaDissolutionRate; //
    Mdouble pureCalciumCarbonateDissolutionRate; //

    Mdouble limitRadius;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Dissolution d;
    d.capRockDissolutionRate = 0; //
    d.pureSilicaDissolutionRate = 1; //
    d.pureCalciumCarbonateDissolutionRate = 1; //
    d.limitRadius = 1e-4;
    d.pressure = 5e3; //pressure on the lid
    d.setXBallsAdditionalArguments(" -v0 -solidf ");
    d.solve();
}
