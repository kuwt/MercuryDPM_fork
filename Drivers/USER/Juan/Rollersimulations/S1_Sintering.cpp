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

/* This program is created to apply sintering over particles and measure the neck growth.
 *
 */

#include "Mercury3D.h"

#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>

class Sintering : public Mercury3D{
public:

    Sintering()
    {
        setName("S0_InitialConditions");
        readRestartFile();

        wallHandler.removeObject(6);

        setRestarted(false);
        setName("S1_Sintering");

        particleNewSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));

        setGravity({0,0,-9.81e-9}); //[mm/micros^2]//To make the units consistent.

    }

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleNewSpecies; //pointer to the particle properties

    void setupInitialConditions() override
    {

        Mdouble HeatCapacity = 500;
        Mdouble ThermalCond = 20;
        Mdouble sinterTime = 0.1;//readFromFile("in","sinterTime",24);


        //particleNewSpecies->setHeatCapacity(HeatCapacity);
        //particleNewSpecies->setThermalConductivity(ThermalCond);
        //particleNewSpecies->setSinterRate(sinterTime);
        //particleNewSpecies->setSinterType(SINTERTYPE::CONSTANT_RATE);

    }

    void actionsAfterTimeStep() override
    {
        //-------->[Set Adhesive Contribution]
        //Set Adhesive properties:
        Mdouble K_adh = particleNewSpecies->getLoadingStiffness();
        Mdouble f_adh_max = 1e-5*K_adh;//5.0*(meanRadius*2.0*meanRadius*2.0)*5.0; ///set adhesive force max based on Bond number [Hao]

        particleNewSpecies->setAdhesionStiffness(K_adh);
        particleNewSpecies->setAdhesionForceMax(f_adh_max);

        if(getKineticEnergy()/getElasticEnergy() < 1e-4)//I will wait again, until particles reach the equilibrium.
        {
            particleNewSpecies->setSinterType(SINTER_APPROACH::FRENKEL);
            particleNewSpecies->setSinterRate(4e-10/0.002);
            particleNewSpecies->setSinterAdhesion(0.013*particleNewSpecies->getLoadingStiffness());

        }


    }

    void printTime() const override
    {
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << "  EneRatio " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << "  KineticEnergy " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                  << "  ElasticEnergy " << std::setprecision(3) << std::left << std::setw(6) << getElasticEnergy()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Sintering problem;

    problem.setTimeMax(0.1);
    problem.solve();

    return 0;
}