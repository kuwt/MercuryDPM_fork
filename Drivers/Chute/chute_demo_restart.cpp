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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "Chute.h"
using namespace std;

class ChuteDemo : public Chute {
public:
    void actionsBeforeTimeStep(){
        random.setRandomSeed(0);
        Chute::actionsBeforeTimeStep();
    }
};

/// \bug When restarting the first time step is not saved, therefore there is a missing time step after a restart
void run_chute_demo(int argc, char *argv[])
{
    // Problem parameters
    ChuteDemo problem;
    problem.setName("chute_demo_restart");
    problem.setTimeMax(0.5);
    
    // Particle properties
    problem.setFixedParticleRadius(0);
    problem.setInflowParticleRadius(0.001);
    problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(2.5e-3, 0.8);
    problem.speciesHandler.getObject(0)->setSlidingDissipation(0.0);
    problem.setTimeStep(); 
    
    // Chute properties
    problem.setChuteAngle(30.0);
    problem.setXMax(0.1);
    problem.setYMax(2.0*problem.getInflowParticleRadius());
    
    // Inflow properties
    problem.setInflowHeight(0.1);
    problem.setInflowVelocity(0.1);
    problem.setInflowVelocityVariance(0.02);

    // Output properties
    problem.setSaveCount(100);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    
    //solve
    problem.solve(argc, argv);
    cout << problem.getTime() << endl;
}

/** This code tests if chute problems restart properly. The simulation itself is similar to the chute_demo.
 **/
int main(int argc, char *argv[])
{
    if (argc==1) 
    {
        //Print description
        cout << endl << "Description: This code tests if chute problems restart properly. The simulation itself is similar to the chute_demo." << endl;
        if(system("./chute_demo_restart.exe    -tmax 0.00999 -name chute_stop")){cout<<"Call failed"<<endl;}
        if(system("./chute_demo_restart.exe    -tmax 0.09999 -name chute_demo_no_restart")){cout<<"Call failed"<<endl;}
        if(system("./chute_demo_restart.exe    -tmax 0.00999")){cout<<"Call failed"<<endl;}
        if(system("./chute_demo_restart.exe -r -tmax 0.09999")){cout<<"Call failed"<<endl;}
    } else run_chute_demo(argc, argv);
    
    //to compare restarted with non-restarted data, use
    //gnuplot> var=3; plot 'chute_demo_restart/chute_demo_restart.ene' u 1:var w l, 'chute_demo_restart/chute_demo_no_restart.ene' u 1:var w l
    //../sc/fpdiff.py ./chute_demo_restart.fstat ./chute_demo_no_restart.fstat
}
