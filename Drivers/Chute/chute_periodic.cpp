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
#include <sys/types.h>
#include <sys/stat.h>
#include "Chute.h"
using namespace std;

class ChutePeriodic : public Chute{
public:
    ///make chute periodic in x and y and calls #add_flow_particles.
    void setupInitialConditions() {
    
        Chute::setupInitialConditions();

        PeriodicBoundary b0;
        b0.set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        
        add_flow_particles();
    }

    ///Do not add or remove particles during the run
    void actionsBeforeTimeStep(){ };
        
    ///Add initial flow particles.
    void add_flow_particles() 
    {
        cout << "Adding flowing particles" << endl;
        //set_HGRID_num_buckets_to_power(particleHandler.getStorageCapacity());
        HGridActionsBeforeTimeLoop();
        HGridActionsBeforeTimeStep();
        writeRestartFile();
        unsigned int N=particleHandler.getNumberOfObjects()+getChuteLength()*getChuteWidth()*InflowHeight/mathsFunc::cubic(getInflowParticleRadius())/8;
        particleHandler.setStorageCapacity((N));
        double H = InflowHeight;
        setZMax(1.2*InflowHeight);
        writeRestartFile();
        //try to find new insertable particles
        while (particleHandler.getNumberOfObjects()<N){
            create_inflow_particle();
            if (IsInsertable(P0)) {
                num_created++;
            } else {
                InflowHeight += .0001* MaxInflowParticleRadius;
            }
        }
        cout << "InflowHeight=" << InflowHeight << endl;
        InflowHeight = H;
        set_HGRID_num_buckets_to_power();
    }
    
    ///creates flow particles in the whole chute
    void create_inflow_particle()
    {
        P0.setRadius(random.get_RN(MinInflowParticleRadius,MaxInflowParticleRadius));
        P0.computeMass(Species);
        P0.setPosition(Vec3D(random.get_RN(getXMin()+2.0*P0.getRadius(),getXMax()),
                              random.get_RN(getYMin()+2.0*P0.getRadius(),getYMax()),
                              random.get_RN(getZMin()+2.0*P0.getRadius(),getInflowHeight())));
        P0.setVelocity(Vec3D(0.0,0.0,0.0));
    }
    
    //~ void printTime() const {
        //~ cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
            //~ << ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax()
            //~ << ", N=" << setprecision(3) << left << setw(6) << particleHandler.getNumberOfObjects()
            //~ << ". " << endl;
    //~ }

};

/** A rough inclined plane with periodic boundary conditions.
 **/
int main(int argc, char *argv[])
{
    //Print description
    cout << endl << "Description: A rough inclined plane with periodic boundary conditions." << endl;

    // Problem parameters
    ChutePeriodic problem;
    problem.setName("chute_periodic");
    problem.setTimeMax(.1);
 
    // Particle properties
    problem.setDensity(2400.0);
    problem.setInflowParticleRadius(0.5e-3);
    problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(4e-4, 0.8);
    problem.speciesHandler.getObject(0)->setSlidingDissipation(problem.get_dissipation());
    problem.speciesHandler.getObject(0)->setSlidingStiffness(2./7.*problem.speciesHandler.getObject(0)->getStiffness());
    problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.5);
    problem.setFixedParticleRadius(problem.getInflowParticleRadius());
    problem.setRoughBottomType(MULTILAYER);
    problem.setTimeStep(); //sets time step to 1/50th of the collision time
    cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; 
    
    // Chute properties
    problem.setChuteAngle(25.0);
    problem.setChuteLength(5e-3);
    problem.setChuteWidth(3e-3);
    problem.setInflowHeight(5e-3);
    problem.makeChutePeriodic();
    
    // Output properties
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1,problem.getTimeMax(),problem.getTimeStep())); //minimize output to the last time step
    problem.set_number_of_saves_data(100); //allow enough data output so the evolution can be viewed in xballs
    problem.set_number_of_saves_ene(100);
    problem.setXBallsAdditionalArguments("-sort -v0 -solidf");
    
    //solve
    problem.solve(argc,argv);
}
