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

#include "MultiDrum.h"

int main(int argc, char *argv[])
{
    //Print description
    std::cout << std::endl << "Description: A quasi-2D moving-bed channel with walls on the left and right boundary." << std::endl;
    
    // Problem parameters
    RotatingDrum problem;
    problem.setName("drum");
//    problem.autoNumber();
    problem.setTimeMax(2000.0);
    problem.setTimeStep(1./(200.0* 50));
    
    
    problem.set_radiusLarge(0.75);
    problem.setDensityRatio(1.0);
    problem.setCoefficientOfRestitution(0.8);
    problem.set_particle_number_volRatio(1.0); //volume ratio of large to small
    //problem.set_particle_numbers(16700, 16700);

//    problem.set_particle_numbers(100,100);
    problem.set_particle_numbers(20,20);


    problem.setChuteAngleAndMagnitudeOfGravity(0.0, 1.0);
    problem.setRPM(0.01);
    
    
    // Chute properties : Simply remove the first line to add side walls.
   //problem.makeChutePeriodic();
    //problem.setYMax(40.5);
    problem.setYMax(10);
    //problem.setDrumRadius(35.5);
    problem.setDrumRadius(10.0);
    
    
    // scalefactor 2 is fine for n-polygon but for stars it most be higher. 3.0 was not higher enought so when to 4;
    /// \bug The scale factor set here gets ignore and the only used one is hard wired in the code. Just track this bug down later.
    problem.setWallParameters(1.0,3.0);

    
    
    //Swap the next two lines to swap between the different type of rought bottoms. Note MONOLAYER_ORDERED does not give enough friction for the circular case.
    //problem.setRoughBottomType(MULTILAYER);
    problem.setRoughBottomType(MONOLAYER_DISORDERED);
    //problem.setRoughBottomType(MONOLAYER_ORDERED);

    
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000, problem.getTimeMax(), problem.getTimeStep()));
    //Colour of mass
    //problem.setXBallsAdditionalArguments("-cmode 8");
    //Color on size
    problem.setXBallsAdditionalArguments("-cmode 7");
    problem.readArguments(argc, argv);
//    problem.setStarShape(5,1,5,2);
    problem.addSegment(5,1, 3);
    problem.addSegment(5,2, 4);
    problem.addSegment(5,1, 3);
//    problem.addSegment(5,2);
//    problem.addSegment(7,1);
//    problem.addSegment(7,4);
//    problem.addSegment(7,3);


    
    if (argc > 4)
    {
        problem.num_restart_large=atoi(argv[3]);
        problem.num_restart_small=atoi(argv[4]);
    }
    else
    {
        problem.num_restart_large=0;
        problem.num_restart_small=0;
    }
    
    std::cout << problem.num_restart_small <<" " <<problem.num_restart_large << std::endl;

    problem.solve();
    
}
