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

#include <iostream>  // ostream
#include <iomanip>   // setw
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "StatisticsVector.h"

///This code (p3dstatistics.cpp) can be used to create statistics from p3d files. 
///It is a modified version of fstatistics.cpp
int main(int argc, char *argv[]) {	

    // if no argument is provided
	if (argc==1) {
        logger(ERROR,"No arguments; use p3statistics.exe filename [-options]");
    }
    else {
        // if you don't use the -help argument and if at least one argument is given
        if (strcmp(argv[1],"-help") != 0)
            logger(INFO,"Getting statistics from %.p3* (or .p4*) files ...",argv[1]);
    }



	//check for '-stattype' option
	StatType T = XYZ;
	for (int i = 2; i<argc; i++) {
		if (!strcmp(argv[i],"-stattype")) {
			if (!strcmp(argv[i+1],"XYZ")) T = XYZ;
			else if (!strcmp(argv[i+1],"RAZ")) T = RAZ;
			else if (!strcmp(argv[i+1],"RA")) T = RA;
			else if (!strcmp(argv[i+1],"RZ")) T = RZ;
			else if (!strcmp(argv[i+1],"AZ")) T = AZ;
			else if (!strcmp(argv[i+1],"R")) T = R;
			else if (!strcmp(argv[i+1],"A")) T = A;
			else if (!strcmp(argv[i+1],"YZ")) T = YZ;
			else if (!strcmp(argv[i+1],"XZ")) T = XZ;
			else if (!strcmp(argv[i+1],"XY")) T = XY;
			else if (!strcmp(argv[i+1],"X")) T = X;
			else if (!strcmp(argv[i+1],"Y")) T = Y;
			else if (!strcmp(argv[i+1],"Z")) T = Z;
			else if (!strcmp(argv[i+1],"O")) T = O;
			else {
                logger(ERROR,"Stattype unknown");
            }
		}
	}
	logger(INFO,"stattype %",T);

	if (T==XY) { // averaging in z-direction
 		StatisticsVector<XY> stats(argc, argv);
	} else if (T==XZ) { // averaging in y-direction
		StatisticsVector<XZ> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.statistics_from_p3();
 	} else if (T==YZ) { // averaging in x-direction
 		StatisticsVector<YZ> stats(argc, argv);
 		stats.setDoPeriodicWalls(false);
 		stats.statistics_from_p3();
 	} else if (T==X) { // averaging in yz-direction
 		StatisticsVector<X> stats(argc, argv);
 		stats.setDoPeriodicWalls(false);
 		stats.statistics_from_p3();
	} else if (T==Y) { // averaging in yz-direction
 		StatisticsVector<Y> stats(argc, argv);
 		stats.setDoPeriodicWalls(false);
 		stats.statistics_from_p3();
	} else if (T==Z) { // averaging in yz-direction
		StatisticsVector<Z> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.statistics_from_p3();
 	} else if (T==O) { // averaging in all directions
 		StatisticsVector<O> stats(argc, argv);
 		stats.setDoPeriodicWalls(false);
 		stats.statistics_from_p3();
	} else if (T==XYZ) { //no averaging
		StatisticsVector<XYZ> stats(argc, argv);
		stats.statistics_from_p3();
	} else {
        logger(ERROR,"stattype not found");
    }

	return 0;
}
	
