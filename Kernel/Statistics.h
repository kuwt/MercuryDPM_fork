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

#ifndef STATISTICS_H
#define STATISTICS_H

#include "StatisticsVector.h"

/*!
 * \brief This is the function that the user should call for statistics (CG).
 * \details 
 * This function offers an interface to the templated class StatisticsVector. It
 * picks the correct StatType template from the user-specified flag. The other
 * flags are processed by StatisticsVector.
 */
void Statistics(int argc, char* argv[])
{
    if (argc > 1 && strcmp(argv[1], "-help"))
        logger(INFO, "\nGet statistics for %", argv[1]);
    
    /* Check for the '-stattype' flag */
    StatType T = O; //default value
    for (unsigned int i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-stattype") || !strcmp(argv[i], "-statType"))
        {
            if (!strcmp(argv[i + 1], "XYZ"))
                T = XYZ;
            else if (!strcmp(argv[i + 1], "RAZ"))
                T = RAZ;
            else if (!strcmp(argv[i + 1], "RA"))
                T = RA;
            else if (!strcmp(argv[i + 1], "RZ"))
                T = RZ;
            else if (!strcmp(argv[i + 1], "AZ"))
                T = AZ;
            else if (!strcmp(argv[i + 1], "R"))
                T = R;
            else if (!strcmp(argv[i + 1], "A"))
                T = A;
            else if (!strcmp(argv[i + 1], "XY"))
                T = XY;
            else if (!strcmp(argv[i + 1], "XZ"))
                T = XZ;
            else if (!strcmp(argv[i + 1], "YZ"))
                T = YZ;
            else if (!strcmp(argv[i + 1], "X"))
                T = X;
            else if (!strcmp(argv[i + 1], "Y"))
                T = Y;
            else if (!strcmp(argv[i + 1], "Z"))
                T = Z;
            else if (!strcmp(argv[i + 1], "O"))
                T = O;
            else
            {
                logger(ERROR, "stattype unknown");
            }
        }
    }
    if (T == XY)
    { // averaging in z-direction
        logger(INFO, "averaging in z-direction");
        StatisticsVector<XY> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == XZ)
    { // averaging in y-direction
        logger(INFO, "averaging in y-direction");
        StatisticsVector<XZ> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == YZ)
    { // averaging in x-direction
        logger(INFO, "averaging in x-direction");
        StatisticsVector<YZ> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == X)
    { // averaging in yz-direction
        logger(INFO, "averaging in yz-direction");
        StatisticsVector<X> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == Y)
    { // averaging in yz-direction
        logger(INFO, "averaging in xz-direction");
        StatisticsVector<Y> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == Z)
    { // averaging in yz-direction
        logger(INFO, "averaging in xy-direction");
        StatisticsVector<Z> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == O)
    { // averaging in all directions
        logger(INFO, "averaging in xyz-direction");
        StatisticsVector<O> stats(argc, argv);
        stats.setDoPeriodicWalls(false);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == RAZ)
    { //no averaging
        logger(INFO, "cylindrical, no averaging");
        StatisticsVector<RAZ> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == RA)
    { //no averaging
        logger(INFO, "cylindrical, Z averaging");
        StatisticsVector<RA> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == RZ)
    { //no averaging
        logger(INFO, "cylindrical, A averaging");
        StatisticsVector<RZ> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == AZ)
    { //no averaging
        logger(INFO, "cylindrical, R averaging");
        StatisticsVector<AZ> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == A)
    { //no averaging
        logger(INFO, "cylindrical, RZ averaging");
        StatisticsVector<A> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == R)
    { //no averaging
        logger(INFO, "cylindrical, AZ averaging");
        StatisticsVector<R> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    else if (T == XYZ)
    { //no averaging
        logger(INFO, "no spatial averaging");
        StatisticsVector<XYZ> stats(argc, argv);
        stats.statistics_from_fstat_and_data();
    }
    
}

#endif

