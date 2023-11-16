//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


// This module loads clump configuration produced by MClump tool

#ifndef CLUMP_INPUT_H
#define CLUMP_INPUT_H

#include <algorithm>
#include <dirent.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include "Math/Matrix.h"
#include <CMakeDefinitions.h>

// Useful containers
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> Double2DVector;
typedef std::vector<std::string> StringVector;

// Helper function for alphabetical sorting of clump names
bool CompareFunction (std::string a, std::string b) {return a < b;}

// Helper function to generate a random double in range (0, MAX)
double RandomDouble(double Max)
{return Max*((double) rand() / (RAND_MAX));}

/*!
* Structure for storing clump instance parameters
*/

struct ClumpData
{
public:
    
	std::string path;	// Path to MClump working directory
	StringVector clumpNames;	// Array of names of Clumps that will be used
	DoubleVector mass;		//  clump mass
	Double2DVector  pebblesX;	//  Pebbles geometry (outer index goes over Clumps, inner - over pebbles)
	Double2DVector  pebblesY;
	Double2DVector  pebblesZ;
	Double2DVector  pebblesR;
	
	Double2DVector  toi;		// Clump tensor of inertia (I11, I12, I13, I21..I33)
	Double2DVector  pd;		// Clump principal directions v1, v2, v3
	
};


/*!
* Loading available Clumps and their names
*/

void LoadConf(ClumpData &a)
{
    // Path to MCLump tool
    a.path = getMercuryDPMSourceDir() + "/Tools/MClump/Clumps/";

    struct dirent *entry;
    DIR *dir = opendir(a.path.c_str());
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_DIR){
            a.clumpNames.push_back(entry->d_name);
        }
    }
    closedir(dir);

    // remove "." and ".." from the list of dirs (position of those in the list is OS-sensitive, hence this code)

    std::sort(a.clumpNames.begin(), a.clumpNames.end(), CompareFunction);//sort the vector
    a.clumpNames.erase(a.clumpNames.begin());
    a.clumpNames.erase(a.clumpNames.begin());

    // Show the names of available Clumps
    for (int i = 0; i<a.clumpNames.size(); i++) std::cout << a.clumpNames[i] << std::endl;

}

/*!
* Loading pebbles of a clump
*/

void LoadPebbles(ClumpData &a)
{
    logger(INFO, "Loading clump pebbles...");

    a.pebblesX.resize(a.clumpNames.size());
    a.pebblesY.resize(a.clumpNames.size());
    a.pebblesZ.resize(a.clumpNames.size());
    a.pebblesR.resize(a.clumpNames.size());

    for (int i = 0; i < a.clumpNames.size(); i++ ){
        std::ifstream infile((a.path + a.clumpNames[i] + "/Clump/Clump.txt").c_str(), std::ios::in | std::ios::binary);
        std::string line{};
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string substring{};
            StringVector val;
            while (std::getline(iss, substring, ',')) val.push_back(substring);
            a.pebblesX[i].push_back(std::stof(val[0]));
            a.pebblesY[i].push_back(std::stof(val[1]));
            a.pebblesZ[i].push_back(std::stof(val[2]));
            a.pebblesR[i].push_back(std::stof(val[3]));
        }
        infile.close();
    }
    logger(INFO, "\t OK\n");
}

/*!
* Load mass of a clump
*/
void LoadMass(ClumpData &a)
{
    logger(INFO, "Loading clump masses...");

    a.mass.resize(a.clumpNames.size());
    
    for (int i = 0; i < a.clumpNames.size(); i++ ){
    	std::ifstream infile((a.path + a.clumpNames[i] + "/Inertia/Mass.txt").c_str(), std::ios::in | std::ios::binary);
    	std::string mass;
        infile >> mass;
        a.mass[i] = std::stof(mass);
        infile.close();
    }
    logger(INFO, "\t OK\n");
}

/*!
* Load tensor of inertia (TOI) of a clump
*/

void LoadTOI(ClumpData &a)
{
    logger(INFO, "Loading clump TOI..");
    a.toi.resize(a.clumpNames.size());

    for (int i = 0; i < a.clumpNames.size(); i++ ){
        std::ifstream infile((a.path + a.clumpNames[i] + "/Inertia/TOI.txt").c_str(), std::ios::in | std::ios::binary);
        std::string line{};
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string substring{};
            StringVector val;
            while (std::getline(iss, substring, ',')) val.push_back(substring);
            for (int k = 0; k < 3; k++) a.toi[i].push_back(std::stof(val[k]));
        }
        infile.close();
    }
    logger(INFO, "\t OK\n");
}

/*!
* Load pre-computed principal directions (PDs) of a clump. Normally MClump tool aligns PDs with the global Cartesian axes.
*/

void LoadPD(ClumpData &a)
{
    logger(INFO, "Loading clump PD..");
    a.pd.resize(a.clumpNames.size());

    for (int i = 0; i < a.clumpNames.size(); i++ ){
        std::ifstream infile((a.path + a.clumpNames[i] + "/Inertia/PD.txt").c_str(), std::ios::in | std::ios::binary);
        std::string line{};
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string substring{};
            StringVector val;
            while (std::getline(iss, substring, ',')) val.push_back(substring);
            for (int k = 0; k < 3; k++) a.pd[i].push_back(std::stof(val[k]));
        }
        infile.close();
    }
    logger(INFO, "\t OK\n");
}


/*!
*  // Main function that loads all the necessary files to initiate Clumps
*/

void LoadClumps(ClumpData &data, bool VERBOSE = false)
{

	logger(INFO, "LOAD CLUMP DATA");
    LoadConf(data);
    LoadPebbles(data);
    LoadMass(data);
    LoadTOI(data);
    LoadPD(data);
    if (VERBOSE) {
        std::cout<<"LOADED CLUMPS"<<std::endl;
        for (int i = 0; i < data.pebblesX.size(); i++) {
            std::cout << data.clumpNames[i] << " mass:" << data.mass[i] << std::endl;
            std::cout << data.clumpNames[i] << " list of pebbles:" << std::endl;
            for (int j = 0; j < data.pebblesX[i].size(); j++) {
                std::cout << "Pebble " << j << ": (" << data.pebblesX[i][j] << "," << data.pebblesY[i][j] << ","
                          << data.pebblesZ[i][j] << ")," << data.pebblesR[i][j] << std::endl;
            }

            std::cout << data.clumpNames[i] << " TOI:" << std::endl;
            std::cout << data.toi[i][0] << "," << data.toi[i][1] << "," << data.toi[i][2] << std::endl;
            std::cout << data.toi[i][3] << "," << data.toi[i][4] << "," << data.toi[i][5] << std::endl;
            std::cout << data.toi[i][6] << "," << data.toi[i][7] << "," << data.toi[i][8] << std::endl;

            std::cout << data.clumpNames[i] << " Principal directions:" << std::endl;
            std::cout << data.pd[i][0] << "," << data.pd[i][1] << "," << data.pd[i][2] << std::endl;
            std::cout << data.pd[i][3] << "," << data.pd[i][4] << "," << data.pd[i][5] << std::endl;
            std::cout << data.pd[i][6] << "," << data.pd[i][7] << "," << data.pd[i][8] << std::endl;
        }
    }

}

/*!
*  Changes PDs of the clump - sufficient for clump to be properly rotated
*/

ClumpData RotateClump(ClumpData data, int clump_index, DoubleVector new_pd)
{
    ClumpData new_data = data;
    new_data.pd[clump_index] = new_pd;
    return new_data;
}






/*!
*  Generate random (and isotropically distributed) principal directions frame
*/
DoubleVector UniformRandomPDs(){

    Vec3D n1, n2, n3, ref;

    // basis vector n1
    double r1 = RandomDouble(2) - 1.0;
    double r2 = RandomDouble(1);

    double theta = acos(r1); // Note that for isotropy of n1 theta is NOT uniform!
    double phi = 2 * M_PI * r2;

    n1 = Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

    // basis vector n2
    r1 = RandomDouble(1);
    r2 = RandomDouble(1);

    theta = acos(2 * r1 - 1); // Note that for isotropy of n1 theta is NOT uniform!
    phi = 2 * M_PI * r2;

    ref = Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    n2 = Vec3D::cross(ref, n1); n2.normalise();
    n3 = Vec3D::cross(n1, n2); n3.normalise();
    return DoubleVector{n1.X, n1.Y, n1.Z, n2.X, n2.Y, n2.Z, n3.X, n3.Y, n3.Z};

}


#endif // CLUMP_INPUT_H
