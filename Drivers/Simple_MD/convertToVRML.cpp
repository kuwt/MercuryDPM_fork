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

#include <iostream>
#include "DPMBase.h"


class Data : public DPMBase {
public:

	void CylinderTransformer(Vec3D pointA, Vec3D pointB, double* length, double* angle, Vec3D* translation, Vec3D* axis)
	{

		*length = GetDistance(pointA, pointB);
			
		// the translation (center of the cylinder) needed is
		*translation = (pointA + pointB) / 2.0f;

		//the initial orientation of the bond is in the y axis
		Vec3D init = Vec3D(0.0f, 1.0f, 0.0f);
			
		//the vector needed is the same as that from a to b
		Vec3D needed = Vec3D(pointB.X - pointA.X, pointB.Y - pointA.Y, pointB.Z - pointA.Z);
			

		Vec3D needed_n, init_n; // normilized

		needed_n	= needed / needed.getLength();
		init_n		= init / init.getLength();

		//so the angle to rotate the bond by is:
		*angle = acos(Vec3D::dot(needed_n, init_n));

		
		//and the axis to rotate by is the cross product of the initial and
		//needed vectors - ie the vector orthogonal to them both
		double vx = init.Y * needed.Z - init.Z * needed.Y;
		double vy = init.Z * needed.X - init.X * needed.Z;
		double vz = init.X * needed.Y - init.Y * needed.X;

		*axis = Vec3D(vx, vy, vz);

	}

	// Export packing configuration to the VRML format.
	void Export_to_VRML_format(string fname)
	{
		ofstream grout(fname.c_str());	

		grout << "\
#VRML V2.0 utf8 \n\n\
Group {\n\
  children [\n\
	WorldInfo {\n\
	  title \"Granular Media\"\n\
	  info [\n\
		\"Packing configuration.\"\n\
		\"Author: Vitaliy Ogarko; vogarko@gmail.com\"\n\
	  ]\n\
	}\n\
	NavigationInfo {\n\
	  type [ \"EXAMINE\", \"ANY\" ]\n\
	}\n\
	Background {\n\
	  skyColor [ 1 1 1 ]\n\
	}\n\
  ]\n\
}\n\n";

/*
	grout << "\
Shape {\n\
	appearance Appearance {\n\
		material Material {\n\
			emissiveColor 0 0 0\n\
		}\n\
	}\n\
	geometry IndexedLineSet {\n\
		coord Coordinate {\n\
			point [\n\
			# Coordinates around the top of the cube\n\
				   0  1.0  1.0,\n\
				 1.0  1.0  1.0,\n\
				 1.0  1.0    0,\n\
				   0  1.0    0,\n\
			# Coordinates around the bottom of the cube\n\
				   0    0  1.0,\n\
				 1.0    0  1.0,\n\
				 1.0    0    0,\n\
				   0    0    0\n\
			]\n\
		}\n\
		coordIndex [\n\
		# top\n\
			0, 1, 2, 3, 0, -1,\n\
		# bottom\n\
			4, 5, 6, 7, 4, -1,\n\
		# vertical edges\n\
			0, 4, -1,\n\
			1, 5, -1,\n\
			2, 6, -1,\n\
			3, 7\n\
		]\n\
	}\n\
}\n\n";
*/


	////// (+) DRAW CUBE -----------------------------------

		Vec3D translation, axis;
		double length, angle;

		Vec3D pointsA[12];
		Vec3D pointsB[12];

		pointsA[0] = Vec3D(0,0,0); pointsB[0] = Vec3D(1,0,0);
		pointsA[1] = Vec3D(0,0,0); pointsB[1] = Vec3D(0,1,0);
		pointsA[2] = Vec3D(0,0,0); pointsB[2] = Vec3D(0,0,1);

		pointsA[3] = Vec3D(1,0,0); pointsB[3] = Vec3D(1,0,1);
		pointsA[4] = Vec3D(1,0,0); pointsB[4] = Vec3D(1,1,0);

		pointsA[5] = Vec3D(0,1,0); pointsB[5] = Vec3D(0,1,1);
		pointsA[6] = Vec3D(0,1,0); pointsB[6] = Vec3D(1,1,0);

		pointsA[7] = Vec3D(0,0,1); pointsB[7] = Vec3D(0,1,1);
		pointsA[8] = Vec3D(0,0,1); pointsB[8] = Vec3D(1,0,1);

		pointsA[9]  = Vec3D(1,1,1); pointsB[9]  = Vec3D(0,1,1);
		pointsA[10] = Vec3D(1,1,1); pointsB[10] = Vec3D(1,1,0);
		pointsA[11] = Vec3D(1,1,1); pointsB[11] = Vec3D(1,0,1);

		double domainsize = max(std::max(getXMax(),getYMax()),getZMax())/4;

		for(int i=0; i<12; i++) {
			pointsA[i].X*=getXMax()/domainsize;
			pointsA[i].Y*=getYMax()/domainsize;
			pointsA[i].Z*=getZMax()/domainsize;
			pointsB[i].X*=getXMax()/domainsize;
			pointsB[i].Y*=getYMax()/domainsize;
			pointsB[i].Z*=getZMax()/domainsize;
		}


		for(int i=0; i<12; i++)
		{

			CylinderTransformer(pointsA[i], pointsB[i], &length, &angle, &translation, &axis);

			if (i == 11) { angle = 0.; } // without this a stick becomes black colored

grout << "\
Transform{\n\
translation " << translation.X << " " << translation.Y << " " << translation.Z << "\n\
rotation " << axis.X << " " << axis.Y << " " << axis.Z << " " << angle << "\n\
children  Shape {\n\
appearance Appearance {\n\
material Material { diffuseColor 0.8 0.8 0.8 }\n\
}\n\
geometry Cylinder {\n\
radius 0.005\n\
height " << length << "\n\
}}}\n\n";


} 

// http://vrmlworks.crispen.org/orient.html >> calculate oriantation for a viewpoint
//~ Camera .5 -.4 1
//~ Look at 3 .2 0
//~ -g: 0 0 1


grout << "\
Viewpoint {\n\
position    -.5 -.4 1\n\
orientation 0.5822204185334514 -0.4590207210231789 -0.6710583893478226 1.7312487661554772\n\
fieldOfView 0.785398\n\
description \"Viewpoint 1\"\n\
}\n\n";

		double maxValue=0;
		for(int i=0; i<get_Nmax(); i++) 
			maxValue = max(maxValue,particleHandler.getObject(i)->getPosition().Z);

		
		for(int i=0; i<get_Nmax(); i+=1) if (particleHandler.getObject(i)->getPosition().X>0.2)//get_Nmax()
		{
	
			grout << "Transform {" << "\n";
			grout << "translation ";
			
			grout << particleHandler.getObject(i)->Position/domainsize << "\n";

			grout << "children Shape {" << "\n";
			grout << "geometry Sphere {" << "\n";
			grout << "radius ";

			grout << particleHandler.getObject(i)->getRadius()/domainsize*2 << "\n";

			grout << "} \n";


			grout << "appearance Appearance {" << "\n";

			//grout << "material Material { diffuseColor 0.915064 0 0.0849365 }" << "\n";
			grout << "material Material { diffuseColor ";
			
			// gets the jet colors for a number between 0 and 8
			double value=particleHandler.getObject(i)->getPosition().Z/maxValue*8;
			Vec3D rgb;
			if (value<1) {
				rgb=Vec3D(0,0,.5+value/2);
			} else if (value<3) {
				rgb=Vec3D(0,-.5+value/2,1);
			} else if (value<5) {
				rgb=Vec3D(-1.5+value/2,1,2.5-value/2);
			} else if (value<7) {
				rgb=Vec3D(1,3.5-value/2,0);
			} else {
				rgb=Vec3D(4.5-value/2,0,0);
			}
			
			grout << rgb << "\n";

			double transparency = 0.0;
			grout << "	transparency " << transparency << "\n}\n";

			grout << "} \n";
			grout << "} \n";
			grout << "} \n";
			grout << "\n";

		} 

		grout.close();

	}


	void load_data_file(const char* datafile) {
		if (!readDataFile(datafile,14)) {
			std::cerr << "input data not found exiting " << std::endl;
			exit(-1);
		}
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Data data;
	//~ data.load_data_file("../../../Chute/run/ini_files/H20A24L0.5M0.5Quick.ini");
	//~ string fname("out.wrl");
	data.load_data_file("../../../Chute/run/Vreman2/Vreman.data.0895");
	string fname("chute.wrl");
	data.Export_to_VRML_format(fname);
	std::cout << "Hello" << std::endl;
}
