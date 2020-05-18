//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPM.org/Team>.
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

/** \brief Test of the binary reader. The files used is STL file with containing one triange
\details This code users BinaryReader to open a binary file called SimpleTrianlge.stl and read in currently the single tringle from this file. The triangle was created using nclab.org and contains a single triangle with normal (0,0,1) and vertex (0,1,0), (1,0,0) and (1,1,0).
 This is currently a UnitTest as does not require selftest data; however, it does require he file SimpleTrianlge.stl to exist
**/

#include<fstream>
#include<iostream>
#include<BinaryReader.h>
#include<Logger.h>
#include<Math/Vector.h>
#include<CMakeDefinitions.h>





int main()
{
   
    std::ifstream STLFile;
    
    std::string directory(getMercurySourceDir());

    BinaryReader STLReader(directory+"/Drivers/ImportTools/ExampleSTLFiles/SimpleTrianlge.stl");
    
    //First read the 80 character header
    std::string header;
    header=STLReader.readString(80);
    //logger(INFO, "Header : %" ,header);
   
    //The next four characters contain at unsigned int which is the number of triangeles
    unsigned int numTriangles = STLReader.readUnsignedInt(4);
    //logger(INFO, "Number of traingles: %", numTriangles);
    
    if (numTriangles!=1)
        logger(FATAL,"Failed to read the correct number of triangles");

    double xTmp,yTmp,zTmp;
    
    for (unsigned int i=0; i<(numTriangles);i++)
    {
        
       
        xTmp = STLReader.readFloat(4);
        yTmp =  STLReader.readFloat(4);
        zTmp = STLReader.readFloat(4);
        
        Vec3D Normal(xTmp,yTmp,zTmp);
        
        if (!(Normal.isEqualTo(Vec3D(0,0,1),1e-10)))
            logger(FATAL,"The normal has been misread");
        
        
        xTmp = STLReader.readFloat(4);
        yTmp =  STLReader.readFloat(4);
        zTmp = STLReader.readFloat(4);
        
        Vec3D Point1(xTmp,yTmp,zTmp);
        
        if (!(Point1.isEqualTo(Vec3D(0,1,0),1e-10)))
            logger(FATAL,"The first vertex has been misread");
        
        xTmp = STLReader.readFloat(4);
        yTmp =  STLReader.readFloat(4);
        zTmp = STLReader.readFloat(4);
        
        Vec3D Point2(xTmp,yTmp,zTmp);
        
        if (!(Point2.isEqualTo(Vec3D(1,0,0),1e-10)))
            logger(FATAL,"The second vertex has been misread");
     
        
        xTmp = STLReader.readFloat(4);
        yTmp =  STLReader.readFloat(4);
        zTmp = STLReader.readFloat(4);
        
         Vec3D Point3(xTmp,yTmp,zTmp);
        
        if (!(Point3.isEqualTo(Vec3D(1,1,0),1e-10)))
            logger(FATAL,"The third vertex has been misread");
        
        //Now ignore (read) the two dummy characters
        STLReader.ignoreChar(2);
    
    }
    
    
}
