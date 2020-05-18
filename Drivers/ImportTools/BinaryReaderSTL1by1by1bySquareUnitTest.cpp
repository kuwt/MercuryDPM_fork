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

/** \brief Test of the STL reader. The files used is STL file with containing 12 triange that a 1 by 1 by 1 square and was created in autocad
**/

#include<fstream>
#include<iostream>
#include<vector>
#include <string>
#include<BinaryReader.h>
#include<Logger.h>
#include<Math/Vector.h>
#include<CMakeDefinitions.h>

class STLTriangle
{
public:
    STLTriangle(){};
    
    STLTriangle(const Vec3D newNormal, const Vec3D newVertex1, const Vec3D newVertex2, const Vec3D newVertex3)
    {
        normal=newNormal;
        vertex1=newVertex1;
        vertex2=newVertex2;
        vertex3=newVertex3;
    }
    
    bool isEqualTo(const STLTriangle answer, double toll)
    {
        if (!normal.isEqualTo(answer.normal,toll)) return false;
        if (!vertex1.isEqualTo(answer.vertex1,toll)) return false;
        if (!vertex2.isEqualTo(answer.vertex2,toll)) return false;
        if (!vertex3.isEqualTo(answer.vertex3,toll)) return false;
        return true;
    }
    
    Vec3D normal;
    Vec3D vertex1;
    Vec3D vertex2;
    Vec3D vertex3;
    
};


class STLReader : public BinaryReader
{
public:
    
    STLReader(std::string fileName) : BinaryReader(fileName)
    {
        
        header_= readString(80);
        numTriangles_ = readUnsignedInt(4);
        
        Triangles_.resize(numTriangles_);
        
        
        for (unsigned int i=0; i<(numTriangles_);i++)
        {
            
            
            Triangles_[i].normal.x() = readFloat(4);
            Triangles_[i].normal.y() = readFloat(4);
            Triangles_[i].normal.z() = readFloat(4);
            
        
            Triangles_[i].vertex1.x() = readFloat(4);
            Triangles_[i].vertex1.y() = readFloat(4);
            Triangles_[i].vertex1.z() = readFloat(4);
            
            Triangles_[i].vertex2.x() = readFloat(4);
            Triangles_[i].vertex2.y() = readFloat(4);
            Triangles_[i].vertex2.z() = readFloat(4);
            
            
            Triangles_[i].vertex3.x() = readFloat(4);
            Triangles_[i].vertex3.y() = readFloat(4);
            Triangles_[i].vertex3.z() = readFloat(4);
        
            
            //Now ignore (read) the two dummy characters
            ignoreChar(2);
            
        }
        
        

        
    }
    
    
    void addTriangle(const Vec3D normal, const Vec3D vertex1, const Vec3D vertex2, const Vec3D vertex3)
    {
        STLTriangle triangleToAdd(normal,vertex1,vertex2,vertex3);
        Triangles_.push_back(triangleToAdd);
    }
    
    void output()
    {
        
        std::cout << header_ << std::endl;
        for (unsigned int i=0;i<(numTriangles_);i++)
        {
            
            std::cout << "Triangle " << i << " has normal "<<Triangles_[i].normal <<" and vertexes" << Triangles_[i].vertex1 << Triangles_[i].vertex2 << Triangles_[i].vertex3 << std::endl;
            
        }
        
        
        
    }
    
    
    STLTriangle getTriangle(const unsigned int num)
    {
        return Triangles_[num];
    }
    
private:
    
    std::string header_;
    unsigned numTriangles_;
    std::vector<STLTriangle> Triangles_;
    
    
    
    

};


int main()
{
    std::string directory(getMercurySourceDir());
    STLReader myReader(directory+"/Drivers/ImportTools/ExampleSTLFiles/Box1x1x1.stl");

    std::vector<STLTriangle> answer;
    
    
    answer.emplace_back(Vec3D(0,1,0),Vec3D(1.60291, 1.36214, 1e-06),Vec3D(0.602907, 1.36214, 1e-06),Vec3D(1.60291, 1.36214, 1));
    
    answer.emplace_back(Vec3D(0, 1, -0),Vec3D(1.60291, 1.36214, 1),Vec3D(0.602907, 1.36214, 1e-06),Vec3D(0.602907, 1.36214, 1));
    
    answer.emplace_back(Vec3D(1, 0, -0),Vec3D(1.60291, 0.362136, 1e-06),Vec3D(1.60291, 1.36214, 1e-06), Vec3D(1.60291, 0.362136, 1));
    
    answer.emplace_back(Vec3D(1, -0, 0),Vec3D(1.60291, 0.362136, 1),Vec3D(1.60291, 1.36214, 1e-06),Vec3D(1.60291, 1.36214, 1));
    
    answer.emplace_back(Vec3D(0, -1, 0),Vec3D(0.602907, 0.362136, 1e-06),Vec3D(1.60291, 0.362136, 1e-06),Vec3D(0.602907, 0.362136, 1));
    
    answer.emplace_back(Vec3D(0, -1, 0),Vec3D(0.602907, 0.362136, 1),Vec3D(1.60291, 0.362136, 1e-06),Vec3D(1.60291, 0.362136, 1));
    
    answer.emplace_back(Vec3D(-1, 1.11022e-16, 0),Vec3D(0.602907, 1.36214, 1e-06),Vec3D(0.602907, 0.362136, 1e-06),Vec3D(0.602907, 1.36214, 1));
    
    answer.emplace_back(Vec3D(-1, 1.11022e-16, 0),Vec3D(0.602907, 1.36214, 1),Vec3D(0.602907, 0.362136, 1e-06),Vec3D(0.602907, 0.362136, 1));
    
    answer.emplace_back(Vec3D(0, -0, 1),Vec3D(0.602907, 0.362136, 1),Vec3D(1.60291, 0.362136, 1),Vec3D(0.602907, 1.36214, 1));
    
    answer.emplace_back(Vec3D(-0, 0, 1),Vec3D(0.602907, 1.36214, 1),Vec3D(1.60291, 0.362136, 1),Vec3D(1.60291, 1.36214, 1));
    
    answer.emplace_back(Vec3D(0, -0, -1),Vec3D(0.602907, 1.36214, 1e-06),Vec3D(1.60291, 1.36214, 1e-06),Vec3D(0.602907, 0.362136, 1e-06));
    
    answer.emplace_back(Vec3D(0, 0, -1),Vec3D(0.602907, 0.362136, 1e-06),Vec3D(1.602907, 1.362136, 1e-06),Vec3D(1.60291, 0.362136, 1e-06));
    
  

    
    for (unsigned int i=0;i<answer.size();i++)
    {
//        std::cout << answer[i].normal << std::endl;
//        std::cout <<myReader.getTriangle(i).normal << std::endl;
//        
//        std::cout << answer[i].vertex1 << std::endl;
//        std::cout <<myReader.getTriangle(i).vertex1 << std::endl;
//        
//        std::cout << answer[i].vertex2 << std::endl;
//        std::cout <<myReader.getTriangle(i).vertex2 << std::endl;
//        
//        std::cout << answer[i].vertex3 << std::endl;
//        std::cout <<myReader.getTriangle(i).vertex3 << std::endl;
        
        
        if (!(answer[i].isEqualTo(myReader.getTriangle(i),1e-4)))
            logger(FATAL,"Triangle % has need read incorrectly",i);
        
        
        
    }
   
   
    
    
}
