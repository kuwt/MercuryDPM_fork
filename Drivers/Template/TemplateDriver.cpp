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

#include <Mercury3D.h>

/*
* This is our problem description. Everything is set up here.
* We inherit from Mercury3D, since this gives full flexibility.
* For more predefined problems (for instance, chutes), please see the
* documentation.
*/
class MyProblem : public Mercury3D
{
public:
 /* We define our own 'setupInitialConditions' function here,
  * which defines all the specifics about our simulation here.
  */
 void setupInitialConditions() override
 {
   //The first step: set any properties which are always true for your system.
   // (for instance, if gravity remains constant, set it here)
   
   
   //Now, decide what Species you need for your system.
   
   
   //Add your walls below, and don't forget to set the species!
   
   
   //Now, either add particles or Insertion/Deletion boundaries
   
   
 }

};

int main(const int argc, const char** argv)
{
  MyProblem problem;
  /* Make sure to add an (unique) name here */
  problem.setName("MyProblemCode");

  //Add any configuration which is modified often below...
  
  /* Now start the simulation */
  problem.solve();
  return 0;
}
