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

#include<iostream>

#include "DPMBase.h"
#include "Mercury2D.h"
#include "Math/ExtendedMath.h"


//#define DEBUG_OUTPUT

class my_problem : public DPMBase{
public:
	
	void setupInitialConditions() override
	{
	
		int N=particleHandler.getSize();
		int N1=static_cast<int>(sqrt(N))+1;
		
		for (int i=0;i<N;i++)
		{
		
			int ix=static_cast<int>(i%N1);
			int iy=static_cast<int>(i/N1);
		
			double x=(getXMax()-getXMin())*(ix+1)/(N1+1);
			double y=(getYMax()-getYMin())*(iy+1)/(N1+1);
			
		
			particleHandler.getObject(i)->setPosition(Vec3D(x,y,0.0));
			particleHandler.getObject(i)->setVelocity(Vec3D(0.1,0.1,0.0));
			particleHandler.getObject(i)->setRadius(0.0001);
		}
	
	}

private:

    void computeInternalForces(BaseParticle* i) override {}
};

class my_problem_HGRID : public Mercury2D{
public:

	my_problem_HGRID(my_problem& other) : DPMBase(other), Mercury2D(other) {}
	
	void setupInitialConditions() override
	{
	
		int N=particleHandler.getSize();
		int N1=static_cast<int>(sqrt(N))+1;
		
		for (int i=0;i<N;i++)
		{
		
			int ix=static_cast<int>(i%N1);
			int iy=static_cast<int>(i/N1);

			
		
			double x=(getXMax()-getXMin())*(ix+1)/(N1+1);
			double y=(getYMax()-getYMin())*(iy+1)/(N1+1);

            particleHandler.getObject(i)->setPosition(Vec3D(x,y,0.0));
            particleHandler.getObject(i)->setVelocity(Vec3D(0.1,0.1,0.0));
            particleHandler.getObject(i)->setRadius(0.0001);
		}	

	}

private:
    void computeInternalForces(BaseParticle* i) override {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::cout << "First test the gamma function is working" <<std::endl;
	std::cout << "---------------------------------------- \n" <<std::endl;
	std::cout << "Gamma(4)   =  " << mathsFunc::gamma(4) << "  : Error " << mathsFunc::gamma(4)-6<<std::endl;
 	std::cout << "Gamma(3.5) =  " << mathsFunc::gamma(3.5) << "  : Error" << mathsFunc::gamma(3.5)-3.323<<std::endl;
 	
 	std::cout << "Second test the chi-squared distribution is working" << std::endl;
 	std::cout << "---------------------------------------------------" << std::endl;
 	
 	std::cout << "\nFirst test chi(1.07,1) = " <<mathsFunc::chi_squared_prob(1.07,1) << "   Error " <<mathsFunc::chi_squared_prob(1.07,1)-0.3 << std::endl;
 	std::cout << "First test chi(1.07,1) = " <<mathsFunc::chi_squared_prob(6.25,3) << "   Error " <<mathsFunc::chi_squared_prob(6.25,3)-0.1 << std::endl;
 	std::cout << "First test chi(1.07,1) = " <<mathsFunc::chi_squared_prob(3.07,6) << "   Error " <<mathsFunc::chi_squared_prob(3.07,6)-0.8 << std::endl;
 	helpers::check(mathsFunc::chi_squared_prob(1.07,1)-0.3,0.0250886,1e-6,"Checking quality of chi-squared test");

	///Start off my solving the default problem
 	my_problem problem;
 	
 	problem.random.setRandomNumberGenerator(RNGType::LINEAR_CONGRUENTIAL_GENERATOR);
 
 	std::cout << "\nSecond test the actually random number generate : Linear Congruential Generator" << std::endl;
 	std::cout << "----------------------------------------------" << std::endl;
    std::cout << "First with the default parameters, prob numbers are from uniform = ";
    helpers::check(problem.random.test(),0.467846,1e-6,"Checking test result");

 	problem.random.setLinearCongruentialGeneratorParmeters(65539,0,1024*1024*1024*2-1);
 	std::cout << "Third test, now with Stefans' the default parameters, prob numbers are from uniform = ";
    helpers::check(problem.random.test(),0.131117,1e-6,"Checking test result");

 	problem.random.setRandomNumberGenerator(RNGType::LAGGED_FIBONACCI_GENERATOR);
 	problem.random.setLinearCongruentialGeneratorParmeters(1103515245,12345,1024*1024*1024);
 	problem.random.setRandomSeed(0);

 	std::cout << "\nForth test the actually random number generate : Lagged Fibonacci Generator" << std::endl;
 	std::cout << "----------------------------------------------" << std::endl;
    std::cout << "First with the default parameters, prob numbers are from uniform = ";
    helpers::check(problem.random.test(),0.617585,1e-6,"Checking test result");
}
