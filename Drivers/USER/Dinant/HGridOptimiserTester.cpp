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

#include "Mercury3D.h"
#include "HGridOptimiser.h"
#include <iostream>
#include <vector>
#include "Species/LinearViscoelasticSpecies.h"

class HGridOptimiserTester : public Mercury3D
{
public:

    HGridOptimiserTester()
    {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

	void setupInitialConditions() override
	{
		if(particleHandler.getNumberOfObjects()!=N)
        {
            particleHandler.clear();
            
            double VP=constants::pi*4.0/3.0;
            L=std::pow(N*VP*DistInt(1,omega)/nu,1.0/3.0);

            setXMin(0);
            setYMin(0);
            setZMin(0);
            setXMax(L);
            setYMax(L);
            setZMax(L);

            particleHandler.setStorageCapacity(2*N);
            SphericalParticle p0;
            p0.setVelocity(Vec3D(0.0,0.0,0.0));

            double V=0;
            
            //Use at least particles with maximum and minimum size
            p0.setRadius(1.0);
            Vec3D position;
            position.Z = random.getRandomNumber(0,getZMax());
            position.Y = random.getRandomNumber(0,getYMax());
            position.X = random.getRandomNumber(0,getXMax());
            p0.setPosition(position);
            particleHandler.copyAndAddObject(p0);
            V+=VP*std::pow(p0.getRadius(),3);
            
            p0.setRadius(omega);
            position.Z = random.getRandomNumber(0,getZMax());
            position.Y = random.getRandomNumber(0,getYMax());
            position.X = random.getRandomNumber(0,getXMax());
            p0.setPosition(position);
            particleHandler.copyAndAddObject(p0);
            V+=VP*std::pow(p0.getRadius(),3);

            //For other particles use a random distribution
            for(unsigned int i=2;i<N;i++)
            {
                p0.setRadius(RandomRadius());
                position.Z = random.getRandomNumber(0,getZMax());
                position.Y = random.getRandomNumber(0,getYMax());
                position.X = random.getRandomNumber(0,getXMax());
                p0.setPosition(position);
                particleHandler.copyAndAddObject(p0);
                V+=VP*std::pow(p0.getRadius(),3);
            }
        }
	}
    
    double RandomRadius()
    {
        double rand=random.getRandomNumber(0,1);
        if(alpha==-1)
        {
            return std::pow(omega,rand);
        }
        else
        {
            return std::pow(rand*(std::pow(omega,1.0+alpha)-1.0)+1.0,1.0/(1.0+alpha));
        }
    }
    
    double DistInt(double s, double e)
    {
        if(omega==1)
        {
            return 1;
        }
        double teller;
        double noemer;
        if(alpha==-1)
        {
            noemer=std::log(omega);
        }
        else
        {
            noemer=(std::pow(omega,1.0+alpha)-1.0)/(1.0+alpha);
        }
        
        if(alpha==-4)
        {
            teller=std::log(e)-std::log(s);
        }
        else
        {
            teller=(std::pow(e,4.0+alpha)-std::pow(s,4.0+alpha))/(4.0+alpha);
        }
        return teller/noemer;
    }

    double L;
    double omega;
    double alpha;
    double nu;
    unsigned int N;
public:
    LinearViscoelasticSpecies* species;
};

void testDerivativeFunctions()
{
		{
		HGridOptimiser optimiser1,optimiser2,optimiser3;
		double omega=10;
		int order=3;
		int p=3;
		
		std::vector<double> coeff;
		std::cout<<"Coefficients are:";
		for(int i=0;i<=order;i++)
		{
			coeff.push_back((1.0+i)/(1.0+order)/(std::pow(omega,i+1)-1));
			std::cout<<" "<<coeff.back();
		}
		std::cout<<std::endl;
		optimiser1.initialisePolyFunc(omega,coeff,100,0);
		optimiser2.initialisePolyFunc(omega,coeff,1000,0);
		optimiser3.initialisePolyFunc(omega,coeff,10000,0);
		double start=omega/std::sqrt(21);
		double end=omega/std::sqrt(11);
		double h=omega/std::sqrt(13);
		
		{
			double val1=optimiser1.pdfInt(start,end,p);
			double val2=optimiser2.pdfInt(start,end,p);
			double val3=optimiser3.pdfInt(start,end,p);
			double analval=0;
			for(int j=0;j<=order;j++)
			{
				analval+=coeff[j]/(p+j+1)*(std::pow(end,p+j+1)-std::pow(start,p+j+1));
			}			
			std::cout<<"T pdfInt("<<start<<","<<end<<","<<p<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}		

		{
			double val1=optimiser1.diffPdfInt(end,p);
			double val2=optimiser2.diffPdfInt(end,p);
			double val3=optimiser3.diffPdfInt(end,p);
			double analval=0;
			for(int j=0;j<=order;j++)
			{
				analval+=coeff[j]*std::pow(end,p+j);
			}			
			std::cout<<"T diffPdfInt("<<end<<","<<p<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}
		
		{
			double val1=optimiser1.expectedCellsIntegral(start,end,p,h);
			double val2=optimiser2.expectedCellsIntegral(start,end,p,h);
			double val3=optimiser3.expectedCellsIntegral(start,end,p,h);
			
			double teller=0;
			double noemer=0;
			for(int j=0;j<=order;j++)
			{
				noemer+=coeff[j]/(1.0+j)*(std::pow(end,j+1)-std::pow(start,j+1));
				double fact=coeff[j]*0.5*h/(p+1);
				for(int i=0;i<=j;i++)
				{
					teller+=fact*(std::pow(end,j-i)*std::pow(2*end/h+2,p+1+i)-std::pow(start,j-i)*std::pow(2*start/h+2,p+1+i));
					fact*=-0.5*h*(j-i)/(p+2+i);
				}
			}
			double analval=teller/noemer;
			std::cout<<"T expectedCellsIntegral("<<start<<","<<end<<","<<p<<","<<h<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}			

		{
			double val1=optimiser1.diffStartExpectedCellsIntegral(start,end,p,h);
			double val2=optimiser2.diffStartExpectedCellsIntegral(start,end,p,h);
			double val3=optimiser3.diffStartExpectedCellsIntegral(start,end,p,h);
			
			double teller=0;
			double diffteller=0;
			double noemer=0;
			double diffnoemer=0;
			for(int j=0;j<=order;j++)
			{
				noemer+=coeff[j]/(1.0+j)*(std::pow(end,j+1)-std::pow(start,j+1));
				diffnoemer+=-coeff[j]*std::pow(start,j);
				double fact=coeff[j]*0.5*h/(p+1);
				for(int i=0;i<=j;i++)
				{
					teller+=fact*(std::pow(end,j-i)*std::pow(2*end/h+2,p+1+i)-std::pow(start,j-i)*std::pow(2*start/h+2,p+1+i));
					diffteller+=-fact*(j-i)*std::pow(start,j-i-1)*            std::pow(2*start/h+2,p+1+i);
					diffteller+=-fact*      std::pow(start,j-i  )*(p+1+i)*2/h*std::pow(2*start/h+2,p+i  );
					fact*=-0.5*h*(j-i)/(p+2+i);
				}
			}
			double analval=(diffteller*noemer-teller*diffnoemer)/(noemer*noemer);
			std::cout<<"T diffStartExpectedCellsIntegral("<<start<<","<<end<<","<<p<<","<<h<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}			
	
		{
			double val1=optimiser1.diffEndExpectedCellsIntegral(start,end,p,h);
			double val2=optimiser2.diffEndExpectedCellsIntegral(start,end,p,h);
			double val3=optimiser3.diffEndExpectedCellsIntegral(start,end,p,h);
			
			double teller=0;
			double diffteller=0;
			double noemer=0;
			double diffnoemer=0;
			for(int j=0;j<=order;j++)
			{
				noemer+=coeff[j]/(1.0+j)*(std::pow(end,j+1)-std::pow(start,j+1));
				diffnoemer+=coeff[j]*std::pow(end,j);
				double fact=coeff[j]*0.5*h/(p+1);
				for(int i=0;i<=j;i++)
				{
					teller+=fact*(std::pow(end,j-i)*std::pow(2*end/h+2,p+1+i)-std::pow(start,j-i)*std::pow(2*start/h+2,p+1+i));
					diffteller+=fact*(j-i)*std::pow(end,j-i-1)*            std::pow(2*end/h+2,p+1+i);
					diffteller+=fact*      std::pow(end,j-i  )*(p+1+i)*2/h*std::pow(2*end/h+2,p+i  );
					fact*=-0.5*h*(j-i)/(p+2+i);
				}
			}
			double analval=(diffteller*noemer-teller*diffnoemer)/(noemer*noemer);
			std::cout<<"T diffEndExpectedCellsIntegral("<<start<<","<<end<<","<<p<<","<<h<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}
		
		{
			double val1=optimiser1.diffHExpectedCellsIntegral(start,end,p,h);
			double val2=optimiser2.diffHExpectedCellsIntegral(start,end,p,h);
			double val3=optimiser3.diffHExpectedCellsIntegral(start,end,p,h);
			
			double teller=0;
			double diffteller=0;
			double noemer=0;
			for(int j=0;j<=order;j++)
			{
				noemer+=coeff[j]/(1.0+j)*(std::pow(end,j+1)-std::pow(start,j+1));
				double fact=coeff[j]*0.5/(p+1);
				for(int i=0;i<=j;i++)
				{
					teller+=fact*std::pow(h,i+1)*(std::pow(end,j-i)*std::pow(2*end/h+2,p+1+i)-std::pow(start,j-i)*std::pow(2*start/h+2,p+1+i));
					diffteller+=fact*(-std::pow(h,i-1)*2*(p+1+i)*(std::pow(end,j-i+1)*std::pow(2*end/h+2,p+i)-std::pow(start,j-i+1)*std::pow(2*start/h+2,p+i))+
					                   std::pow(h,i)*(i+1)*(std::pow(end,j-i)*std::pow(2*end/h+2,p+1+i)-std::pow(start,j-i)*std::pow(2*start/h+2,p+1+i)));
					fact*=-0.5*(j-i)/(p+2+i);
				}
			}
			double analval=diffteller/noemer;
			std::cout<<"T diffHExpectedCellsIntegral("<<start<<","<<end<<","<<p<<","<<h<<") this should be "<<analval<<"="<<val1<<" "<<val2<<" "<<val3<<" diff="<<val1-analval<<" "<<val2-analval<<" "<<val3-analval<<std::endl;
		}			
	}
}

void basicOptimiser()
{
	HGridOptimiserTester problem;
    HGridOptimiser optimiser;
    
    problem.omega=10;
    problem.alpha=-2;
    problem.nu=0.7;
    problem.N=1000000;
    
    problem.setupInitialConditions();
    
    optimiser.initialise(problem,100,1);
    
    {
		std::vector<double> hGridCellSizes;
		std::vector<double> dfdx;
		optimiser.getOptimalDistribution(hGridCellSizes, 10, BOTTOMUP, 2);
		optimiser.calculateDiffWork(hGridCellSizes, dfdx, BOTTOMUP,1);   
		optimiser.calcDfDx(hGridCellSizes, dfdx, BOTTOMUP, 1);
		
		/*{
			optimiser.calculateDiffWork(hGridCellSizes, dfdx, BOTTOMUP, 1);
			double max=optimiser.checkLimit(hGridCellSizes,dfdx,1);
			int steps=1000;
			std::cout.precision(15);
			std::cout<<"Y=["<<optimiser.calculateWork(hGridCellSizes, BOTTOMUP,0);
			for(int i=0;i<steps;i++)
			{
				optimiser.applyStep(hGridCellSizes,dfdx,1e-4*max/steps,0);
				std::cout<<","<<optimiser.calculateWork(hGridCellSizes, BOTTOMUP,0);
			}
			std::cout<<"]"<<std::endl;
			optimiser.calculateWork(hGridCellSizes, BOTTOMUP,1);
			optimiser.applyStep(hGridCellSizes,dfdx,max/steps,0);
			optimiser.calculateWork(hGridCellSizes, BOTTOMUP,1);
		}*/
		//optimiser.getOptimalDistribution(hGridCellSizes, 10, TOPDOWN, 1);
		//optimiser.calculateDiffWork(hGridCellSizes, dfdx, TOPDOWN,1);
	}
}

int main(int argc UNUSED, char *argv[] UNUSED)
{
	basicOptimiser();

	
	/*HGridOptimiserTester problem;
    HGridOptimiser optimiser;
    
    
    problem.omega=10;
    problem.alpha=-1;
    problem.nu=0.7;
    problem.N=180;
    problem.NLevels=7;
    
    problem.setupInitialConditions();
    
    optimiser.Initialise(problem,0);*/
    
    
    /*{
		std::vector<double> hGridCellSizes;
		std::vector<double> dfdx;
		optimiser.getOptimalDistribution(hGridCellSizes, 10, BOTTOMUP, 1);
		//optimiser.calculateDiffWork(hGridCellSizes, dfdx, BOTTOMUP,1);   
		optimiser.getOptimalDistribution(hGridCellSizes, 10, TOPDOWN, 1);
		optimiser.calculateDiffWork(hGridCellSizes, dfdx, TOPDOWN,1);
		optimiser.calcDfDx(hGridCellSizes, dfdx, TOPDOWN, 1);
	}*/
	
	/*{   //Constant
		std::cout<<"This should be 1="<<optimiser.pdfInt(1,problem.omega,0)<<std::endl;
		double start=1.0;
		double end=0.5*(problem.omega+1);
		for(unsigned int p=0;p<5;p++)
		{
			double val=optimiser.pdfInt(start,end,p);
			double analval=(pow(end,p+1)-pow(start,p+1))/(p+1)/(problem.omega-1);
			std::cout<<"This should be "<<analval<<"="<<val<<" diff="<<val-analval<<std::endl;
		}
	}*/
	/*{   //Linear
		std::cout<<"This should be 1="<<optimiser.pdfInt(1,problem.omega,0)<<std::endl;
		double start=1.0;
		double end=0.5*(problem.omega+1);
		double a=1.0/(std::pow(problem.omega,2)-1);
		double b=0.5/(problem.omega-1);
		for(unsigned int p=0;p<5;p++)
		{
			double val=optimiser.pdfInt(start,end,p);
			double analval=(pow(end,p+2)-pow(start,p+2))*a/(p+2)+(pow(end,p+1)-pow(start,p+1))*b/(p+1);
			std::cout<<"This should be "<<analval<<"="<<val<<" diff="<<val-analval<<std::endl;
		}
	}*/
	
	//optimiser.checkPdfInt(3.6999,1);
    //optimiser.checkPdfInt(3.7001,1);
    /*{
        std::vector<double> dfdx;
        std::vector<double> hGridCellSizes;
        hGridCellSizes.clear();
        double rMin=problem.particleHandler.getSmallestParticle()->getRadius();;
        double rMax=problem.particleHandler.getLargestParticle()->getRadius();
        for (unsigned int i=0; i<problem.NLevels; i++)
        {
            hGridCellSizes.push_back(2.0*(rMin+(rMax-rMin)*i/problem.NLevels));
        }
        optimiser.calculateWork(hGridCellSizes, TOPDOWN,1);    
        //optimiser.calculateDiffWork(hGridCellSizes, dfdx, TOPDOWN,1);    
        //optimiser.calcDfDx(hGridCellSizes, dfdx, TOPDOWN,1);
    }*/
    
    /*{
        HGridMethod method=BOTTOMUP;
        std::cout<<"Test aditional step"<<std::endl;
        std::vector<double> dfdx;
        double nW,W;
        double epsilon=1e-1;
        double totalstep=0;
        //int loc=1;
        optimiser.calculateDiffWork(hGridCellSizes, dfdx, method,0);
        double max=optimiser.checkLimit(hGridCellSizes,dfdx,0);
        W=optimiser.calculateWork(hGridCellSizes, method,0);
        std::cout.precision(15);
        std::cout<<W<<std::endl;

        while(epsilon>1e-16)
        {
            optimiser.applyStep(hGridCellSizes,dfdx,-epsilon,0);
            totalstep-=epsilon;
            nW=optimiser.calculateWork(hGridCellSizes, method,0);
            while(nW<W)
            {
                optimiser.applyStep(hGridCellSizes,dfdx,-epsilon,0);
                totalstep-=epsilon;
                W=nW;
                nW=optimiser.calculateWork(hGridCellSizes, method,0);
            }
            optimiser.applyStep(hGridCellSizes,dfdx,epsilon,0);
            totalstep+=epsilon;
            epsilon/=10;
        }
        std::cout.precision(15);
        std::cout<<"step="<<totalstep<<" work="<<W<<std::endl;
        std::cout.precision(5);
        optimiser.calculateDiffWork(hGridCellSizes, dfdx, method,0); 
    }*/
}
