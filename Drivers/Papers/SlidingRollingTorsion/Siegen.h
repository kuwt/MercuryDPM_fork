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

#include "DPMBase.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <iomanip>
#include "Species/LinearViscoelasticFrictionSpecies.h"

class Siegen : public DPMBase{
public:

	//We define the particle properties for the siegen experiments and insert one particle
	Siegen(double NormalForce_=1000e-6) : DPMBase() {
		// User input: problem name
		setName("Siegen");
		NormalForce=NormalForce_;
		
		// User input: material parameters
		// radius of the particle
		Mdouble RadiusPrime = 10e-6; //m
		//calculate rel overlap
		Mdouble E1=179e9;
        Mdouble E2=71e9;
        Mdouble nu1=0.17;
        Mdouble nu2=0.17;
        Mdouble Estar = 1/( (1-nu1*nu1)/E1 + (1-nu2*nu2)/E2 );

        Mdouble G1 = E1/2/(1+nu1);
        Mdouble G2 = E2/2/(1+nu2);
        Mdouble Gstar = 1/( (2-nu1)/G1 + (2-nu2)/G2 );
        //Mdouble Estar = 52.3488827e9;
		//Mdouble Poisson = 0.17;
		//Mdouble Gstar = 11.871e9;
		std::cout << "Estar=" << Estar << std::endl;
		std::cout << "Gstar=" << Gstar << std::endl;
		Gstar = 1.2138e9;
		// relative overlap for given normal force
		Mdouble Overlap=0;
		for (int i=0; i<10; i++) {
			Overlap = pow(3./4.*NormalForce/Estar/sqrt(RadiusPrime+Overlap),2./3.);
		}
		relOverlap = Overlap/RadiusPrime; //wrt RadiusPrime
		relOverlap = 0.001; //wrt RadiusPrime
		Mdouble Radius = RadiusPrime*(1+relOverlap);
		std::cout << "relOverlap=" << relOverlap << std::endl;
		std::cout << "overlap=" << relOverlap*RadiusPrime << std::endl;
		std::cout << "contact radius=" << sqrt(relOverlap)*RadiusPrime << std::endl;
		std::cout << "rel contact radius=" << sqrt(relOverlap) << std::endl;
		//density of SiO2

        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(2648/mathsFunc::cubic(1+relOverlap)); //kg/m^3
        //(dynamic) sliding friction
		//setSlidingFrictionCoefficient(.23);
        species->setSlidingFrictionCoefficient(0.23); //rough
		//species->setSlidingFrictionCoefficientStatic(.7); //rough
		//rolling friction
        species->setRollingFrictionCoefficient(.00103);
		//torsion friction
        species->setTorsionFrictionCoefficient(.005);
		//restitution coefficient
		Mdouble Restitution=0.99;
                
		//insert particle and get mass

		setParticleDimensions(3);
        setSystemDimensions(3);
		BaseParticle* p0 = particleHandler.copyAndAddObject(SphericalParticle());
        p0->setRadius(Radius);
        Mdouble Mass = p0->getMass();

        //define material properties
		//set normal force
        //setStiffnessAndRestitutionCoefficient(NormalForce/(relOverlap*RadiusPrime), Restitution, Mass);
        species->setStiffnessAndRestitutionCoefficient(NormalForce/relOverlap/RadiusPrime, Restitution, Mass);
        //set dt accordingly
        Mdouble tc = species->getCollisionTime(Mass);
        //setTimeStep(tc/50.);
        setTimeStep(7.399163453192103488e-08);
        std::cout << "Mass=" << Mass << std::endl;
		std::cout << "tc=" << std::setprecision(19)<< tc << std::endl;
		std::cout << "dt=" << getTimeStep() << std::endl;
		std::cout << "tmax=" << getTimeMax() << std::endl;
		//set elastic tangential sliding force
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc,Restitution,Restitution, Mass);
		//set elastic tangential rolling/torsion force
		Mdouble factor=Gstar/Estar; //0.4 is for equal oscillation frequency
		factor=0.4; //0.4 is for equal oscillation frequency
        species->setRollingStiffness(factor*species->getStiffness());
        species->setTorsionStiffness(factor*species->getStiffness());

//		if (false) { //Hertz-Mindlin
//            species->setStiffness(Estar);
//            species->setSlidingStiffness(Gstar);
//            Mdouble disp=0.0062; //en=0.99
//            //Mdouble disp=0.0382; //en=0.88
//            species->setDissipation(disp);
//            species->setSlidingDissipation(disp);
//			std::cout << "E*=" << species->getStiffness() << ", G*/E*=" << species->getSlidingStiffness()/species->getStiffness() << std::endl;
//			//setRollingFrictionCoefficient(0);
//			//setTorsionFrictionCoefficient(0);
//			//setSlidingFrictionCoefficient(0);
//            ///\todo TWnow reimplement HERTZ_MINDLIN_DERESIEWICZ
//            //species->setForceType(ForceType::HERTZ_MINDLIN_DERESIEWICZ);
//		}

        species->setSlidingDissipation(sqrt(factor)*species->getDissipation());
        species->setTorsionDissipation(sqrt(factor)*species->getDissipation());
        species->setRollingDissipation(sqrt(factor)*species->getDissipation());
		
		//set geometry
		setGravity(Vec3D(0,0,0));
		setXMax(Radius);
		setXMin(-getXMax());
		setYMax(2.*Radius);
		setYMin(0);
		setZMax(Radius);
		setZMin(-getZMax());
	}
	
	void setupInitialConditions() {}

	Mdouble relOverlap;
	Mdouble NormalForce;
	Mdouble RelLoopSize;
	Mdouble Angle;
	Mdouble TangentialVelocity;
	Mdouble LoopTime;
    LinearViscoelasticFrictionSpecies* species;
};
