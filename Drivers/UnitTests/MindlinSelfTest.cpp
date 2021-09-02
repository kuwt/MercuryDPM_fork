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

#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Species/HertzianViscoelasticSlidingFrictionSpecies.h"
#include "DPMBase.h"

int main(int argc, char *argv[])
{
	Mdouble displacement=0.015;
	Mdouble tangentialDisplacement=0.0025;
	Mdouble velocity=0.01;
	Mdouble radius=0.25;

	//runs a loading - unloading - reloading test for normal and tangential forces, plus tests the objectiveness
	//- for Hertz-Mindlin
	HertzianViscoelasticMindlinSpecies species;
	species.setDensity(6./constants::pi);
    species.setEffectiveElasticModulusAndPoissonRatio(1e5, 0.3);
	species.setDissipation(2.0);
//	species.setSlidingDissipation(4.0/7.0);
	species.setSlidingFrictionCoefficient(0.1);

	helpers::loadingTest(&species, displacement, velocity, radius, "MindlinSelfTestLoading");
	helpers::normalAndTangentialLoadingTest(&species, displacement, tangentialDisplacement, velocity, radius,
											"MindlinSelfTestNormalAndTangentialLoading");
    helpers::objectivenessTest(&species, displacement, tangentialDisplacement, velocity, radius, "MindlinSelfTestFrameIndependence");

	/*
	 * //- for Hertz-SlidingFriction
	HertzianViscoelasticSlidingFrictionSpecies species2;
	species2.setDensity(6./constants::pi);
	species2.setEffectiveElasticModulus(1e5);
	species2.setDissipation(2.0);
	species2.setSlidingStiffness(300.0);
	species2.setSlidingFrictionCoefficient(1.0);
*/
	//Note: for correct signs, plot:
	//"MindlinSelfTestTangentialLoading.fstat" using 8:($10*$14)


	//helpers::loadingTest(&species2, displacement, velocity, radius);
    //helpers::normalAndTangentialLoadingTest(&species2, displacement, tangentialDisplacement, velocity, radius);
    //helpers::objectivenessTest(&species2, displacement, tangentialDisplacement, velocity, radius);
}
