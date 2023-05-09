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
#include "Walls/InfiniteWall.h"
#include "Species/Species.h"
#include "Species/HertzianBSHPViscoelasticSpecies.h"

/*
    In this file particles with the species HertzianBSHPViscoelasticSpecies are 
    colliding with given velocities. The coefficient of restitution is determined
    from the resulting velocities and checked against the theoretical value calculated
    with Pade approximations. See https://doi.org/10.1103/PhysRevE.84.021302
    https://doi.org/10.1103/PhysRevE.84.021302
*/

class HertzianBSHPInteractionTwoParticleElasticCollision : public Mercury3D
{
    /*!
     * Calculates the coefficient of restitution of the HertzianBSHPViscoelasticNormalInteraction
     * based on pade approximations. For a description see https://doi.org/10.1103/PhysRevE.84.021302
     */
    Mdouble padeCoefficientOfRestitution(Mdouble v, Mdouble A, Mdouble meff, Mdouble Reff, Mdouble Eeff)
    {
        Mdouble rho = 4. / 3. * Eeff * std::sqrt(Reff);
        Mdouble beta = 3. / 2. * A * std::pow(rho/meff, 2. /5.);
        Mdouble vStar = std::pow(beta, 1./2.) * std::pow(v, 1. / 10.);
        
        std::vector<Mdouble> aCoeff = {1.0, 1.07232, 0.574198, 0.141552};
        std::vector<Mdouble> bCoeff = {1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585};
        
        Mdouble enumerator = 0;
        Mdouble denominator = 0;
        
        for (unsigned int i=0; i<4; i++)
        {
            enumerator += aCoeff[i] * std::pow(vStar, i);
        }
        
        for (unsigned int i=0; i<7; i++)
        {
            denominator += bCoeff[i] * std::pow(vStar, i);
        }
        
        Mdouble epsilon = enumerator / denominator;
        return epsilon;

    }

    void setupInitialConditions() override {
        setMax({0.1,0.1,0.1});
        setMin({-0.1,-0.1,-0.1});
        
        setGravity({0.0,0.0,0.0});
        Mdouble relativeVelocity = initialRelativeVelocity_;
        
        Mdouble radius = 1e-3;
        Mdouble ySpacing = 8*radius;
        Mdouble xSpacing = 4*radius;
        
        SphericalParticle p0, p1;
        p0.setSpecies(speciesHandler.getObject(0));
        p1.setSpecies(speciesHandler.getObject(0));

        p0.setRadius(1e-3);
        p1.setRadius(1e-3);
        
        Mdouble Reff = p0.getRadius()*p1.getRadius() / (p0.getRadius()+p1.getRadius());
        Mdouble Eeff = dynamic_cast<HertzianBSHPViscoelasticSpecies*>(speciesHandler.getObject(0))->getEffectiveElasticModulus();
        Mdouble meff = p0.getMass()*p1.getMass() / (p0.getMass()+p1.getMass());
        Mdouble A = dynamic_cast<HertzianBSHPViscoelasticSpecies*>(speciesHandler.getObject(0))->getDissipation();
        
        for (unsigned int i = 0; i < initialRelativeVelocities_.size(); i++)
        {
            p0.setPosition(Vec3D(xSpacing/2, ySpacing*i, 0.0));
            p1.setPosition(Vec3D(-xSpacing/2, ySpacing*i, 0.0));
            
            p0.setVelocity(Vec3D(-initialRelativeVelocities_[i]/2, 0.0, 0.0));
            p1.setVelocity(Vec3D(initialRelativeVelocities_[i]/2, 0.0, 0.0));

            particleHandler.copyAndAddObject(p0);
            particleHandler.copyAndAddObject(p1);
            
            PadeCor_.push_back(padeCoefficientOfRestitution(initialRelativeVelocities_[i], A, meff, Reff, Eeff));
        }
        
        wallHandler.clear();
    }
    
    void actionsAfterSolve() override
    {
        {
            for (unsigned int i = 0; i < initialRelativeVelocities_.size(); i++)
            {
                Mdouble relativeVelocity = particleHandler.getObjectById(2*i+0)->getVelocity().X - particleHandler.getObjectById(2*i+1)->getVelocity().X;
                logger(INFO, "Initial relative velocity is %", initialRelativeVelocities_[i]);
                logger(INFO, "    Numerical cor is %", relativeVelocity / initialRelativeVelocities_[i]);
                logger(INFO, "    Cor should be %", PadeCor_[i]);
                
                Mdouble numericalCOR = relativeVelocity / initialRelativeVelocities_[i];
                
                Mdouble relativeDeviation = std::fabs(numericalCOR-PadeCor_[i])/PadeCor_[i];
                logger.assert_always(relativeDeviation < 0.01, "Error: Determined coefficient of restitution for vel % deviates from Pade approximation.", initialRelativeVelocities_[i]);
            }
            
        }
    }
    
private:
    Mdouble cor_;
    Mdouble initialRelativeVelocity_ = 1.0;
    std::vector<Mdouble> initialRelativeVelocities_ = {0.1, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0};
    std::vector<Mdouble> PadeCor_;
};

int main(int argc UNUSED, char* argv[] UNUSED)
{

    HertzianBSHPInteractionTwoParticleElasticCollision twoParticleElasticCollisionProblem;
    twoParticleElasticCollisionProblem.setName("HertzianBSHPInteractionTwoParticleElasticCollision");

    HertzianBSHPViscoelasticSpecies species;
    species.setDensity(2000);
    species.setEffectiveElasticModulus(1e8);
    species.setDissipation(1e-5);
    
    twoParticleElasticCollisionProblem.speciesHandler.copyAndAddObject(species);

    twoParticleElasticCollisionProblem.setTimeMax(0.1);
    twoParticleElasticCollisionProblem.setSaveCount(1e5);
    twoParticleElasticCollisionProblem.setTimeStep(1e-6);
    twoParticleElasticCollisionProblem.fStatFile.setFileType(FileType::NO_FILE);
    twoParticleElasticCollisionProblem.solve();
}
