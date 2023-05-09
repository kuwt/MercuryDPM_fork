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

#include <Mercury3D.h>
#include <Species/HertzianViscoelasticMindlinRollingTorsionSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/IntersectionOfWalls.h>

class TwoParticleContact : public Mercury3D
{
protected:
    Mdouble radius_ = 1.5e-3;
    HertzianViscoelasticMindlinRollingTorsionSpecies* species_;

public:

    void setupInitialConditions() override
    {
        //set general properties
        setMin(Vec3D(-2, -1, -1) * radius_);
        setMax(Vec3D(2, 1, 1) * radius_);
    
        //define material/contact properties
        HertzianViscoelasticMindlinRollingTorsionSpecies s;
        s.setDensity(950);
    
        Mdouble poisson = 0.4;
        Mdouble shearModulus = 1e8;
        Mdouble restitutionCoeff = 1e-4;
        s.setEffectiveElasticModulusAndRestitutionCoefficient(
                shearModulus * (3 * poisson + 2 * shearModulus) / (poisson + shearModulus), restitutionCoeff);
    
        s.setSlidingFrictionCoefficient(.5);
        s.setRollingFrictionCoefficient(.1);
        s.setTorsionFrictionCoefficient(0);
        s.setEffectiveShearModulus(shearModulus);
        //s.setRollingStiffness(s.getElasticModulus()*radius_/300.);
        //s.setTorsionStiffness(s.getElasticModulus()*radius_/300.);
        //s.setSlidingDissipation(s.getDissipation());
        //s.setRollingDissipation(s.getDissipation());
        //s.setTorsionDissipation(s.getDissipation());
        species_ = speciesHandler.copyAndAddObject(s);

        //two particles
        Mdouble relativeVelocity = 0.1;
        SphericalParticle p(species_);
        p.setRadius(radius_);
        p.setPosition(Vec3D(-1,0,0)*radius_);
        p.setVelocity(Vec3D(1,0,0)*relativeVelocity);
        particleHandler.copyAndAddObject(p);

        p.setPosition(-p.getPosition());
        p.setVelocity(-p.getVelocity());
        particleHandler.copyAndAddObject(p);

        //time-stepping properties
        //setTimeStep(0.02*species_->getCollisionTime(2.0*radius_,species_->getDensity(),relativeVelocity));
        setTimeStep(1.4e-06); //10% of Rayleigh time
        setTimeMax(500.0*getTimeStep());
        setSaveCount(1);
    }

    void test() {
        setName("TwoParticleContact");
        solve();

        helpers::writeToFile("TwoParticleContact.gnu","set logscale y\n"
                "p 'TwoParticleContact.ene' u 1:(sqrt($3*1e7))");

        helpers::check(particleHandler.getLastObject()->getAngularVelocity().Y,0,1e-2,
                       "TwoParticleContact: Particle has wrong angular velocity");
    }
};

class RollingContact : public TwoParticleContact
{
public:

    void setupInitialConditions() override
    {
        TwoParticleContact::setupInitialConditions();

        particleHandler.clear();

        //one particle
        SphericalParticle p(species_);
        p.setRadius(radius_);
        p.setPosition(Vec3D(0,0,0)*radius_);
        particleHandler.copyAndAddObject(p);

        //walls in z-direction
        wallHandler.copyAndAddObject(InfiniteWall(Vec3D(0,0,-1),getMin(),species_));

        //gravity
        Mdouble angle = 0.8*constants::pi/3.0;
        setGravity(Vec3D(sin(angle), 0, -cos(angle))*10);

        //time-stepping properties
        setTimeMax(50000.0*getTimeStep());
        setSaveCount(50);
    }

    void printTime() const override
    {
        Mdouble v = particleHandler.getLastObject()->getVelocity().X;
        Mdouble w = particleHandler.getLastObject()->getAngularVelocity().Y;
        Mdouble r = particleHandler.getLastObject()->getRadius();
        logger(INFO,"t % v % w % vrel %",getTime(),v,w,v-r*w);
    }

    void test()
    {
        setName("RollingContact");
        solve();

        //plot position
        helpers::writeToFile("RollingContactPosition.gnu","set logscale y\n"
                "p 'RollingContact.ene' u 6:8");
        //plot energy conservation
        helpers::writeToFile("RollingContact.gnu",""
                "p 'RollingContact.ene' u 1:($3+$4) t 'Ekin', '' u 1:(-0.81*$2) w l t 'Epot'");

        helpers::check(particleHandler.getLastObject()->getVelocity().Y,0,1e-2,
                       "RollingContact: Particle has wrong velocity");
    }
};


int main()
{
    TwoParticleContact dpm;
    dpm.test();

    RollingContact dpm2;
    dpm2.test();

    return 0;
}
