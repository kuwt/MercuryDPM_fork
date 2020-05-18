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

#include "Species/LinearViscoelasticSpecies.h"
#include "DPMBase.h"
using constants::pi;

void test1() {
    std::cout << "Test 1:\n"
        "We apply a torque T around the z-axis to a motionless\n"
        "particle of inertia I for a time t using a time step dt.\n"
        "The result is a final angular velocity omega=T/I*t\n"
        "and rotation angle alpha=T/I*t^2/2.\n"
        "\n"
        "The input/output values should be:\n"
        "   T = 0.1*pi*(0 0 1) Nm\n"
        "   I = 0.1 kg m^2\n"
        "   t = 1 s\n"
        "   dt = 1e-5 s\n"
        "   omega=pi*(0 0 1) rad/s\n"
        "   alpha = pi/2*(0 0 1) rad\n"
        << std::endl;

    //compute mass (requires species and handler to be set)
    DPMBase D;
    LinearViscoelasticSpecies* S = D.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    BaseParticle* P = D.particleHandler.copyAndAddObject(SphericalParticle(S));
    D.setDimension(3); ///\todo shouldn't this be default?
    S->setDensity(6.0/constants::pi);
    P->setRadius(0.5);
    S->computeMass(P);
    
    //std::cout << "P " << *P << std::endl;
    P->setForce({1,0,0});
    P->setTorque(0.1*constants::pi*Vec3D(0,0,1));
    Mdouble timeMax = 1.0;
    Mdouble timeStep = 1e-5;
    for (Mdouble time = 0; time<timeMax; time+=timeStep) {
        P->accelerate(P->getForce() * P->getInvMass() * timeStep);
        P->move(P->getVelocity() * timeStep);   
        P->angularAccelerate(P->getOrientation().rotateInverseInertiaTensor(P->getInvInertia())*P->getTorque() * timeStep);
        P->rotate(P->getAngularVelocity() * timeStep);
    }
    
    std::cout << "Results:\n"
        << "   omega= " << P->getAngularVelocity() << std::endl
        << "   alpha= " << P->getOrientation().getEuler() << std::endl;
    
    //    std::cout << P->getAngularVelocity() - Vec3D(0,0,pi) << std::endl;
    //    std::cout << P->getOrientation().getEuler() - Vec3D(0,0,pi/2) << std::endl;
    
    //check results (with numerical error values added))
    if (!mathsFunc::isEqual(P->getAngularVelocity(), Vec3D(0,0,pi+3.14159e-05), 1e-10))
    {
        logger(ERROR, "angular velocity is %, but should be %", P->getAngularVelocity(), pi/2);
    }
    if (!mathsFunc::isEqual(P->getOrientation().getEuler(), {0,0,pi/2+4.71241e-05}, 1e-10))
    {
        logger(ERROR, "orientation is %, but should be %", P->getOrientation().getEuler().Z, pi/2);
    }
}

void test2() {
    std::cout << "Test 2:\n"
        "We apply a torque T around the z-axis to a motionless\n"
        "particle of inertia I for a time t using a time step dt.\n"
        "The result is a final angular velocity omega=inv(I)*T*t\n"
        "and rotation angle alpha=inv(I)*T*t^2/2.\n"
        "\n"
        "The input/output values should be:\n"
        "   T = 0.1*pi*(0 0 1) Nm\n"
        "   I = 0.1*(2 0 0;0 2 0;0 0 1) kg m^2\n"
        "   t = 1 s\n"
        "   dt = 1e-5 s\n"
        "   omega=pi*(0 0 1) rad/s\n"
        "   alpha = pi/2*(0 0 1) rad\n"
        << std::endl;

    //compute mass (requires species and handler to be set)
    DPMBase D;
    LinearViscoelasticSpecies* S = D.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    BaseParticle* P = D.particleHandler.copyAndAddObject(SphericalParticle(S));
    D.setDimension(3); ///\todo shouldn't this be default?
    S->setDensity(6.0/constants::pi);
    P->setRadius(0.5);
    S->computeMass(P);
    P->setInertia(MatrixSymmetric3D(2,1,0,
                                      2,0,
                                        1)*0.1);
    
    //std::cout << "P " << *P << std::endl;
    P->setForce({1,0,0});
    P->setTorque(0.1*constants::pi*Vec3D(0,0,1));
    Mdouble timeMax = 1.0;
    Mdouble timeStep = 1e-5;
    for (Mdouble time = 0; time<timeMax; time+=timeStep) {
        P->accelerate(P->getForce() * P->getInvMass() * timeStep);
        P->move(P->getVelocity() * timeStep);   
        P->angularAccelerate(P->getOrientation().rotateInverseInertiaTensor(P->getInvInertia())*P->getTorque() * timeStep);
        P->rotate(P->getAngularVelocity() * timeStep);
    }
    
    std::cout << "Results:\n"
        << "   omega= " << P->getAngularVelocity() << std::endl
        << "   alpha= " << P->getOrientation().getEuler() << std::endl;
    
    //    std::cout << P->getAngularVelocity() - Vec3D(0,0,pi) << std::endl;
    //    std::cout << P->getOrientation().getEuler() - Vec3D(0,0,pi/2) << std::endl;
    
    //check results (with numerical error values added))
    if (!mathsFunc::isEqual(P->getAngularVelocity(), Vec3D(0,0,pi+3.14159e-05), 1e-10))
    {
        logger(ERROR, "angular velocity is %, but should be %", P->getAngularVelocity(), pi/2);
    }
    if (!mathsFunc::isEqual(P->getOrientation().getEuler(), {0,0,pi/2+4.71241e-05}, 1e-10))
    {
        logger(ERROR, "orientation is %, but should be %", P->getOrientation().getEuler().Z, pi/2);
    }
}

int main(int argc, char** argv) {
    test1();
    test2();
    return 0;
}
