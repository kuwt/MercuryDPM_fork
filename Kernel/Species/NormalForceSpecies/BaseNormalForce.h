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

#ifndef MERCURY_BASENORMALFORCE_H
#define MERCURY_BASENORMALFORCE_H
#include "Species/BaseForce.h"

class BaseNormalForce : public BaseForce
{
public:

    BaseNormalForce() {
        constantRestitution_ = false;
    }

    BaseNormalForce(const BaseNormalForce& p) {
        constantRestitution_ = p.constantRestitution_;
    }

    /*!
     * Accesses constantRestitution_.
     */
    bool getConstantRestitution() const {
        return constantRestitution_;
    }

    /*!
     * Sets constantRestitution_.
     */
    void setConstantRestitution(bool constantRestitution) {
        constantRestitution_ = constantRestitution;
        logger(INFO,"You are using a constant restitution contact model; make sure you understand the implications this has on the stiffness and dissipation used by the contact model!");
        //logger(INFO,"You are using a constant restitution contact model; if use setCollisionTime or setRestitutionCoefficient to set the parameters of your contact model, then these functions should be called AFTER calling setConstantRestitution!");
        //For more info, see http://docs.mercurydpm.org/Trunk/d0/db6/ConstantRestitution.html
    }

private:

    /*!
     * If constantRestitution_ is true, the elastic and dissipative force is multiplied by the harmonic mean mass, making restitution and collision time independent of the particle mass. This is set to false by default.
     */
    bool constantRestitution_;
};

#endif //MERCURY_BASENORMALFORCE_H
