//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "ParticleInfoTab.h"

ParticleInfoTab::ParticleInfoTab(int particleIndex, DPMBase* problem)
{
    this->particleIndex=particleIndex;
    this->problem=problem;

    positionLabel.set_text("Position");
    velocityLabel.set_text("Velocity");
    x.set_text("x");
    y.set_text("y");
    z.set_text("z");
    positionX.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().X));
    positionY.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().Y));
    positionZ.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().Z));
    velocityX.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().X));
    velocityY.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().Y));
    velocityZ.set_text(ToString(problem->particleHandler.getObject(particleIndex)->getPosition().Z));

    add(grid);
    grid.set_row_spacing(10);
    grid.set_column_spacing(20);

    grid.attach(x,1,0,1,1);
    grid.attach(y,2,0,1,1);
    grid.attach(z,3,0,1,1);
    grid.attach(positionLabel,0,1,1,1);
    grid.attach(velocityLabel,0,2,1,1);
    grid.attach(positionX,1,1,1,1);
    grid.attach(positionY,2,1,1,1);
    grid.attach(positionZ,3,1,1,1);
    grid.attach(velocityX,1,2,1,1);
    grid.attach(velocityY,2,2,1,1);
    grid.attach(velocityZ,3,2,1,1);
    show_all_children();
}

//junk destructor
ParticleInfoTab::~ParticleInfoTab()
{
}
