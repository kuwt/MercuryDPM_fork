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

#include "ChangeDomainDialog.h"

void ChangeDomainDialog::getDomainMinMax(Vec3D& domainMin, Vec3D& domainMax)
{
    domainMin.X=ToDouble(domainMinX.get_text());
    domainMin.Y=ToDouble(domainMinY.get_text());
    domainMin.Z=ToDouble(domainMinZ.get_text());
    domainMax.X=ToDouble(domainMaxX.get_text());
    domainMax.Y=ToDouble(domainMaxY.get_text());
    domainMax.Z=ToDouble(domainMaxZ.get_text());                    
}

ChangeDomainDialog::ChangeDomainDialog(const Vec3D& domainMin, const Vec3D& domainMax)
{
    set_title("Change domain");
    add_button("OK",Gtk::RESPONSE_OK);
    add_button("CANCEL",Gtk::RESPONSE_CANCEL);
    domainMinXLabel.set_text("Min X");
    domainMinYLabel.set_text("Min Y");
    domainMinZLabel.set_text("Min Z");
    domainMaxXLabel.set_text("Max X");
    domainMaxYLabel.set_text("Max Y");
    domainMaxZLabel.set_text("Max Z");
    domainMinX.set_max_length(10);
    domainMinY.set_max_length(10);
    domainMinZ.set_max_length(10);
    domainMaxX.set_max_length(10);
    domainMaxY.set_max_length(10);
    domainMaxZ.set_max_length(10);
    domainMinX.set_text(ToString(domainMin.X));
    domainMinY.set_text(ToString(domainMin.Y));
    domainMinZ.set_text(ToString(domainMin.Z));
    domainMaxX.set_text(ToString(domainMax.X));
    domainMaxY.set_text(ToString(domainMax.Y));
    domainMaxZ.set_text(ToString(domainMax.Z));

    get_content_area()->add(grid);

    grid.attach(domainMinXLabel,0,0,1,1);
    grid.attach(domainMinX,     1,0,1,1);
    grid.attach(domainMaxXLabel,2,0,1,1);
    grid.attach(domainMaxX,     3,0,1,1);
    grid.attach(domainMinYLabel,0,1,1,1);
    grid.attach(domainMinY,     1,1,1,1);
    grid.attach(domainMaxYLabel,2,1,1,1);
    grid.attach(domainMaxY,     3,1,1,1);
    grid.attach(domainMinZLabel,0,2,1,1);
    grid.attach(domainMinZ,     1,2,1,1);
    grid.attach(domainMaxZLabel,2,2,1,1);
    grid.attach(domainMaxZ,     3,2,1,1);

    show_all_children();
}
