//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
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
#include "CG/TimeAveragedCG.h"

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mercury3D cg;

    //cg.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg.readRestartFile("TriolietDietFeeder");
    cg.dataFile.setCounter(0);
    cg.readNextDataFile();
    for (auto& p: cg.particleHandler) {
        if (p->getPosition().Y<0)
            p->setId(0);
        else
            p->setId(1);
    }
    auto p = *cg.particleHandler.getLastObject();
    //logger(ERROR,"% % %",idLeft,idRight,n);

    TimeAveragedCG<CGCoordinates::XYZ,CGFunctions::Heaviside> c0;
    Mdouble h = 0.1;
    c0.setNX(5.078/h);
    c0.setNY(2.658/h);
    c0.setNZ(1.995/h);
    c0.selectSpecies(0);
    //c0->setWidth(cbrt(0.25/pi));//so heaviside should give height 3
    c0.setWidth(0.15);
    //create statistics of the particles which are at x<0 in the restart file
    c0.setSelectedParticle([] (const BaseInteractable* p) { return (p->getId()==0)?false:true; });
    cg.cgHandler.copyAndAddObject(c0);
    //create statistics of the particles which areat x>=0 in the restart file
    c0.setSelectedParticle([] (const BaseInteractable* p) { return (p->getId()==1)?false:true; });
    cg.cgHandler.copyAndAddObject(c0);
//    //create general statistics
//    c0.setSelectedParticle([] (const BaseInteractable* p) { return true; });
//    cg.cgHandler.copyAndAddObject(c0);

    cg.cgHandler.evaluateDataFiles(false);

    //helpers::more(cg.cgHandler.getLastObject()->statFile.getName(),10);

    logger(INFO,"Run 'TriolietDietFeeder3dStatistics.m' in Matlab to view the output");

    helpers::writeToFile("TriolietDietFeeder3dStatistics.m","%% load stat file into matlab\n"
     "raw=importdata('TriolietDietFeeder.0.stat',' ',2);\n"
     "t = raw.data(:,1);\n"
     "x = raw.data(:,2);\n"
     "y = raw.data(:,3);\n"
     "z = raw.data(:,4);\n"
     "volumeFractionLeft = raw.data(:,5);\n"
     "\n"
     "raw=importdata('TriolietDietFeeder.1.stat',' ',2);\n"
     "volumeFractionRight = raw.data(:,5);\n"
     "\n"
     "volumeFraction = volumeFractionLeft + volumeFractionRight;\n"
     "\n"
     "% get data into grid form (this is written general, could be optimized)\n"
     "xUnique = unique(x);\n"
     "yUnique = unique(y);\n"
     "zUnique = unique(z);\n"
     "[X,Y,Z] = meshgrid(xUnique,yUnique,zUnique);\n"
     "VolumeFraction = griddata(x,y,z,volumeFraction,X,Y,Z);\n"
     "VolumeFractionLeft = griddata(x,y,z,volumeFractionLeft,X,Y,Z);\n"
     "VolumeFractionRight = griddata(x,y,z,volumeFractionRight,X,Y,Z);\n"
     "\n"
     "%% plot volume fraction for different z\n"
     "figure(1)\n"
     "iz=round(size(Z,3)/2);\n"
     "subplot(3,1,1); contourf(X(:,:,iz),Y(:,:,iz),VolumeFraction(:,:,iz),'EdgeColor','none'); \n"
     "title('Full stats'); caxis([0 .7]); colorbar; axis image\n"
     "subplot(3,1,2); contourf(X(:,:,iz),Y(:,:,iz),VolumeFractionLeft(:,:,iz),'EdgeColor','none'); \n"
     "title('Left particles'); caxis([0 .7]); colorbar; axis image\n"
     "subplot(3,1,3); contourf(X(:,:,iz),Y(:,:,iz),VolumeFractionRight(:,:,iz),'EdgeColor','none'); \n"
     "title('Right particles'); caxis([0 .7]); colorbar; axis image\n"
     "\n"
     "\n"
     "%% plot volume fraction for different z\n"
     "figure(2)\n"
     "v=VideoWriter('VolumeFraction.avi')\n"
     "v.FrameRate=3\n"
     "open(v)\n"
     "xlabel('x [m]')\n"
     "ylabel('y [m]')\n"
     "caxis([0 .7]) \n"
     "colorbar\n"
     "set(gca,'nextplot','replacechildren')\n"
     "for iz = 1:length(zUnique)\n"
     "    contourf(X(:,:,iz),Y(:,:,iz),VolumeFraction(:,:,iz),'EdgeColor','none')\n"
     "    legend(['z=' num2str(zUnique(iz),2) 'm'])\n"
     "    frame=getframe;\n"
     "    writeVideo(v,frame);\n"
     "end\n"
     "close(v)");
}

