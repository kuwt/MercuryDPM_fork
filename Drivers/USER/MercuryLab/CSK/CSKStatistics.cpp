#include "CSK.h"
#include <sstream>
#include <chrono>
#include <string>
#include <iostream>
#include <CG/CG.h>
#include <CG/TimeAveragedCG.h>

/**
 * Analyses the mixing stage
 */

int main()
{
    //Load a specific case study by changing the name variable
    std::string name = "2earlyMixRestartTest-rpm36";
    CSK csk(name);
    csk.setRPM(36);//not sure if needed
    csk.setName(name);

    //Create spatially averaged statistics
    auto cg0 = csk.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg0->statFile.setName(name+".O.stat");

    //Create time- and y-averaged statistics for all particles in the slice |y|<0.05
    auto cg1 = csk.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XZ,CGFunctions::Heaviside>());
    cg1->setMin({-.9,-.05,0});
    cg1->setMax({ .9, .05,0.6});
    cg1->setNX(100);
    cg1->setNZ(50);
    cg1->setWidth(0.02);
    cg1->setSelectedParticle([] (const BaseInteractable* p) { return fabs(p->getPosition().Y)<0.05; });
    cg1->statFile.setName(name+".XZ.stat");

    //Create time- and z-averaged statistics for all particles in the slice |z-0.4|<0.05
    auto cg2 = csk.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XY,CGFunctions::Heaviside>());
    cg2->setMin({-.9,-.9,0.35});
    cg2->setMax({ .9, .9,0.45});
    cg2->setNX(70);
    cg2->setNY(70);
    cg2->setWidth(0.02);
    cg2->setSelectedParticle([] (const BaseInteractable* p) { return fabs(p->getPosition().Z-0.4)<0.05; });
    cg2->statFile.setName(name+".XY.stat");

    //Create time- and y-averaged statistics for species-1 particles in the slice |y|<0.05
    auto cg3 = csk.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XZ,CGFunctions::Heaviside>());
    cg3->setMin({-.9,-.05,0});
    cg3->setMax({ .9, .05,0.6});
    cg3->setNX(100);
    cg3->setNZ(50);
    cg3->setWidth(0.02);
    cg3->setSelectedParticle([] (const BaseInteractable* p) { return (p->getIndSpecies()==0) && fabs(p->getPosition().Y)<0.05; });
    cg3->statFile.setName(name+".XZ.Species0.stat");

    //Set the time interval on which all statistics are evaluated (by default, no timeMax is set).
    //To test-run, choose a small time interval.
    for (auto cg: csk.cgHandler)
    {
        cg->setTimeMin(2);
        //cg->setTimeMax(8);
    }

    //Now do the evaluation
    csk.cgHandler.evaluateDataFiles(false);

    //Write a gnuplot script to view the spatially-avaraged output
    helpers::writeToFile(name+".O.gnu", "p '2earlyMixRestartTest-rpm36.O.stat' u 1:5 w l");

    //Write a matlab file to view the other output
    helpers::writeToFile(name+".m","%% load time- and y-averaged statistics for all particles in the slice |y|<0.05\n"
                          "\n"
                          "raw=importdata('2earlyMixRestartTest-rpm36.XZ.stat',' ',2);\n"
                          "t = raw.data(:,1);\n"
                          "x = raw.data(:,2);\n"
                          "z = raw.data(:,3);\n"
                          "volFraction = raw.data(:,4);\n"
                          "rho = raw.data(:,5);\n"
                          "velX = raw.data(:,6);\n"
                          "velY = raw.data(:,7);\n"
                          "velZ = raw.data(:,8);\n"
                          "% get data into grid form (this is written general, could be optimized)\n"
                          "xUnique = unique(x);\n"
                          "zUnique = unique(z);\n"
                          "[X,Z] = meshgrid(xUnique,zUnique);\n"
                          "volumeFraction = griddata(x,z,volFraction,X,Z);\n"
                          "density = griddata(x,z,rho,X,Z);\n"
                          "velocityX = griddata(x,z,velX,X,Z);\n"
                          "velocityY = griddata(x,z,velY,X,Z);\n"
                          "velocityZ = griddata(x,z,velZ,X,Z);\n"
                          "speed = sqrt(velocityX.^2+velocityY.^2+velocityZ.^2);\n"
                          "\n"
                          "%% plot density and overlay with momentum field \n"
                          "% note, momentum is better than velocity as it shows mass transfer rate)\n"
                          "% note, why is the volume fraction so high at teh base of the flow\n"
                          "figure(1)\n"
                          "h=pcolor(X,Z,density); \n"
                          "set(h,'EdgeColor','none')\n"
                          "xlabel('x')\n"
                          "ylabel('z')\n"
                          "%caxis([0 1]); \n"
                          "colormap((1:-.01:0.5)'*[1 1 1])\n"
                          "colorbar; axis image\n"
                          "ylim([0 .6])\n"
                          "xlim([-.9 .9])\n"
                          "hold on\n"
                          "quiver(X,Z,volumeFraction.*velocityX,volumeFraction.*velocityZ,5); \n"
                          "hold off\n"
                          "saveas(gcf,'VerticalConvection.jpg')\n"
                          "\n"
                          "%% load time- and y-averaged statistics for species-0 particles in the slice |y|<0.05\n"
                          "\n"
                          "raw=importdata('2earlyMixRestartTest-rpm36.XZ.Species0.stat',' ',2);\n"
                          "volFraction0 = raw.data(:,4);\n"
                          "rho0 = raw.data(:,5);\n"
                          "velX0 = raw.data(:,6);\n"
                          "velY0 = raw.data(:,7);\n"
                          "velZ0 = raw.data(:,8);\n"
                          "volumeFraction0 = griddata(x,z,volFraction0,X,Z);\n"
                          "density0 = griddata(x,z,rho0,X,Z);\n"
                          "velocityX0 = griddata(x,z,velX0,X,Z);\n"
                          "velocityY0 = griddata(x,z,velY0,X,Z);\n"
                          "velocityZ0 = griddata(x,z,velZ0,X,Z);\n"
                          "speed0 = sqrt(velocityX0.^2+velocityY0.^2+velocityZ0.^2);\n"
                          "\n"
                          "%% plot density of species-0 \n"
                          "% note, momentum is better than velocity as it shows mass transfer rate)\n"
                          "% note, why is the volume fraction so high at teh base of the flow\n"
                          "figure(2)\n"
                          "h=pcolor(X,Z,volumeFraction0); \n"
                          "% uncomment to plot ratio of species-0 volFraction vs total volFraction\n"
                          "% h=pcolor(X,Z,volumeFraction0./volumeFraction); \n"
                          "set(h,'EdgeColor','none')\n"
                          "xlabel('x')\n"
                          "ylabel('z')\n"
                          "colormap(jet)\n"
                          "colorbar; axis image\n"
                          "saveas(gcf,'Index.jpg')\n"
                          "\n"
                          "\n"
                          "\n"
                          "\n"
                          "%% load time- and z-averaged statistics for all particles in the slice |z-0.4|<0.05\n"
                          "\n"
                          "raw=importdata('2earlyMixRestartTest-rpm36.XY.stat',' ',2);\n"
                          "t = raw.data(:,1);\n"
                          "x = raw.data(:,2);\n"
                          "y = raw.data(:,3);\n"
                          "volFraction = raw.data(:,4);\n"
                          "rho = raw.data(:,5);\n"
                          "velX = raw.data(:,6);\n"
                          "velY = raw.data(:,7);\n"
                          "velZ = raw.data(:,8);\n"
                          "% get data into grid form (this is written general, could be optimized)\n"
                          "xUnique = unique(x);\n"
                          "yUnique = unique(y);\n"
                          "[X,Y] = meshgrid(xUnique,yUnique);\n"
                          "volumeFraction = griddata(x,y,volFraction,X,Y);\n"
                          "density = griddata(x,y,rho,X,Y);\n"
                          "velocityX = griddata(x,y,velX,X,Y);\n"
                          "velocityY = griddata(x,y,velY,X,Y);\n"
                          "velocityZ = griddata(x,y,velZ,X,Y);\n"
                          "speed = sqrt(velocityX.^2+velocityY.^2+velocityZ.^2);\n"
                          "\n"
                          "%% plot density and overlay with momentum field \n"
                          "% note, momentum is better than velocity as it shows mass transfer rate)\n"
                          "% note, why is the volume fraction so high at teh base of the flow\n"
                          "figure(3)\n"
                          "h=pcolor(X,Y,density); \n"
                          "set(h,'EdgeColor','none')\n"
                          "xlabel('x')\n"
                          "ylabel('y')\n"
                          "colormap((1:-.01:0.5)'*[1 1 1])\n"
                          "colorbar; axis image\n"
                          "ylim([-.9 .9])\n"
                          "xlim([-.9 .9])\n"
                          "hold on\n"
                          "quiver(X,Y,velocityX,velocityY,1); \n"
                          "XOut = X\n"
                          "XOut(X.^2+Y.^2<0.7)=nan\n"
                          "quiver(XOut,Y,-sqrt(X.^2+Y.^2)*36/60.*Y,sqrt(X.^2+Y.^2)*36/60.*X,1); \n"
                          "hold off\n"
                          "saveas(gcf,'HorizontalConvection.jpg')\n"
                          ""
     );

    return 0;
}
