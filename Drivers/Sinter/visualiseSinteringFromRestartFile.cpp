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

#include<iostream>
#include "Mercury3D.h"
#include "Interactions/NormalForceInteractions/SinterInteraction.h"
using namespace std;

int main(int argc, char *argv[])
{
    if (argc<2) {
        cout << "Please enter the name of the simulation you want to restart and, optionally, the name of the simulation after restart" << endl;
        exit(-1);
    } else {
        cout << "restart data: " << argv[1] << ".restart" << endl;
    }

    Mercury3D dpm;
    dpm.setName(argv[1]);
    dpm.readRestartFile();
    dpm.setRestarted(false);
    dpm.setTimeMax(dpm.getTimeStep());
    dpm.restartFile.setFileType(FileType::NO_FILE);

    logger(INFO,"Writing file %.tex",std::string(argv[1]));
    ofstream os(std::string(argv[1]) + ".tex" );
    os << "\\documentclass[11pt]{standalone}\n"
     "\\usepackage{tikz,graphics,multirow,multicol,tabularx, tabu,xcolor}\n"
     "\n"
     "\\begin{document}\n"
     "\\begin{tikzpicture}[\n"
     "pf/.style ={fill=yellow!50,draw=none},%fill particle\n"
     "pd/.style ={draw=black} %draw particle\n"
     "]\n";

    os << "%draw particles' inside\n";
    for (auto p : dpm.particleHandler)
    {
        const Vec3D& pos = p->getPosition()*1e6;
        const Mdouble& r = p->getRadius()*1e6;
        os << "\\draw[pf] (" + std::to_string(pos.X)
              +"," + std::to_string(pos.Z)
              + ") circle (" + std::to_string(r) +");\n";
    }

    os << "%draw particles' outside\n";
    for (auto p : dpm.particleHandler)
    {
        const Vec3D& pos = p->getPosition()*1e6;
        const Mdouble& r = p->getRadius()*1e6;
        os << "\\draw[pd] (" + std::to_string(pos.X)
              +"," + std::to_string(pos.Z)
              + ") circle (" + std::to_string(r) +");\n";
    }

    os << "%draw plastic overlap' outside\n";
    for (auto i : dpm.interactionHandler)
    {
        auto c = dynamic_cast<SinterInteraction*>(i);
        const Vec3D& cp = c->getContactPoint()*1e6;
        const Vec3D& n = c->getNormal();
        const Vec3D t = Vec3D(n.Z,0,-n.X);
        const Mdouble& po = 0.5*c->getPlasticOverlap()*1e6;
        const Mdouble& o = 0.5*c->getOverlap()*1e6;
        const BaseParticle* p = dynamic_cast<BaseParticle*>(c->getP());
        logger.assert_always(p!= nullptr,"not particle");
        const Mdouble& r = p->getRadius()*1e6;
        const Mdouble x = sqrt(2.0*po*r);
        const Vec3D a0 = cp+x*t+(o-po)*n, a1 = cp+x*t+0.5*r*n, a2 = cp-x*t+0.5*r*n, a3 = cp-x*t+(o-po)*n;
        os << "\\draw[pf] ("+ std::to_string(a0.X) +","+ std::to_string(a0.Z)
              +") -- ("+ std::to_string(a1.X) +","+ std::to_string(a1.Z)
              +") -- ("+ std::to_string(a2.X) +","+ std::to_string(a2.Z)
              +") -- ("+ std::to_string(a3.X) +","+ std::to_string(a3.Z)
              +") -- cycle;\n";
        os << "\\draw[pd] ("+ std::to_string(a0.X) +","+ std::to_string(a0.Z)
              +") -- ("+ std::to_string(a3.X) +","+ std::to_string(a3.Z) + ");\n";
        const Vec3D b0 = cp+x*t-(o-po)*n, b1 = cp+x*t-0.5*r*n, b2 = cp-x*t-0.5*r*n, b3 = cp-x*t-(o-po)*n;
        os << "\\draw[pf] ("+ std::to_string(b0.X) +","+ std::to_string(b0.Z)
              +") -- ("+ std::to_string(b1.X) +","+ std::to_string(b1.Z)
              +") -- ("+ std::to_string(b2.X) +","+ std::to_string(b2.Z)
              +") -- ("+ std::to_string(b3.X) +","+ std::to_string(b3.Z)
              +") -- cycle;\n";
        os << "\\draw[pd] ("+ std::to_string(b0.X) +","+ std::to_string(b0.Z)
              +") -- ("+ std::to_string(b3.X) +","+ std::to_string(b3.Z) + ");\n";
    }
    os << "\\end{tikzpicture}\n"
     "\\end{document}";

    //dpm.setParticlesWriteVTK(true);
//	dpm.setWallsWriteVTK(FileType::MULTIPLE_FILES);
//	dpm.solve(argc-1, argv+1);
    return 0;
}
