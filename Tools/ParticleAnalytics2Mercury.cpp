//Copyright (c) 2015, The MercuryDPM Developers Team. All rights reserved.
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

#include <iostream>
#include <fstream>
#include <Logger.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearViscoelasticReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include "Mercury3D.h"
using constants::pi;

/*!
 * \brief This gives functionality to read information from binary formats like STL etc.
 * This class is complete stand-alone and is tested with one any reference to other MecuryDPM code except Vections and Logger.
 */
class FileReader
{
public:

    /*!
     * \brief Default constuction, requires to users to prove the name of the file that will be opened.
     */
    explicit FileReader(std::string name) {
        pFile_.open(name+".p3p");
        if (pFile_) {
            version_ = Version::P3;
            wFile_.open(name+".p3w");
            cFile_.open(name+".p3c");
            if (!wFile_) logger(ERROR, "Could not open %", name + ".p3w");
            if (!cFile_) logger(ERROR, "Could not open %", name + ".p3c");
        } else {
            version_ = Version::P4;
            pFile_.open(name+".p4p");
            wFile_.open(name+".p4w");
            cFile_.open(name+".p4c");
            if (!pFile_) logger(ERROR, "Could not open % or %", name + ".p3p",name + ".p4p");
            if (!wFile_) logger(ERROR, "Could not open %", name + ".p4w");
            if (!cFile_) logger(ERROR, "Could not open %", name + ".p4c");
        }
        logger(INFO,"Input files opened");

        dpm.setName(name);
    }

    bool read() {
        dpm.interactionHandler.clear();
        dpm.particleHandler.clear();
        dpm.wallHandler.clear();

        std::string line;

        //read first line p3p
        std::getline(pFile_, line);
        if (!line.compare(""))
            return false;

        //read second line p3p
        Mdouble time;
        unsigned N;
        pFile_ >> time >> N;
        dpm.setTime(time);
        dpm.particleHandler.setStorageCapacity(N);
        std::getline(pFile_, line);

        //read third line p3p
        std::getline(pFile_, line);

        //read next lines p3p
        logger(INFO,"Reading % particles at time %",N,time);
        {
            SphericalParticle p; ///\todo taking this line of of the for loop gives a huge speed improvement; why?
            for (unsigned i = 0; i < N; ++i) {
                unsigned id, species;
                Mdouble volume, mass;
                Vec3D position, velocity;
                pFile_ >> id >> species >> volume >> mass >> position >> velocity;
                std::getline(pFile_, line);
                if (version_==Version::P3) species--;
                //logger(INFO,"% % % % % %",id, species, volume, mass, position, velocity);
                // add new species as necessary
                while (dpm.speciesHandler.getNumberOfObjects() <= species) {
                    LinearViscoelasticReversibleAdhesiveSpecies species;
                    species.setDensity(mass / volume);
                    //use an infinite interaction radius
                    species.setAdhesionForceMax(1e20);
                    species.setAdhesionStiffness(1);
                    logger(INFO, "Adding species of density %", species.getDensity());
                    dpm.speciesHandler.copyAndAddObject(species);
                }
                p.setSpecies(dpm.speciesHandler.getObject(species));
                p.setRadius(cbrt(0.75 / pi * volume));
                p.setPosition(position);
                p.setVelocity(velocity);
                dpm.particleHandler.copyAndAddObject(p)->setId(id);
            }
        }
        //logger(INFO,"Read % particles",dpm.particleHandler.getNumberOfObjects());

        //read first line p3c
        std::getline(cFile_, line);
        if (!line.compare(""))
            return false;

        //read second line p3c
        cFile_ >> time >> N;
        if (fabs(dpm.getTime()/time-1)>0.01)
            logger(ERROR,"Timesteps in p3c and p3p do not agree: % %",dpm.getTime(),time);
        dpm.interactionHandler.setStorageCapacity(N);
        std::getline(cFile_, line);

        //read third line p3c
        std::getline(cFile_, line);

        //read next lines p3c
        logger(INFO,"Reading % contacts",N,time);
        for (unsigned i=0; i<N; ++i) {
            unsigned id1, id2;
            Vec3D force, contact;
            cFile_ >> id1 >> id2;
            BaseParticle* p1 = dpm.particleHandler.getObjectById(id1);
            BaseParticle* p2 = dpm.particleHandler.getObjectById(id2);
            logger.assert(p1!=nullptr,"Particle % does not exist",id1);
            logger.assert(p2!=nullptr,"Particle % does not exist",id1);
            Vec3D P1ToP2 = p2->getPosition()-p1->getPosition();
            BaseInteraction* c = p1->getInteractionWith(p2,time,&dpm.interactionHandler);
            logger.assert(c!= nullptr,"Particle-particle interaction % % does not exist",p1,p2);
            c->setDistance(P1ToP2.getLength());
            c->setNormal(P1ToP2/c->getDistance());
            c->setOverlap(c->getDistance()-p1->getRadius()-p2->getRadius());
            if (version_==Version::P3) {
                cFile_ >> force;
                contact = p1->getPosition()-P1ToP2*((p1->getRadius()-0.5*c->getOverlap())/c->getDistance());
            } else {
                cFile_ >> contact >> force;
            }
            std::getline(cFile_, line);
            c->setContactPoint(contact);
            c->setForce(force);
            if (i%(N/10)==0) {std::cout << "\r " << std::round((double)i/N*100) << '%'; std::cout.flush();}
        }
        std::cout << '\n';

        //read first line p3w
        std::getline(wFile_, line);
        if (!line.compare(""))
            return false;

        //read second line p3w
        wFile_ >> time >> N;
        if (fabs(dpm.getTime()/time-1)>0.01)
            logger(ERROR,"Timesteps in p3w and p3p do not agree");
        dpm.interactionHandler.setStorageCapacity(N);
        std::getline(wFile_, line);

        //read third line p3w
        std::getline(wFile_, line);

        //create wall
        InfiniteWall wall;
        wall.setSpecies(dpm.speciesHandler.getObject(0));
        auto w = dpm.wallHandler.copyAndAddObject(wall);

        //read next lines p3w
        logger(INFO,"Reading % wall contacts",N,time);
        for (unsigned i=0; i<N; ++i) {
            unsigned id;
            Vec3D force, contact, particleToContact;
            wFile_ >> id;
            BaseParticle* p = dpm.particleHandler.getObjectById(id);
            logger.assert(p!=nullptr,"Particle % does not exist",id);
            BaseInteraction* c = w->getInteractionWith(p,time,&dpm.interactionHandler);
            logger.assert(c!= nullptr,"Particle-wall interaction % % does not exist",p,w);
            if (version_==Version::P3) {
                wFile_ >> force >> particleToContact;
                contact = p->getPosition()-particleToContact;
            } else {
                wFile_ >> contact >> force;
                particleToContact = p->getPosition()-contact;
            }
            std::getline(wFile_, line);
            c->setContactPoint(contact);
            c->setDistance(particleToContact.getLength());
            c->setNormal(particleToContact/c->getDistance());
            c->setOverlap(c->getDistance()-p->getRadius());
            c->setForce(force);
        }

        logger(INFO,"Writing output files");
        for (const auto p : dpm.particleHandler) {
            dpm.setMin(Vec3D::min(dpm.getMin(),p->getPosition()));
            dpm.setMax(Vec3D::max(dpm.getMax(),p->getPosition()));
        }
        for (const auto c : dpm.interactionHandler) {
            dpm.setMin(Vec3D::min(dpm.getMin(),c->getContactPoint()));
            dpm.setMax(Vec3D::max(dpm.getMax(),c->getContactPoint()));
        }
        dpm.forceWriteOutputFiles();
        return true;
    }

    /*!
     * \brief Destructor, simple closes the file
     */
    ~FileReader() {
        pFile_.close();
        wFile_.close();
        cFile_.close();
        logger(INFO,"Input files closed");
    }

private:

    /// Pointers for the input files.
    std::ifstream pFile_, cFile_, wFile_;
    /// The version number of the particle analytics files
    enum class Version {P3, P4} version_;

    Mercury3D dpm;
};

int main(int argc, char** argv)
{
    helpers::writeToFile("XRT","xrt");

    //Check to see if we actually received two arguments
    if (argc < 2) {
    //We didn't. Print a usage and exit the program.
    logger(FATAL,"This program converts Particle Analytics (.p3* or p4*) to MercuryDPM files.\n"
                 //"These file can then be used in MercuryCG to analyse your data.\n"
                 "Usage: Call the executable with the base name as argument.\n"
                 "E.g. to convert name.p* call\n"
                 "  ./ParticleAnalytics2Mercury name\n", argv[0]);
    }
    FileReader fileReader(argv[1]);
    while (fileReader.read());

}
