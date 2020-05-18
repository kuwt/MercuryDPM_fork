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

#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/IntersectionOfWalls.h"

///Particles of a single Species
class ChuteWithContraction : public Chute
{
public:

    ChuteWithContraction(std::string restart_file)
            : Chute()
    {
        std::cout << "ChuteWithContraction" << std::endl;

        loadPeriodicBox(restart_file);

        //dimensions 700x130x10
        //contraction at 300-500x13

        //corrected to dimensions 260x130x4
        //contraction at 60-260x13
        double xContractionStart = 60;
        double xContractionEnd = 140;
        double yContractionStart = 50;
        double yContractionEnd = 25;

        //defines how often the inflow has to be copied
        double lengthPeriodicChute = getXMax(); //length of restarted periodic chute (typically 20); will also be length of inflow
        double widthPeriodicChute = getYMax(); //width of restarted periodic chute (typically 10);
        int numRepetitions = round(xContractionEnd / getXMax()); //length 13x20=260
        int numRepetitionsInWidth = round(yContractionStart / getYMax()); //width 13x10=130

        double sumVel = 0;
        double num = 0;
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            if (!inflowParticle_.isFixed())
            {
                sumVel += inflowParticle_.getVelocity().X;
                num++;
            } // changed P->isfixed()
        setInflowVelocity(sumVel / num);

        //creates enough particles to fill contraction
        //set_Nmax(get_N()*numRepetitionsInWidth*(numRepetitions+2.)); //not needed anymore
        //remove flowing particles
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            if (!inflowParticle_.isFixed())
            {
                particleHandler.removeObject(i);
                i--;
            } // changed P->isFixed()
        //extend contraction (species 1) in width and height
        int N = particleHandler.getNumberOfObjects();
        for (int k = 0; k < numRepetitionsInWidth; k++)
        {
            for (int j = 0; j < numRepetitions; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    if (inflowParticle_.isFixed())
                    {
                        particleHandler.addObject(particleHandler.getObject(i)); //Particles.push_back(Particles[i]);
                        particleHandler.getLastObject()->move(Vec3D(lengthPeriodicChute * j, 0.0, 0.0)); //Particles.back().Position.X+=lengthPeriodicChute*j;
                        particleHandler.getLastObject()->move(Vec3D(0.0, widthPeriodicChute * k, 0.0)); //Particles.back().Position.Y+=widthPeriodicChute*k;
                    }
                }
            }
        }
        std::cout << "Particles added: " << particleHandler.getNumberOfObjects() << std::endl;

        setXMax((numRepetitions + 1) * lengthPeriodicChute);
        setYMax(numRepetitionsInWidth * widthPeriodicChute);

//<<<<<<< .mine
//	///Sets variable values for particles that are created at the inflow
//	void create_inflow_particle()
//	{
//		//~ cout << "create_inflow_particle" << endl;
//        	inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
//		inflowParticle_.computeMass();
//
//                double x,y,z,vx,vy,vz;
//                Vec3D posX,posY,posZ,velX,velY,velZ;
//
//                // X coord of the particle position
//		x = (getXMin() + inflowParticle_.getRadius()*random.getRandomNumber(1.,0.5));
//		posX=Vec3D(x,0.,0.);//P0.Position.X = get_xmin() + P0.Radius*random(1.,5.);
//                // Y coord of the particle position
//		if ((getYMax()-getYMin())==2.0* getMaxInflowParticleRadius())
//		  {
//		    y=getYMin()+ inflowParticle_.getRadius();
//		    posY=Vec3D(0.0,y,0.0);//P0.Position.Y = get_ymin() + P0.Radius;
//                  }
//		else
//                  {
//                    y=random.getRandomNumber(getYMin() + inflowParticle_.getRadius(),getYMax() - inflowParticle_.getRadius());
//                    posY=Vec3D(0.0,y,0.0);//P0.Position.Y = random(get_ymin() + P0.Radius, get_ymax() - P0.Radius);
//                  }
//                //
//                z=random.getRandomNumber(getZMin()+ getFixedParticleRadius() + inflowParticle_.getRadius(),getZMax()- inflowParticle_.getRadius());
//		posZ=Vec3D(0.0,0.0,z);//P0.Position.Z = random(get_zmin() + FixedParticleRadius + P0.Radius, get_zmax() - P0.Radius);
//		//
//                Vec3D pos = posX + posY + posZ;
//		inflowParticle_.setPosition(pos);
//                //
//                // particle velocity
//                vx= getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance()) + getInflowVelocity();
//                vy= getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());
//                vz= getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());
//                //
//                Vec3D vel = Vec3D(vx,vy,vz);
//                inflowParticle_.setVelocity(vel);
//		//P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
//		//P0.Velocity.Y = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance);
//		//P0.Velocity.Z = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance);
//	}
//=======
        PeriodicBoundary* perw = static_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
        perw->set(Vec3D(0.0, 1.0, 1.0), getYMin(), getYMax());
        //set_NWallPeriodic(1);
        //WallsPeriodic[0].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
//>>>>>>> .r819

        //set_NWall(get_NWall()+1);
        //Walls[get_NWall()-1].set(Vec3D( -1.0, 0.0, 0.0), -getXMin());

        //fix_hgrid();
        //set_HGRID_num_buckets_to_power(get_Nmax()); // automated
        //~ Walls.resize(1);
        //~ Walls.back().set(Vec3D(0,0,-1),getInflowParticleRadius());
        set_symmetric_contraction(xContractionStart, xContractionEnd, (yContractionStart - yContractionEnd) / 2.);
        // remove particles in contraction wall
        double dist;
        Vec3D normal;
        bool touch;
        int counter = 0;
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            for (int k = wallHandler.getNumberOfObjects() - 2; k < wallHandler.getNumberOfObjects(); k++)
            {
                touch = wallHandler.getObject(k)->getDistanceAndNormal(*particleHandler.getObject(i), dist, normal);
                if (touch)
                {
                    counter++;
                    particleHandler.removeObject(i);
                    i--;
                    break;
                }
            }
        }
        std::cout << "Particles deleted: " << counter << std::endl;

        setName("ChuteWithStdInflow");

        std::cout << "ChuteWithContraction" << std::endl;
    }

    ///Do not add bottom
    void setup_particles_initial_conditions()
    {
        write(std::cout,false);
        setTimeMax(10000. * getTimeStep());
    }

    ///loads chute data from restart file
    void loadPeriodicBox(std::string const restart_file)
    {
        std::cout << "loadPeriodicBox(" << restart_file << ")" << std::endl;

        //load particles and chute setup into periodic chute
        setName(restart_file.c_str());
        readRestartFile();
        //we don't want to treat the data as restarted, i.e. we start the program at time=0 and create new output files
        setRestarted(false);

        //keep file name but create files in the local directory, i.e. remove folder
        size_t found = restart_file.find_last_of("/\\");
        setName(restart_file.substr(found + 1).c_str());

        std::cout << "loadPeriodicBox end" << std::endl;
    }

    ///Sets variable values for particles that are created at the inflow
    void create_inflow_particle()
    {
        //~ cout << "create_inflow_particle" << endl;
        inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        inflowParticle_.computeMass(Species);

        double x, y, z, vx, vy, vz;
        Vec3D posX, posY, posZ, velX, velY, velZ;

        // X coord of the particle position
        x = (getXMin() + inflowParticle_.getRadius() * random.getRandomNumber(1., 0.5));
        posX = Vec3D(x, 0., 0.); //P0.Position.X = get_xmin() + P0.Radius*random(1.,5.);
        // Y coord of the particle position
        if ((getYMax() - getYMin()) == 2.0 * getMaxInflowParticleRadius())
        {
            y = getYMin() + inflowParticle_.getRadius();
            posY = Vec3D(0.0, y, 0.0); //P0.Position.Y = get_ymin() + P0.Radius;
        }
        else
        {
            y = random.getRandomNumber(getYMin() + inflowParticle_.getRadius(), getYMax() - inflowParticle_.getRadius());
            posY = Vec3D(0.0, y, 0.0); //P0.Position.Y = random(get_ymin() + P0.Radius, get_ymax() - P0.Radius);
        }
        //
        z = random.getRandomNumber(getZMin() + getFixedParticleRadius() + inflowParticle_.getRadius(), getZMax() - inflowParticle_.getRadius());
        posZ = Vec3D(0.0, 0.0, z); //P0.Position.Z = random(get_zmin() + FixedParticleRadius + P0.Radius, get_zmax() - P0.Radius);
        //
        Vec3D pos = posX + posY + posZ;
        inflowParticle_.setPosition(pos);
        //
        // particle velocity
        vx = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance()) + getInflowVelocity();
        vy = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());
        vz = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());
        //
        Vec3D vel = Vec3D(vx, vy, vz);
        inflowParticle_.setVelocity(vel);
        //P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
        //P0.Velocity.Y = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance);
        //P0.Velocity.Z = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance);
    }

    ///Remove particles if they fall below a certain height (allows them to become supercritical)
    void cleanChute()
    {
        //clean outflow every 100 time steps
        static int count = 0, maxcount = 100;
        if (count > maxcount)
        {
            count = 0;
            // delete all outflowing particles
            for (unsigned int i = 0; i < particleHandler.getNumberOfObjects();)
            {
                //~ if (P->Position.Z<get_zmin()-10*getInflowParticleRadius()){
                if (inflowParticle_.getPosition().X > getXMax())
                {
                    //~ cout << "Remove particle" << endl;
                    particleHandler.removeObject(i);
                }
                else
                    i++;
            }
        }
        else
            count++;
    }

    void actions_before_time_step()
    {
        //add_particles(); need to check if it works without this
        cleanChute();
    }

    ///add some particular output
    void cout_time()
    {
    }

    void set_symmetric_contraction(double x_min, double x_max, double delta_y)
    {
        IntersectionOfWalls wall;
        //Walls.resize(Walls.size()+1);
        //back wall
        Vec3D normalIntoWall = Vec3D(-1, 0, 0);
        Vec3D point = Vec3D(x_max, 0, 0);
        wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall, point));
        //slanted wall
        double delta_x = x_max - x_min;
        normalIntoWall = Vec3D(delta_y, -delta_x, 0) / sqrt(mathsFunc::square(delta_x) + mathsFunc::square(delta_y));
        point = Vec3D(x_min, 0, 0);
        wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall, point));

        //Walls.resize(Walls.size()+1);
        //back wall
        normalIntoWall = Vec3D(-1, 0, 0);
        point = Vec3D(x_max, getChuteWidth(), 0);
        wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall, point));
        //slanted wall
        delta_x = x_max - x_min;
        normalIntoWall = Vec3D(delta_y, delta_x, 0) / sqrt(mathsFunc::square(delta_x) + mathsFunc::square(delta_y));
        point = Vec3D(x_min, getChuteWidth(), 0);
        wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall, point));
        wallHandler.copyAndAddObject(wall);
    }

    SphericalParticle inflowParticle_;


};

int main(int argc, char* argv[])
{
    ChuteWithContraction problem("../ini/H10A28L2M0.5B0.5");
    problem.setTimeMax(1e4);
    //problem.setSaveCount(1./problem.getTimeStep());
    problem.setSaveCount(500);
    problem.setXBallsAdditionalArguments("-v0 -solidf -h 600 -w 1400 ");
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.solve(argc, argv);
}
