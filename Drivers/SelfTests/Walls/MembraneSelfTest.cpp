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
#include <Mercury3D.h>
#include <Species/HertzianViscoelasticMindlinSpecies.h>
#include <Walls/Membrane/Membrane.h>
#include <CMakeDefinitions.h>

#include <climits>

/*
* This is a Demo for the Class Membrane utilizing the MeshTriangle walls.
* In this demo, we have a square elastic membrane fixed at its edges. We put 
* a few particles on it and simulate until the system is relaxed.
*/
class MembraneSelfTest : public Mercury3D
{
public:
    
     void setupInitialConditions() override
    {
        // Setting the dimensions of the simulation box
        setXMin(-1);
        setYMin(-1);
        setZMin(-1);
        setXMax(1);
        setYMax(1);
        setZMax(1);
        
        // Initialize the species
        initializeSpecies();
        
        // Create the membrane
        setUpMembrane(membrane_);
        
        // Place a few particles
        SphericalParticle p0;
        p0.setRadius(particleRadius_);
        p0.setSpecies(particleSpecies_);
        
        
        Vec3D pos = Vec3D(0.5, 1, 1.2)*p0.getRadius();
        p0.setPosition(pos);
        particleHandler.copyAndAddObject(p0);
        
        setRotation(1);
    }
    
    /*!
     * \details This function initializes the species needed.
     */
    void initializeSpecies()
    {
        // Declare a base species
        HertzianViscoelasticMindlinSpecies species;

        // Initiallize the particle species
        // For Mindlin
        species.setEffectiveElasticModulusAndPoissonRatio(particleElasticModulus_, particlePoissonRatio_);
        species.setEffectiveElasticModulusAndRestitutionCoefficient(particleElasticModulus_, particleCor_);
        species.setDensity(particleDensity_);
        species.setSlidingFrictionCoefficient(particleSlidingFrictionCoefficient_);
        species.setSlidingDissipation(species.getDissipation());
        particleSpecies_ = speciesHandler.copyAndAddObject(species);

        logger(INFO, "Particle material properties:\n    Young's Modulus [Pa] %\n    Poisson's ratio %\n    COR % \n    Friction Coefficient %",
                    particleElasticModulus_, particlePoissonRatio_, particleCor_, particleSlidingFrictionCoefficient_);

        // Initialize the membrane species
        // For Mindlin
        species.setEffectiveElasticModulusAndPoissonRatio(membraneElasticModulus_, membranePoissonRatio_);
        species.setEffectiveElasticModulusAndRestitutionCoefficient(membraneElasticModulus_, membraneCor_);
        species.setDensity(membraneDensity_);
        species.setSlidingFrictionCoefficient(membraneSlidingFrictionCoefficient_);
        species.setSlidingDissipation(species.getDissipation());
        membraneSpecies_ = speciesHandler.copyAndAddObject(species);
        membraneParticleSpecies_ = speciesHandler.copyAndAddObject(species);
        
        logger(INFO, "Membrane material properties:\n    Young's Modulus [Pa] %\n    Poisson's ratio %\n    COR % \n    Friction Coefficient %",
                    membraneElasticModulus_, membranePoissonRatio_, membraneCor_, membraneSlidingFrictionCoefficient_);


        auto mixedSpecies = speciesHandler.getMixedObject(particleSpecies_, membraneParticleSpecies_);
        
        // Deactivate interactions between the membrane particles and other particles
        // For Mindlin
        mixedSpecies->setEffectiveElasticModulus(0);
        mixedSpecies->setSlidingFrictionCoefficient(0);
            
        
    }
    
    /*!
     * \details This function creates a membrane from the supplied stl file.
     */
    void setUpMembrane(Membrane &membrane)
    {
        // Create the membrane and set some importan attributes
        membrane = Membrane(membraneSpecies_, membraneParticleSpecies_);
        membrane.setElasticModulusAndThickness(membraneElasticModulus_, thickness_);
        membrane.setCriticalDampingCoefficient(critDampCoeff_);
        membrane.setKeAndKd(Ke_, Kd_);
        
        membrane.setDPMBase(this);
        
        membraneParticleRadius_ = particleRadius_*0.2;
        
        membrane.setParticleRadius(membraneParticleRadius_);
        
        logger(INFO, "membraneParticleRadius %", membraneParticleRadius_);
        logger(INFO, "Set membrane Kn=%, xi=%, Ke_=%, Kd_=%", 
               membrane.getKn(), membrane.getCriticalDampingCoefficient(),
               membrane.getKe(), membrane.getKd());

        // Construct the Mercury particles using the Icosahedron
        SphericalParticle p0;

        // Set initial attributes
        p0.setSpecies(membraneParticleSpecies_);
        p0.setRadius(membrane.getParticleRadius());
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        membrane.loadFromSTL(p0, getMercuryDPMSourceDir() + "/Drivers/SelfTests/Walls/MembraneSelfTest.stl");
                    
        fixMembraneEdges();
    }
    
    /*!
     * \details this function fixes the particles at the edges of the membrane.
     */
    void fixMembraneEdges()
    {
        std::vector<BaseParticle*> vertexParticles = membrane_.getVertexParticles();
        
        Mdouble eps = 1e-8;
        Vec3D pos;
        
        Mdouble minY=getMax().Y;
        Mdouble maxY=getMin().Y;
        for (auto p0: vertexParticles)
        {
            pos = p0->getPosition();
            minY = std::min(pos.Y, minY);
            maxY = std::max(pos.Y, maxY);
        }
        
        for (auto p0: vertexParticles)
        {
            pos = p0->getPosition();
        
            if ( std::fabs(pos.Y-minY) < eps || std::fabs(pos.Y-maxY) < eps)
                 {
                     p0->fixParticle();
                     logger(INFO, "Fized particle at %", pos);
                 }
        }
        
        // membrane_.updateEdgeMass();
    }
    
    void computeAdditionalForces() override
    {
        // Apply the forces due to the mass spring system
        membrane_.computeAdditionalForces();
    }
    
    void write(std::ostream& os, bool writeAllParticles) const override
    {
        Mercury3D::write(os, writeAllParticles);
        os << " membrane " << membrane_;
        
        os << " particleSpecies " << particleSpecies_->getId();
        os << " membraneSpecies " << membraneSpecies_->getId();
        os << " membraneParticleSpecies " << membraneParticleSpecies_->getId();
    }
    
    void read(std::istream& is, ReadOptions opt) override
    {
        
        std::string dummy;
        unsigned int Id;
        
        Mercury3D::read(is, opt);
        
        membrane_.setDPMBase(this);
        is >> dummy >> membrane_;
        
        is >> dummy >> Id;
        particleSpecies_ = dynamic_cast<HertzianViscoelasticMindlinSpecies*>(speciesHandler.getObjectById(Id));
        is >> dummy >> Id;
        membraneSpecies_ = dynamic_cast<HertzianViscoelasticMindlinSpecies*>(speciesHandler.getObjectById(Id));
        is >> dummy >> Id;
        membraneParticleSpecies_ = dynamic_cast<HertzianViscoelasticMindlinSpecies*>(speciesHandler.getObjectById(Id));
    }
    
private:
    // Membrane Definitions
    Mdouble membraneElasticModulus_ = 1.25e6;
    Mdouble membranePoissonRatio_ = 1/3.0; // Has to be 1/3.0 as the mass spring system is constrained to that
    Mdouble membraneCor_ = 0.66;
    Mdouble membraneSlidingFrictionCoefficient_ = 1.2;
    Mdouble membraneDensity_ = 2000;
    Mdouble Ke_ = 1e-3;
    Mdouble Kd_ = 1e-4;
    Mdouble critDampCoeff_ = 0.15;
    
    Mdouble thickness_ = 0.0003;
    
    Mdouble membraneParticleRadius_;
    
    // Particle Definitions
    Mdouble particleRadius_ = 4e-3;
    Mdouble particleElasticModulus_ = 1e7;
    Mdouble particlePoissonRatio_ = 0.245;
    Mdouble particleCor_ = 0.1;
    Mdouble particleSlidingFrictionCoefficient_ = 0.5;
    Mdouble particleDensity_ = 2000;
    
    // Pointers to the species.
    HertzianViscoelasticMindlinSpecies* particleSpecies_ = nullptr;
    HertzianViscoelasticMindlinSpecies* membraneSpecies_ = nullptr;
    HertzianViscoelasticMindlinSpecies* membraneParticleSpecies_ = nullptr;
    
    Membrane membrane_;

};

int main(int argc, char** argv)
{
  MembraneSelfTest problem;


  problem.setTimeMax(0.2);
  problem.setGravity(Vec3D(0.0,0.0,-9.81));

  /* Make sure to add an (unique) name here */
  problem.setName("MembraneSelfTest");

  // problem.setFileType(FileType::ONE_FILE);
  
  problem.dataFile.setFileType(FileType::NO_FILE);
  problem.restartFile.setFileType(FileType::ONE_FILE);
  problem.fStatFile.setFileType(FileType::NO_FILE);
  problem.eneFile.setFileType(FileType::NO_FILE);

  //Add any configuration which is modified often below...
  problem.setTimeStep(1e-5);

  problem.setSaveCount(floor(1.0/problem.getTimeStep()));
  
  if (helpers::readFromCommandLine(argc, argv, "-vtk"))
  {
      helpers::removeFromCommandline(argc, argv, "-vtk", 0);
      problem.setSaveCount(floor(0.01/problem.getTimeStep()));
      problem.setParticlesWriteVTK(true);
      // Write one file per timestep
      problem.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
  }
  
  problem.solve(argc,argv);
  return 0;
}
