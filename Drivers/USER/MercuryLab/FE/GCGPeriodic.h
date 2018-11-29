#ifndef MERCURY_GCGPERIODIC_H
#define MERCURY_GCGPERIODIC_H

#include <random>
#include <Walls/Screw.h>
#include <Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include "Mercury3D.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/BasicIntersectionOfWalls.h"
#include "addSpecies.h"

using constants::pi;
using mathsFunc::cubic;

class GCGPeriodic : public Mercury3D
{
public:
    /**
     * Simulates only one section of the GCG blade; two derived classes exist that add extra functionality:
     *  - GCGBlender, which  simulates the full blender with in- and outflow
     *  - GCGBlenderWithBelt, which adds a belt beneath the outflow
     *
     * All defining parameters are passed into the constructor.
     *
     * The dimensions of the system are:
     *  - |y|,|z| < casingRadius
     *  - |x| < bladeDistance,
     *
     * The geometry is split into several walls:
     *  - casing: cylinder forming the outer walls
     *    parameters:
     *     - casingRadius: inner radius of the casing
     *     - casingThickness: thickness of the casing (only needed for visualisation)
     *  - shaft: cylinder forming the inner shaft onto which the blades are attached
     *    parameters:
     *    - shaftRadius: outer radius of the shaft
     *  - blade: four blades that will be positioned on the shaft; their shape is a flat plate shaped like a slice of a circle
     *    parameters:
     *    - bladeRadius: outer radius of the blades
     *    - bladeDistance: distance between the blades (blades will be attached at x = +/- bladeDistance/2)
     *    - bladeInclination: angle by which the blade normal is inclined with respect to the x-axis
     *    - bladeOpening: arc angle of the slice
     *    - bladeThickness: width of the blade
     *  The values of these parameters are fixed and not part of the constructor arguments.
     *
     *  @param particlesType speciesType describing the material/contact properties of the particle interaction
     *  @param wallType      speciesType describing the contact properties of the particle-wall interaction
     *  @param solidFraction     volume fraction of particles in the mixer: volume of particles/(volume inside casing - volume inside shaft and blades)
     *  @param rpm             rotation rate of the mixer (in revolutions per minute)
     *  @param particleRadius  mean particle radius (in m)
     *  @param particleRadius2 mean particle radius of species 2 (in m)
     *  @param volFracSpec1 volume fraction of species 1
     *  @param polydispersity  relative standard deviation of the particle radius, assuming normal distribution for species 1
     *  @param polydispersity2  relative standard deviation of the particle radius2, assuming normal distribution for species 2
     */
    GCGPeriodic(SpeciesType fillerType, SpeciesType apiType, SpeciesType wallType, Mdouble solidFraction, Mdouble rpm, Mdouble particleRadius, Mdouble particleRadius2, Mdouble volFracSpec1, Mdouble polydispersity, Mdouble polydispersity2, Mdouble shaftRadius)
    : fillerType(fillerType), apiType(apiType), solidFraction_(solidFraction), particleRadius_(particleRadius), particleRadius2_(particleRadius2), volFracSpec1_(volFracSpec1), polydispersity_(polydispersity), polydispersity2_(polydispersity2), shaftRadius(shaftRadius)
    {
        logger(INFO, "Mean particle radius species 1: % +/- % m", particleRadius, polydispersity * particleRadius);
        logger(INFO, "Mean particle radius species 2: % +/- % m", particleRadius2,
               polydispersity2 * particleRadius2); /// ADDED by M
    
        //set name, gravity, dimensions
        setGeometry();
        setName("GCGPeriodic");
        setGravity({0,0,-9.8});
        addSpecies();
        addWalls(rpm);
    }
    
    void addSpecies () {
        //add species
        logger(INFO,"Adding four species for filler (0), api (1), walls (2) and belt (3)");
        //set the filler species
        Mdouble minParticleRadius1, minParticleRadius2;
        if (hasPSD(fillerType)) {
            particleRadius_ = getMedianParticleRadius(getPSD(fillerType));
            minParticleRadius1 = getMinParticleRadius(getPSD(fillerType));
        } else {
            minParticleRadius1 = (1. - 1.5 * polydispersity_) * particleRadius_;
        }
        auto filler = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(addSingleSpecies(fillerType, particleRadius_, minParticleRadius1, volFracSpec1_!=0));
        //set the api species
        if (hasPSD(apiType)) {
            particleRadius2_ = getMedianParticleRadius(getPSD(apiType));
            minParticleRadius2 = getMinParticleRadius(getPSD(apiType));
        } else {
            minParticleRadius2 = particleRadius2_*(1-1.5*polydispersity2_);
        }
        auto api = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(addSingleSpecies(apiType, particleRadius2_, minParticleRadius2, volFracSpec1_!=1));
        //add mixed property
        speciesHandler.getMixedObject(0, 1)->mixAll(filler, api);
        //set the species of the walls
        auto wall = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.copyAndAddObject(filler));
        auto fillerWall = speciesHandler.getMixedObject(filler,wall);
        auto apiWall = speciesHandler.getMixedObject(api,wall);
        fillerWall->mixAll(filler, filler);
        fillerWall->setSlidingFrictionCoefficient(0.44);
        fillerWall->setRollingFrictionCoefficient(0.11);
        fillerWall->setCohesionStiffness(0.15*fillerWall->getUnloadingStiffnessMax());
        apiWall->mixAll(api, api);
        apiWall->setSlidingFrictionCoefficient(0.44);
        apiWall->setRollingFrictionCoefficient(0.11);
        apiWall->setCohesionStiffness(0.15*apiWall->getUnloadingStiffnessMax());

        //set the species of the belt
        auto belt = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.copyAndAddObject(filler));
        auto fillerBelt = speciesHandler.getMixedObject(filler,belt);
        auto apiBelt = speciesHandler.getMixedObject(api,belt);
        fillerBelt->mixAll(filler, filler);
        fillerBelt->setRestitutionCoefficient(0.1, filler->getMassFromRadius(minParticleRadius1));
        apiBelt->mixAll(api, api);
        apiBelt->setRestitutionCoefficient(0.1, api->getMassFromRadius(minParticleRadius2));
    }

    
    void addWalls(Mdouble rpm) {
        // angular velocity of axis, screw and blades
        const Vec3D angularVelocity = {2. * pi * rpm / 60., 0, 0};

        //set name, gravity, dimensions
        setName("GCGPeriodic");
        setGravity({0, 0, -9.8});
        setMax({bladeDistance, casingRadius, casingRadius});
        setMin(-getMax());
        setXBallsAdditionalArguments("-s 5.5 -oh 0 -o 100 -h 530 -w 1500");
    
        //set the species of the walls
        logger.assert(speciesHandler.getNumberOfObjects()>=2,"Number of species too small");
        auto wSpecies = speciesHandler.getObject(1);

        //save every 0.02s
        setSaveCount(0.02 / getTimeStep());

        // add casing
        AxisymmetricIntersectionOfWalls casing;
        casing.setSpecies(wSpecies);
        casing.setPosition({0, 0, 0});
        casing.setAxis({1, 0, 0});
        casing.addObject({1, 0, 0}, {casingRadius, 0, 0});
        AxisymmetricIntersectionOfWalls casingRender(casing);
        casingRender.addObject({-1, 0, 0}, {casingRadius + casingThickness, 0, 0});
        casingRender.addObject({0, 0, -1}, {0, 0, bladeDistance});
        casingRender.addObject({0, 0, 1}, {0, 0, -bladeDistance});
        casingRender.clear(); //uncomment to make casing invisible
        casing.addRenderedWall(casingRender.copy());
        wallHandler.copyAndAddObject(casing);

        // add shaft
        AxisymmetricIntersectionOfWalls shaft;
        shaft.setSpecies(wSpecies);
        shaft.setPosition({0, 0, 0});
        shaft.setAxis({1, 0, 0});
        shaft.addObject({-1, 0, 0}, {shaftRadius, 0, 0});
        shaft.setAngularVelocity(angularVelocity);
        AxisymmetricIntersectionOfWalls shaftRender(shaft);
        shaftRender.addObject({0, 0, -1}, {0, 0, bladeDistance});
        shaftRender.addObject({0, 0, 1}, {0, 0, -bladeDistance});
        shaft.addRenderedWall(shaftRender.copy());
        wallHandler.copyAndAddObject(shaft);
    
        addBladePair(bladeInclination, {-bladeDistance / 2, 0, 0}, {0, 0, 0}, angularVelocity, wSpecies);
        addBladePair(-bladeInclination, {bladeDistance / 2, 0, 0}, {2, 0, 0}, angularVelocity, wSpecies);

        logger(INFO,"\nAdded % walls", wallHandler.getNumberOfObjects());

        PeriodicBoundary p;
        p.set({1,0,0},getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(p);
    }
    
    void writeParaviewSnapshotScript () {
        helpers::writeToFile(getName()+".snapshot.py",
                             "from paraview.simple import *\n"
                             "import "+getName()+"\n"
                                                 "\n"
                                                 "# get last timestep\n"
                                                 "animationScene1 = GetAnimationScene()\n"
                                                 "animationScene1.GoToLast()\n"
                                                 "\n"
                                                 "# Turn camera into sideway view\n"
                                                 "# https://www.paraview.org/Wiki/ParaView_and_Python#Control_the_camera\n"
                                                 "camera=GetActiveCamera()\n"
                                                 "camera.Elevation(-90)\n"
                                                 "\n"
                                                 "# get active view\n"
                                                 "renderView1 = GetActiveViewOrCreate('RenderView')\n"
                                                 "# make background white\n"
                                                 "renderView1.Background = [1.0, 1.0, 1.0]\n"
                                                 "\n"
                                                 "# if using pvpython, use Interact() to keep the view open\n"
                                                 "Interact()\n"
                                                 "# save screenshot\n"
                                                 "SaveScreenshot('"+getName()+".png', magnification=4, quality=100, view=renderView1)"
        );
    }

protected:
    
    /**
     * Adds a pair of blades with
     *  - bladeInclination; inclination of the blade
     *  - position: the position of the blade center
     *  - rotation: the angular position of the blade
     *  - angularVelocity: rotation speed of the blade
     *  - species: species of the blade
     */
    void addBladePair(Mdouble bladeInclination, Vec3D position, Vec3D rotation, Vec3D angularVelocity,
                      ParticleSpecies* species)
    {
        //logger(WARN,"bi % pos % rot % avel % bo %",bladeInclination, position, rotation, angularVelocity, bladeOpening);
        //GCG bi 0.349066 pos -0.015 0 0 rot 0 0 0 avel 12.5664 0 0
        // add blade
        IntersectionOfWalls bladeSlice;
        bladeSlice.setSpecies(species);
        const Mdouble ci = cos(bladeInclination*pi/180.0);
        const Mdouble si = sin(bladeInclination*pi/180.0);
        const Mdouble co = cos(bladeOpening*pi/180.0 / 2);
        const Mdouble so = sin(bladeOpening*pi/180.0 / 2);
        Matrix3D rotI = {ci, si, 0, -si, ci, 0, 0, 0, 1};
        Matrix3D rotO = {1, 0, 0, 0, co, -so, 0, so, co};
        Matrix3D rotOinv = {1, 0, 0, 0, co, so, 0, -so, co};
        bladeSlice.addObject(rotI * (rotO * Vec3D(0, 1, 0)), {0, 0, 0}); //narrow
        bladeSlice.addObject(rotI * (rotOinv * Vec3D(0, -1, 0)), {0, 0, 0});
        bladeSlice.addObject(rotI * Vec3D(-1, 0, 0), rotI * 0.5 * bladeThickness * Vec3D(1, 0, 0)); //wide
        bladeSlice.addObject(rotI * Vec3D(1, 0, 0), rotI * 0.5 * bladeThickness * Vec3D(-1, 0, 0));
        AxisymmetricIntersectionOfWalls bladeAngle;
        bladeAngle.setSpecies(species);
        bladeAngle.setAxis({ci, -si, 0});
        bladeAngle.addObject({-1, 0, 0}, {bladeRadius, 0, 0});
        BasicIntersectionOfWalls blade;
        blade.setSpecies(species);
        blade.add(bladeSlice);
        blade.add(bladeAngle);
        blade.setAngularVelocity(angularVelocity);
        IntersectionOfWalls bladeRender(bladeSlice);
        for (int i = -25; i <= 25; i++)
        {
            const Mdouble a = bladeOpening * i / 50.0;
            const Mdouble co = cos(a*pi/180.0);
            const Mdouble so = sin(a*pi/180.0);
            Matrix3D rotO = {1, 0, 0, 0, co, -so, 0, so, co};
            const Vec3D n = rotI * (rotO * Vec3D(0, 0, -1));
            bladeRender.addObject(n, -bladeRadius * n);
        }
        //bladeRender.setAngularVelocity(angularVelocity);
        auto bladeLeft = wallHandler.copyAndAddObject(blade);
        bladeRender.setPrescribedPosition([position](Mdouble time){return position;});
        bladeRender.setPrescribedOrientation([bladeLeft](Mdouble time){return bladeLeft->getOrientation();});
        bladeLeft->addRenderedWall(bladeRender.copy());
        bladeLeft->setPosition(position);
        bladeLeft->rotate(rotation);
        auto bladeRight = wallHandler.copyAndAddObject(blade);
        bladeRender.setPrescribedOrientation([bladeRight](Mdouble time){return bladeRight->getOrientation();});
        bladeRight->addRenderedWall(bladeRender.copy());
        bladeRight->setPosition(position);
        // rotate(quarterturn) rotates object by 90 degrees
        const Vec3D quarterTurn = {2, 0, 0};
        bladeRight->rotate(rotation);
        bladeRight->rotate(quarterTurn);
        bladeRight->rotate(quarterTurn);
    }
    
    void addBlade(Mdouble bladeInclination, Vec3D position, unsigned rotation, Vec3D angularVelocity,
                      ParticleSpecies* species)
    {
        // check if blade is a pin
        if (std::isnan(bladeInclination)) {
            AxisymmetricIntersectionOfWalls w;
            w.setSpecies(species);
            //w.setAngularVelocity(angularVelocity);
            w.setOrientationViaNormal(Vec3D(0,0,1));
            w.addObject(Vec3D(-1,0,0), Vec3D(bladeThickness,0,0));
            w.addObject(Vec3D(0,0,-1), Vec3D(0,0,bladeRadius-shaftRadius-casingRadius));
            w.addObject(Vec3D(0,0,1), Vec3D(0,0,-casingRadius));
            w.setPosition(position);
            wallHandler.copyAndAddObject(w);
            return;
        }

        //uncomment to aid visualisation of forward-backward pattern
        //if (bladeInclination>0) bladeThickness=10e-3;
        //else bladeThickness = 2e-3;

        // add blade
        IntersectionOfWalls bladeSlice;
        bladeSlice.setSpecies(species);
        const Mdouble ci = cos(bladeInclination*pi/180.0);
        const Mdouble si = sin(bladeInclination*pi/180.0);
        const Mdouble co = cos(bladeOpening*pi/180.0 / 2);
        const Mdouble so = sin(bladeOpening*pi/180.0 / 2);
        Matrix3D rotI = {ci, si, 0, -si, ci, 0, 0, 0, 1};
        Matrix3D rotO = {1, 0, 0, 0, co, -so, 0, so, co};
        Matrix3D rotOinv = {1, 0, 0, 0, co, so, 0, -so, co};
        bladeSlice.addObject(rotI * (rotO * Vec3D(0, 1, 0)), {0, 0, 0}); //narrow
        bladeSlice.addObject(rotI * (rotOinv * Vec3D(0, -1, 0)), {0, 0, 0});
        bladeSlice.addObject(rotI * Vec3D(-1, 0, 0), rotI * 0.5 * bladeThickness * Vec3D(1, 0, 0)); //wide
        bladeSlice.addObject(rotI * Vec3D(1, 0, 0), rotI * 0.5 * bladeThickness * Vec3D(-1, 0, 0));
        AxisymmetricIntersectionOfWalls bladeAngle;
        bladeAngle.setSpecies(species);
        bladeAngle.setAxis({ci, -si, 0});
        bladeAngle.addObject({-1, 0, 0}, {bladeRadius, 0, 0});
        BasicIntersectionOfWalls blade;
        blade.setSpecies(species);
        blade.add(bladeSlice);
        blade.add(bladeAngle);
        blade.setAngularVelocity(angularVelocity);
        IntersectionOfWalls bladeRender(bladeSlice);
        for (int i = -25; i <= 25; i++)
        {
            const Mdouble a = bladeOpening * i / 50.0;
            const Mdouble co = cos(a*pi/180.0);
            const Mdouble so = sin(a*pi/180.0);
            Matrix3D rotO = {1, 0, 0, 0, co, -so, 0, so, co};
            const Vec3D n = rotI * (rotO * Vec3D(0, 0, -1));
            bladeRender.addObject(n, -bladeRadius * n);
        }
        //bladeRender.setAngularVelocity(angularVelocity);
        auto bladeLeft = wallHandler.copyAndAddObject(blade);
        bladeRender.setPrescribedPosition([position](Mdouble time){return position;});
        bladeRender.setPrescribedOrientation([bladeLeft](Mdouble time){return bladeLeft->getOrientation();});
        bladeLeft->addRenderedWall(bladeRender.copy());
        bladeLeft->setPosition(position);
        const Vec3D quarterTurn = {-2, 0, 0};
        for (unsigned i=0; i<=rotation; ++i) {
            bladeLeft->rotate(quarterTurn);
        }
    }
    
    
    ParticleSpecies* addSingleSpecies(SpeciesType type, Mdouble medianParticleRadius, Mdouble minParticleRadius, bool adjustTimeStep=false)
    {
        return ::addSingleSpecies(type, medianParticleRadius, minParticleRadius, *this, adjustTimeStep);
    }
    
    void printTime() const override 
    {
        std::cout << "time " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " timeMax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " particles " << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfRealObjectsLocal()
                  << std::endl;
    }

    Mdouble generateRadius() {
        //set up normal distribution
        static std::mt19937 gen(0);
        static std::normal_distribution<> d(particleRadius_, polydispersity_*particleRadius_);
        static const Mdouble radiusMin = (1.-1.5*polydispersity_)*particleRadius_;
        static const Mdouble radiusMax = (1.+1.5*polydispersity_)*particleRadius_;
        Mdouble chosenRad = d(gen);
        while (chosenRad > radiusMax || chosenRad < radiusMin){
            chosenRad = d(gen);
        }
        return chosenRad;
    }

    Mdouble generateRadius2() { /// ADDED by M; does this need to be changed as polydispersity looks weird in histogram plot
        static std::mt19937 gen(0);
        static std::normal_distribution<> d(particleRadius2_, polydispersity2_*particleRadius2_);
        static const Mdouble radiusMin = (1.-1.5*polydispersity2_)*particleRadius2_;
        static const Mdouble radiusMax = (1.+1.5*polydispersity2_)*particleRadius2_;
        Mdouble chosenRad2 = d(gen);
        while (chosenRad2 > radiusMax || chosenRad2 < radiusMin){
            chosenRad2 = d(gen);
        }
        return chosenRad2;
    }

    Mdouble getPeriodicMixerVolume() {
        const Mdouble casingVolume = constants::pi*mathsFunc::square(casingRadius)*2.0*bladeDistance;
        const Mdouble shaftVolume = constants::pi*mathsFunc::square(shaftRadius)*2.0*bladeDistance;
        const Mdouble bladeVolume = constants::pi*mathsFunc::square(bladeRadius)*bladeThickness*bladeOpening/360.;
        return casingVolume - shaftVolume -4.0*bladeVolume;
    }

    //remove particles according to (max.) outflow rate
    void setupInitialConditions() override
    {
        //add filler particles
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        const Mdouble mixerVolume = getPeriodicMixerVolume();
        p.setRadius(generateRadius());
        for (Mdouble fillVolume = volFracSpec1_*solidFraction_*mixerVolume; fillVolume>0; fillVolume-=p.getVolume()) {
            Mdouble r2Min = mathsFunc::square(shaftRadius+p.getRadius());
            Mdouble r2Max = mathsFunc::square(casingRadius-p.getRadius());
            Mdouble r2, y, z;
            do
            {
                y = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                z = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                r2 = y*y + z*z;
            } while (r2>r2Max || r2<r2Min);
            const Mdouble x = random.getRandomNumber(getXMin(),getXMax());
            p.setPosition({x, y, z});
            particleHandler.copyAndAddObject(p);
            p.setRadius(generateRadius());
        }
        logger(INFO,"Added % filler particles",particleHandler.getSize());
        int s1Added = particleHandler.getSize();

        p.setRadius(generateRadius2());
        for (Mdouble fillVolumeS2 = (1.0-volFracSpec1_)*solidFraction_*mixerVolume; fillVolumeS2>0; fillVolumeS2-=p.getVolume()) {
            Mdouble r2Min = mathsFunc::square(shaftRadius+p.getRadius());
            Mdouble r2Max = mathsFunc::square(casingRadius-p.getRadius());
            Mdouble r2, x, y;
            do
            {
                x = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                y = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                r2 = x*x + y*y;
            } while (r2>r2Max || r2<r2Min);
            const Mdouble z = random.getRandomNumber(getXMin(),getXMax());
            p.setPosition({x, y, z});
            particleHandler.copyAndAddObject(p);
            p.setRadius(generateRadius2());
        }
        logger(INFO,"Added % API particles",particleHandler.getSize()-s1Added);
        /// See above
    }

    void write(std::ostream& os, bool writeAllParticles) const override
    {
        Mercury3D::write(os, writeAllParticles);
        os << "fillerType " << fillerType;
        os << " apiType " << apiType;
        os << " particleRadius " << particleRadius_;
        os << " particleRadius2 " << particleRadius2_;
        os << " volFracSpec1 " << volFracSpec1_;
        os << " polydispersity " << polydispersity_;
        os << " polydispersity2 " << polydispersity2_;
        os << " solidFraction " << solidFraction_;
        os << " casingRadius " << casingRadius;
        os << " casingThickness " << casingThickness;
        os << " shaftRadius " << shaftRadius;
        os << " bladeRadius " << bladeRadius;
        os << " bladeDistance " << bladeDistance;
        os << " bladeInclination " << bladeInclination;
        os << " bladeOpening " << bladeOpening;
        os << " bladeThickness " << bladeThickness << '\n';
    }

    void read(std::istream& is)
    {
        Mercury3D::read(is);

        //read one line
        std::stringstream line;
        helpers::getLineFromStringStream(is, line);
        std::string dummy;
        line >> dummy >> fillerType;
        line >> dummy >> apiType;
        line >> dummy >> particleRadius_ >> dummy >> particleRadius2_ >> dummy >> volFracSpec1_ >> dummy >> polydispersity_ >> dummy >> polydispersity2_;
        line >> dummy >> solidFraction_ >> dummy >> casingRadius >> dummy >> casingThickness >> dummy >> shaftRadius >> dummy >> bladeRadius;
        line >> dummy >> bladeDistance >> dummy >> bladeInclination >> dummy >> bladeOpening >> dummy >> bladeThickness;
    }
    
    GCGPeriodic() {
        setGeometry();
    }
    
    GCGPeriodic(SpeciesType fillerType, SpeciesType apiType, Mdouble particleRadius, Mdouble particleRadius2, Mdouble volFracSpec1, Mdouble polydispersity, Mdouble polydispersity2, Mdouble bladeInclination, Mdouble shaftRadius_)
            : fillerType(fillerType), apiType(apiType), particleRadius_(particleRadius), particleRadius2_(particleRadius2), volFracSpec1_(volFracSpec1), polydispersity_(polydispersity), polydispersity2_(polydispersity2), solidFraction_(0), shaftRadius(shaftRadius_)
    {
        setGeometry();
        bladeInclination = bladeInclination;
    }
    
    void setGeometry (Mdouble scaleFactor = 1.0) {
        casingRadius = 35e-3*scaleFactor; //< inner radius of the casing
        casingThickness = 0.02*casingRadius; //< thickness of the casing (only needed for visualisation)
        shaftRadius *= scaleFactor; //< outer radius of the shaft
        bladeRadius = casingRadius - 1.5e-3*scaleFactor; //< outer radius of the blades
        bladeDistance = 30e-3*scaleFactor; //< distance between the blades
        bladeInclination = 20.0; //< angle by which the blade normal is inclined with respect to the x-axis
        bladeOpening = 60.0; //< arc angle of the slice
        bladeThickness = 2.5e-3*scaleFactor; //< width of the blade
    }

    void writeEneTimeStep(std::ostream& os) const override {
        const static long int width = os.precision() + 6;

        if ((eneFile.getCounter() == 1 || eneFile.getFileType() == FileType::MULTIPLE_FILES ||
             eneFile.getFileType() == FileType::MULTIPLE_FILES_PADDED) && !getAppend()) {
            os << std::setw(width+1)
               << "time " << std::setw(width+1)
               << "gravitEnergy " << std::setw(width+1) //gravitational potential energy
               << "tKineticEnergy " << std::setw(width+1) //translational kinetic energy
               << "rKineticEnergy " << std::setw(width+1) //rotational kE
               << "elasticEnergy " << std::setw(width+1)
               << "centerOfMassX " << std::setw(width+1)
               << "centerOfMassY " << std::setw(width+1)
               << "centerOfMassZ " << std::setw(width+1)
               << "mass " << std::setw(width+1)
               << "torqueX " << std::setw(width+1)
               << "torqueY " << std::setw(width+1)
               << "torqueZ\n";
        }

        const Mdouble m = particleHandler.getMass();
        const Vec3D com = m==0?Vec3D(0,0,0):particleHandler.getMassTimesPosition()/m;

        Vec3D bladeTorque = {0,0,0};
        for (auto w : wallHandler) {
            if (!w->getAngularVelocity().isZero()) {
                //bladeForce += w->getForce();
                bladeTorque += w->getTorque() + Vec3D::cross(w->getPosition(), w->getForce());
            }
        }

        //Ensure the numbers fit into a constant width column: for this we need the precision given by the operating system,
        //plus a few extra characters for characters like a minus and scientific notation.
        os << std::setw(width) << getTime()
        << " " << std::setw(width) << -Vec3D::dot(getGravity(), com)
        << " " << std::setw(width) << particleHandler.getKineticEnergy()
        << " " << std::setw(width) << particleHandler.getRotationalEnergy()
        << " " << std::setw(width) << getElasticEnergy()
        << " " << std::setw(width) << com.X
        << " " << std::setw(width) << com.Y
        << " " << std::setw(width) << com.Z
        << " " << std::setw(width) << m
        << " " << std::setw(width) << bladeTorque.X
        << " " << std::setw(width) << bladeTorque.Y
        << " " << std::setw(width) << bladeTorque.Z
        << std::endl;
    };

    void outputXBallsDataParticle(unsigned int i, unsigned int format, std::ostream& os) const override {
        const auto p = particleHandler.getObject(i);
        os << p->getPosition() << ' '
           << p->getVelocity() << ' '
           << p->getRadius() << ' '
           //<< p->getOrientation().getEuler() << ' '
           << "0 0 0 "
           << p->getAngularVelocity() << ' '
           << p->getSpecies()->getId() << ' '
           << p->getId() << '\n';
        //logger(ERROR,"S");
    };

    SpeciesType fillerType;
    SpeciesType apiType;
    Mdouble particleRadius_ = 0; //< median particle radius filler
    Mdouble particleRadius2_ = 0; //< median particle radius api
    Mdouble volFracSpec1_ = 1; //< volume fraction filler
    Mdouble polydispersity_ = 0; //< polydispersity filler (only used for certain species) 
    Mdouble polydispersity2_ = 0; //< polydispersity api (only used for certain species) 
    Mdouble solidFraction_ = 0; //<determines how many particles should be inserted into the mixer at the beginning
    Mdouble casingRadius; //< inner radius of the casing
    Mdouble casingThickness; //< thickness of the casing (only needed for visualisation)
    Mdouble shaftRadius; //< outer radius of the shaft
    Mdouble bladeRadius; //< outer radius of the blades
    Mdouble bladeDistance; //< distance between the blades
    Mdouble bladeInclination; //< angle by which the blade normal is inclined with respect to the x-axis
    Mdouble bladeOpening; //< arc angle of the slice
    Mdouble bladeThickness; //< width of the blade
};

#endif //MERCURY_GCG_H
