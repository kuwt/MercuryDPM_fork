#ifndef MERCURY_GCG_H
#define MERCURY_GCG_H

#include <Boundaries/DeletionBoundary.h>
#include "GCGPeriodic.h"
#include <iomanip>
#include <CG/TimeSmoothedCG.h>
#include <CG/TimeAveragedCG.h>
#include <Walls/TriangleWall.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <MercuryTime.h>
#include <Walls/BasicUnionOfWalls.h>
#include "BladePattern.h"

class GCG : public GCGPeriodic
{
public:

    GCG(SpeciesType fillerType, SpeciesType apiType, Mdouble fillRate, Mdouble rpm, Mdouble beltVelocity, Mdouble volFracSpec1_, BladePattern bladePattern, Mdouble scaleFactor, Mdouble bladeInclination, Mdouble shaftRadius=15e-3)
            : GCGPeriodic(fillerType, apiType, getMedianParticleRadius(getPSD(fillerType)),getMedianParticleRadius(getPSD(apiType)),volFracSpec1_, 0, 0, bladeInclination, shaftRadius), fillRate(fillRate), bladeNumber(getBladeNumber(bladePattern))
    {
        logger(INFO, "Fill rate % [m^3/s]", fillRate);
        //set name, gravity, dimensions
        setGeometry(scaleFactor);
        setName("GCG");
        setGravity({0, 0, -9.8});
        addSpecies();
        addWalls(rpm, beltVelocity);
        addBlades(bladePattern);
        addInsertionBoundaries();
    }

    
    GCG(SpeciesType fillerType, SpeciesType apiType, Mdouble fillRate, Mdouble rpm, Mdouble beltVelocity, Mdouble volFracSpec1_, bool longScrew, Mdouble scaleFactor, Mdouble bladeInclination, Mdouble shaftRadius=15e-3)
    : GCGPeriodic(fillerType, apiType, getMedianParticleRadius(getPSD(fillerType)),getMedianParticleRadius(getPSD(apiType)),volFracSpec1_, 0, 0, bladeInclination, shaftRadius), fillRate(fillRate), bladeNumber(longScrew?16:7)
    {
        logger(INFO, "Fill rate % [m^3/s]", fillRate);
        //set name, gravity, dimensions
        setGeometry(scaleFactor);
        setName("GCG");
        setGravity({0, 0, -9.8});
        addSpecies();
        addWalls(rpm, beltVelocity);
        addBlades(longScrew?(BladePattern::LongScrew):(BladePattern::ShortScrew));
        addInsertionBoundaries();
    }
    
    /**
     * Simulates the full GCG mixer; thus, the following changes are applied to GCGPeriodic:
     *  - the periodic boundary is removed
     *  - the casing is replaced by a casing with in-and outflow
     *  - particles are not immediately inserted, but over time, based on a (maximum) fillRate
     *  - if belt is specified, then a moving rim is inserted below the outflow
     *
     * All defining parameters are passed into the constructor.
     *
     * The dimensions of the system are:
     *  - |y| < casingRadius
     *  - outflowHeight < |z| < inflowHeight
     *  - |x| < casingLength/2,
     *
     * New geometric features:
     *  - casing: The cylinder with holes forming the casing
     *    parameters:
     *     - casingLength: length of the whole mixer
     *     - capLength: distance of in/outflow from the casing's cap
     *  - Inflow: the inflow hopper
     *    parameters:
     *     - inflowRadius
     *     - inflowHeight
     *     - fillRate
     *  - Outflow:
     *    Parameters
     *     - outflowHeight
     *
     *  @param particlesType speciesType describing the material/contact properties of the particle interaction
     *  @param wallType      speciesType describing the contact properties of the particle-wall interaction
     *  @param solidFraction     volume fraction of particles in the mixer: volume of particles/(volume inside casing - volume inside shaft and blades)
     *  @param rpm             rotation rate of the mixer (in revolutions per minute)
     *  @param particleRadius  mean particle radius (in m)
     *  @param particleRadius2 mean particle radius of species 2 (in m)
     *  @param volFracSpec1_ volume fraction of species 1
     *  @param polydispersity  relative standard deviation of the particle radius, assuming normal distribution for species 1
     *  @param polydispersity2  relative standard deviation of the particle radius, assuming normal distribution for species 2

     */
    GCG(SpeciesType particleType, Mdouble fillRate, Mdouble rpm, Mdouble beltVelocity, Mdouble particleRadius, Mdouble particleRadius2, Mdouble volFracSpec1_, Mdouble polydispersity, Mdouble polydispersity2, bool longScrew, Mdouble scaleFactor, Mdouble bladeInclination, Mdouble shaftRadius=15e-3)
            : GCGPeriodic(particleType, particleType, particleRadius, particleRadius2, volFracSpec1_, polydispersity, polydispersity2, bladeInclination, shaftRadius), fillRate(fillRate), bladeNumber(longScrew?16:7)
    {
        logger(INFO, "Mean particle radius1: % +/- % m", particleRadius, polydispersity * particleRadius);
        logger(INFO, "Mean particle radius2: % +/- % m", particleRadius2,
               polydispersity2 * particleRadius2); /// ADDED by M
        logger(INFO, "Fill rate % [m^3/s]", fillRate);
    
        //set name, gravity, dimensions
        setGeometry(scaleFactor);
        setName("GCG");
        setGravity({0,0,-9.8});
        addSpecies();
        addWalls(rpm, beltVelocity);
        addBlades(longScrew?(BladePattern::LongScrew):(BladePattern::ShortScrew));
        addInsertionBoundaries();
    }
    
    void addInsertionBoundaries() {
        CubeInsertionBoundary fillerInsertion;
        BaseParticle particle(speciesHandler.getObject(0));
        Vec3D pos = Vec3D(capLengthIn + hopperOffset, 0, casingRadius + exitHeight + hopperHeight);
        Vec3D min = pos + Vec3D(-0.9 * hopperRadius, -0.1 * hopperRadius, -0.2 * hopperHeight);
        Vec3D max = pos + Vec3D(-0.7 * hopperRadius, +0.1 * hopperRadius, 0);
        fillerInsertion.set(&particle, 100, min, max, Vec3D(0, 0, 0), Vec3D(0, 0, 0), particleRadius_*(1-polydispersity_), particleRadius_*(1+polydispersity_));
        fillerInsertion.setDistribution(CubeInsertionBoundary::Distribution::Normal_1_5);
        if (hasPSD(fillerType)) fillerInsertion.setPSD(getPSD(fillerType));
        fillerInsertion.setVolumeFlowRate(volFracSpec1_ * fillRate);
        fillerInsertionBoundary = boundaryHandler.copyAndAddObject(fillerInsertion);
    
        CubeInsertionBoundary apiInsertion;
        particle.setSpecies(speciesHandler.getObject(1));
        min = pos + Vec3D(0.5 * hopperRadius, -0.1 * hopperRadius, -0.2 * hopperHeight);
        max = pos + Vec3D(0.7 * hopperRadius, +0.1 * hopperRadius, 0);
        apiInsertion.set(&particle, 100, min, max, Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1, 1);
        apiInsertion.setDistribution(CubeInsertionBoundary::Distribution::Normal_1_5);
        if (hasPSD(apiType)) apiInsertion.setPSD(getPSD(apiType));
        apiInsertion.setVolumeFlowRate((1 - volFracSpec1_) * fillRate);
        apiInsertionBoundary = boundaryHandler.copyAndAddObject(apiInsertion);
    }

    void useVariableFillRates(std::vector<Mdouble> fillerFlowRate, Mdouble fillerSamplingInterval, std::vector<Mdouble> apiFlowRate, Mdouble apiSamplingInterval)
    {
        // add initial volume is necessary
        for (auto& v : fillerFlowRate) v += fillerInsertionBoundary->getInitialVolume();
        for (auto& v : apiFlowRate) v += apiInsertionBoundary->getInitialVolume();
        // set variable flow rates
        fillerInsertionBoundary->setVariableVolumeFlowRate(fillerFlowRate, fillerSamplingInterval);
        apiInsertionBoundary->setVariableVolumeFlowRate(apiFlowRate, apiSamplingInterval);
    }

    void addBlades (BladePattern pattern) {
        //set angular velocity and species
        Vec3D angularVelocity;
        for (auto w : wallHandler) {
            if (!w->getAngularVelocity().isZero()) {
                angularVelocity = w->getAngularVelocity();
                break;
            }
        }
        auto wSpecies = speciesHandler.getObject(speciesHandler.getSize()-2);

        const Mdouble center = 0.5*(inflowScrewLength+casingLength-outflowScrewLength);
        const Mdouble last = center + 0.5*bladeNumber*bladeDistance;
    
        std::string patternString  = getBladePatternString(pattern);
        Vec3D position = Vec3D(last, 0, 0);
        unsigned rotation = 0;
        unsigned count = 0;
        for (char patternChar : patternString) {
            if (patternChar!=' ')
            {
                Mdouble inclination;
                if (patternChar == 'f') inclination = bladeInclination;
                else if (patternChar == 'b') inclination = -bladeInclination;
                else if (patternChar == 'p') inclination = NaN;
                else logger(ERROR,"Blade pattern % does not exist", patternChar);
                logger(INFO, "inserting blade: type % inclination % position % rotation %", patternChar, inclination, position.X, rotation);
                addBlade(inclination, position, rotation, angularVelocity, wSpecies);
            }
            ++count;
            position.X -= 0.5*bladeDistance;
            if (count>=2*bladeNumber) {
                position.X = last;
                count = 0;
                ++rotation;
            }
        }
//            std::vector<double> sgn = {1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1};
//        for (int i = 0; i < bladeNumber; ++i)
//        {
//            addBladePair(sgn[i] * bladeInclination, {last - i * bladeDistance, 0, 0},
//                         {(i % 2 == 0) ? 2.0 : 0, 0, 0}, angularVelocity, wSpecies);
//        }
    }
    
    void addWalls(Mdouble rpm, Mdouble beltVelocity) {
        // angular velocity of axis, screw and blades
        const Vec3D angularVelocity = {2. * pi * rpm / 60., 0, 0};
    
        //set dimensions
        setDomain({0,-casingRadius,-beltHeight},{casingLength,casingRadius,casingRadius+exitHeight+hopperHeight});
        // set xballs arguments
        setXBallsAdditionalArguments("-v0 -solidf -cmode 8 -s 4.7 -oh 0 -o 175 -h 1000 -w 1500");

        //set the species of the walls
        logger.assert(speciesHandler.getNumberOfObjects()>=3,"Number of species too small");
        auto wSpecies = speciesHandler.getObject(speciesHandler.getSize()-2);
        auto bSpecies = speciesHandler.getLastObject();

        //save every 0.02s
        setSaveCount(0.02 / getTimeStep());

        // add mixer
        AxisymmetricIntersectionOfWalls mixer;
        mixer.setSpecies(wSpecies);
        mixer.setPosition({0,0,0});
        mixer.setAxis({1,0,0});
        mixer.addObject({1,0,0},{casingRadius,0,0});
        //mixer.addObject({0,0,1},{0.5*casingRadius,0,0});
        //wallHandler.copyAndAddObject(mixer);

        // add caps
        IntersectionOfWalls cap;
        cap.setSpecies(wSpecies);
        cap.addObject({1,0,0},{casingLength,0,0});
        cap.addObject({0,0,1},{0,0,-casingRadius});
        AxisymmetricIntersectionOfWalls capRender;
        capRender.setSpecies(wSpecies);
        capRender.setAxis({1,0,0});
        capRender.addObject({0,0,1},{0,0,casingLength});
        capRender.addObject({-1,0,0},{casingRadius,0,0});
        cap.addRenderedWall(capRender.copy());
        wallHandler.copyAndAddObject(cap);
        logger(INFO,"Wall %: cap",wallHandler.getSize()-1);

        // add caps
        InfiniteWall cap2;
        cap2.setSpecies(wSpecies);
        cap2.set({-1,0,0},{0,0,0});
        AxisymmetricIntersectionOfWalls capRender2;
        capRender2.setSpecies(wSpecies);
        capRender2.setAxis({1,0,0});
        capRender2.addObject({0,0,-1},{0,0,0});
        capRender2.addObject({-1,0,0},{casingRadius,0,0});
        cap2.addRenderedWall(capRender2.copy());
        wallHandler.copyAndAddObject(cap2);
        logger(INFO,"Wall %: cap2",wallHandler.getSize()-1);

        // add wall above hopper
        InfiniteWall hopperTop;
        hopperTop.setSpecies(wSpecies);
        hopperTop.set({0,0,1},{0,0,casingRadius+exitHeight+hopperHeight});
        wallHandler.copyAndAddObject(hopperTop);
        logger(INFO,"Wall %: hopperTop",wallHandler.getSize()-1);

        double n = 16;
        Mdouble t = std::tan(pi/n);
        Vec3D L = {capLengthIn,0,casingRadius+exitHeight}; // center low
        Vec3D H = {capLengthIn+hopperOffset,0,casingRadius+exitHeight+hopperHeight}; // center high
        Vec3D AL = {1,t,0}; // distance A to center low at 0-rotation
        Vec3D BL = {1,-t,0}; // distance B to center low at 0-rotation
        AL *= inflowRadius/AL.getLength();
        BL *= inflowRadius/BL.getLength();
        Vec3D CH = {1,t,0}; // distance C to center high at 0-rotation
        Vec3D DH = {1,-t,0}; // distance D to center high at 0-rotation
        CH *= hopperRadius/CH.getLength();
        DH *= hopperRadius/DH.getLength();

        TriangleWall hopper;
        hopper.setSpecies(wSpecies);
        for (double i=0; i<n; ++i)
        {
            if (2*i<n) hopper.setVTKVisibility(false);
            else hopper.setVTKVisibility(true);

            Mdouble angle = 2.0*pi*i/n, s = std::sin(angle), c = std::cos(angle);
            Matrix3D o = {c,s,0,-s,c,0,0,0,1};
            Vec3D A = L+o*AL;
            Vec3D B = L+o*BL;
            Vec3D C = H+o*CH;
            Vec3D D = H+o*DH;
            hopper.setVertices(A,B,C);
            wallHandler.copyAndAddObject(hopper);
            logger(INFO,"Wall %: hopper ABC %",wallHandler.getSize()-1,i);
            hopper.setVertices(B,C,D);
            wallHandler.copyAndAddObject(hopper);
            logger(INFO,"Wall %: hopper BCD %",wallHandler.getSize()-1,i);
        }
        
        // add inflow
        AxisymmetricIntersectionOfWalls inflow;
        inflow.setSpecies(wSpecies); //do we need species here?
        inflow.setPosition({capLengthIn,0,0});
        inflow.setAxis({0,0,1});
        // cylindrical wall
        inflow.addObject({1,0,0},{inflowRadius,0,0});
        // end of cylinder at z=0
        inflow.addObject({0,0,1},{0,0,0});
        // end of cylinder at z=0
        inflow.addObject({0,0,-1},{0,0,casingRadius+exitHeight});
        AxisymmetricIntersectionOfWalls inflowRender = inflow;
        //render only above z=0.8*mixerHeight
        inflowRender.addObject({0,0,1},{0,0,casingRadius});
        //render a wall thickness of 0.1*inflowRadius
        inflowRender.addObject({-1,0,0},{1.1*inflowRadius,0,0});

        // add outflow
        AxisymmetricIntersectionOfWalls outflow;
        outflow.setSpecies(wSpecies); //do we need species here?
        outflow.setPosition({casingLength-capLengthOut,0,0});
        outflow.setAxis({0,0,1});
        outflow.addObject({1,0,0},{inflowRadius,0,0});
        outflow.addObject({0,0,-1},{0,0,0});
        outflow.addObject({0,0,1},{0,0,-casingRadius-exitHeight});
        AxisymmetricIntersectionOfWalls outflowRender = outflow;
        outflowRender.addObject({0,0,-1},{0,0,-casingRadius});
        outflowRender.addObject({-1,0,0},{1.01*inflowRadius,0,0});

        //intersect mixer and inflow
        BasicIntersectionOfWalls mixerInflow;
        mixerInflow.setSpecies(wSpecies);
        mixerInflow.add(mixer);
        mixerInflow.add(inflow);
        //mixerInflow.addRenderedWall(mixer.copy());
        //mixerInflow.addRenderedWall(inflowRender.copy());
        //wallHandler.copyAndAddObject(mixerInflow);
        //logger(INFO,"Wall %: mixerInflow",wallHandler.getSize()-1);

        //intersect mixer and outflow
        BasicIntersectionOfWalls mixerOutflow;
        mixerOutflow.setSpecies(wSpecies);
        mixerOutflow.add(mixer);
        mixerOutflow.add(outflow);
        //mixerOutflow.addRenderedWall(outflowRender.copy());
        //wallHandler.copyAndAddObject(mixerOutflow);
        //logger(INFO,"Wall %: mixerOutflow",wallHandler.getSize()-1);

        BasicUnionOfWalls unionMixer;
        unionMixer.setSpecies(wSpecies);
        unionMixer.add(mixerOutflow);
        unionMixer.add(mixerInflow);
        unionMixer.addRenderedWall(inflowRender.copy());
        wallHandler.copyAndAddObject(unionMixer);
        logger(INFO,"Wall %: unionMixer",wallHandler.getSize()-1);

        // add shaft
        AxisymmetricIntersectionOfWalls shaft;
        shaft.setSpecies(wSpecies);
        shaft.setPosition({0, 0, 0});
        shaft.setAxis({1, 0, 0});
        shaft.addObject({-1, 0, 0}, {shaftRadius, 0, 0});
        shaft.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(shaft);
        logger(INFO,"Wall %: shaft",wallHandler.getSize()-1);

        //Screw(Vec3D start, Mdouble l, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness);
        Screw inflowScrew({0,0,0}, inflowScrewLength, bladeRadius, -inflowScrewLength/screwPitch, rpm/60.0, bladeThickness/2.0, ScrewType::singleHelix);
        inflowScrew.setSpecies(bSpecies); //make the helix rubber-like to reduce backsplatter
        inflowScrew.setPosition({0,0,0});
        inflowScrew.setOrientationViaNormal({1,0,0});
        wallHandler.copyAndAddObject(inflowScrew);
        logger(INFO,"Wall %: inflowScrew",wallHandler.getSize()-1);

        Screw outflowScrew({0,0,0}, outflowScrewLength, bladeRadius, outflowScrewLength/screwPitch, -rpm/60.0, bladeThickness/2.0, ScrewType::singleHelix);
        outflowScrew.setSpecies(wSpecies);
        outflowScrew.setPosition({casingLength,0,0});
        outflowScrew.setOrientationViaNormal({-1,0,0});
        const Vec3D quarterTurn = {2, 0, 0};
        outflowScrew.rotate(quarterTurn);
        wallHandler.copyAndAddObject(outflowScrew);
        logger(INFO,"Wall %: outflowScrew",wallHandler.getSize()-1);

        // add belt
        if (beltVelocity!=0)
        {
            InfiniteWall belt;
            belt.setSpecies(bSpecies);
            belt.setVelocity({beltVelocity,0,0});
            belt.set({0, 0, -1},{0,0,-beltHeight});
            wallHandler.copyAndAddObject(belt);
            logger(INFO,"Wall %: belt",wallHandler.getSize()-1);

            DeletionBoundary deletionBoundary;
            deletionBoundary.set({1,0,0},1.2*casingLength);
            boundaryHandler.copyAndAddObject(deletionBoundary);
            logger(INFO,"Boundary %: deletionBoundary",boundaryHandler.getSize()-1);
        } else {
            DeletionBoundary deletionBoundary;
            deletionBoundary.set({0,0,-1},beltHeight);
            boundaryHandler.copyAndAddObject(deletionBoundary);
            logger(INFO,"Boundary %: deletionBoundary",boundaryHandler.getSize()-1);
        }

        logger(INFO,"Added % walls", wallHandler.getNumberOfObjects());
    }
    
    /**
     * Constructor used for restarting simulations
     * @param restart name of the restart file
     * @param timeMax new final time for the restarted simulation
     */
    GCG(std::string restart, Mdouble timeMax, std::string newName = "")
    : GCGPeriodic()
    {
        /// Can this be changed to automatically fill in zero for everything?
        /// Changed to 5 zeros instead of 2 zeros
        // restart from .restart file
        readRestartFile(restart);
        if (timeMax <= getTime()) {
            logger(ERROR,"restarted simulation % is already at time %, higher than the requested timeMax value of %",getName(),getTime(), getTimeMax());
        }
        // change timeMax and, if necessary, name
        setTimeMax(timeMax);
        logger(INFO,"restarted simulation %, set timeMax to %",getName(),getTimeMax());
        if (!newName.empty()) {
            setName(newName);
        }
        //the visualisation of the blades requires special treatment, as lambda functions cannot be restarted
        for (auto wall : wallHandler)
        {
            auto blade = dynamic_cast<BasicIntersectionOfWalls*>(wall);
            if (blade!=nullptr && !blade->getPosition().isZero())
            {
                Vec3D position = wall->getPosition();
                wall->getRenderedWall(0)->setPrescribedPosition([position](Mdouble time) { return position; });
                wall->getRenderedWall(0)->setPrescribedOrientation([blade](Mdouble time) { return blade->getOrientation(); });
            }
        }
        
        fillerInsertionBoundary = dynamic_cast<CubeInsertionBoundary*>(boundaryHandler.getObject(boundaryHandler.getSize()-2));
        apiInsertionBoundary = dynamic_cast<CubeInsertionBoundary*>(boundaryHandler.getLastObject());
        
        logger.assert(fillerInsertionBoundary!= nullptr,"fillerInsertionBoundary not well-defined");
        logger.assert(apiInsertionBoundary!= nullptr,"apiInsertionBoundary not well-defined");
        // after restart, solve
        solve();
    }

    void initialFillBasedOnMeanResidenceTime(Mdouble meanResidenceTime) {
        // meanResidenceTime [s]
        // fillRate [m^3/s]
        initialFillVolume = fillRate*meanResidenceTime;
        const Mdouble mixerVolume = getPeriodicMixerVolume()*bladeNumber;
        const Mdouble solidFraction = initialFillVolume/mixerVolume;
        logger(INFO,"Blender will be initially filled with % ml of particles (solid fraction of %) to reach steady state quickly",initialFillVolume*1e6, solidFraction);
    }
    
    /**
     * Allows for some preliminary statistical output without resorting to output files:
     *  - XZT statistics, smoothed in time and space with w=2*d (Heaviside), wt=0.2 (Gaussian)
     *  - XT statistics, smoothed in time and space with w=2*d (Heaviside), wt=0.2 (Gaussian)
     *  - global T statistics
     */
    void addCG () {
        addCG3D();
        addCG2D();
        addCG1D();
        addCG0D();
        helpers::writeToFile(
                getName()+".m",
                "addpath('../../../MercuryCG')\n"
                "clear variables\n"
                "O=readMercuryCG('GCGSelfTest.O.stat');\n"
                "X=readMercuryCG('GCGSelfTest.X.stat');\n"
                "XZ=readMercuryCG('GCGSelfTest.XZ.stat');\n"
                "XYZ=readMercuryCG('GCGSelfTest.XYZ.stat');\n"
                "O1=readMercuryCG('GCGSelfTest.O.S1.stat');\n"
                "X1=readMercuryCG('GCGSelfTest.X.S1.stat');\n"
                "XYZ1=readMercuryCG('GCGSelfTest.XYZ.S1.stat');\n"
                "\n"
                "figure(1)\n"
                "plot(O.t,O.VelocityX)\n"
                "xlabel('t'), ylabel('v_x')\n"
                "\n"
                "figure(2)\n"
                "plot(X.x,X.VelocityX)\n"
                "legend(num2str(X.t(1,:)','t=%3.2g'))\n"
                "xlabel('x'), ylabel('v_x')\n"
                "\n"
                "figure(3)\n"
                "h = pcolor(XZ.x(:,:,end),XZ.z(:,:,end),XZ.VelocityX(:,:,end));\n"
                "set(h,'LineStyle','none')\n"
                "xlabel('x'), ylabel('z'), title('v_x'), axis equal tight; colorbar\n"
                "\n"
                "figure(4)\n"
                "%ArntzOtterBrielsBussmannBeeftinkBoom2008.pdf\n"
                "x = XYZ.ParticleNumber;\n"
                "x1 = XYZ1.ParticleNumber;\n"
                "x2 = x-x1;\n"
                "lx = log(x); lx(x==0)=0;\n"
                "lx1 = log(x1); lx1(x1==0)=0;\n"
                "lx2 = log(x2); lx2(x2==0)=0;\n"
                "S = mean(mean(mean((x1.*lx1+x2.*lx2).*x,1),2),3);\n"
                "SS = mean(mean(mean(x.*lx.*x,1),2),3);\n"
                "SM = mean(mean(mean(x.*(lx-log(2)).*x,1),2),3);\n"
                "phi = (S-SS)./(SM-SS);\n"
                "plot(squeeze(XZ.t(1,1,:)),squeeze(phi));\n"
                "xlabel('t'), ylabel('S'), title('Mixing entropy')\n"
                "\n"
                "figure(5)\n"
                "ene = importdata('GCGSelfTest.ene')\n"
                "ene.t = ene.data(:,1);\n"
                "ene.Torque = ene.data(:,12:14);\n"
                "plot(ene.t,ene.Torque)\n"
                "legend('X','Y','Z')\n"
                "xlabel('t'), ylabel('\\tau'), title('Torque')\n");
        logger(INFO,"Run "+getName()+".m to view basic cg output");
        
        helpers::writeToFile(getName()+"_Torque.gnu",
                             "set xlabel 'time'\n"
                             "set ylabel 'torque'\n"
                             "p 'GCGSelfTest.ene' u 1:12\n");
        logger(INFO,"Run 'gnuplot "+getName()+"_Torque.gnu' to view torque over time.");
        
    }

private:
    
    void printTime() const override
    {
        static Time time;
        logger(INFO,"time %\tparticles %\tfill [ml/s] % %\t[#] % %\tCPU:SIM-time %",
               std::round(getTime()*100)/100,
               particleHandler.getNumberOfRealObjectsLocal(),
               fillerInsertionBoundary->getVolumeOfParticlesInserted()/getTime()*1e6,
               apiInsertionBoundary->getVolumeOfParticlesInserted()/getTime()*1e6,
               fillerInsertionBoundary->getNumberOfParticlesInserted(),
               apiInsertionBoundary->getNumberOfParticlesInserted(),
               time.toctic());
    }

    void setupInitialConditions() override
    {
        if (apiInsertionBoundary== nullptr)
        {
            //add particles
            BaseParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            const Mdouble mixerVolume = getPeriodicMixerVolume();
    
            //p.setRadius(generateRadius());
            for (Mdouble fillVolume = volFracSpec1_ * initialFillVolume; fillVolume > 0; fillVolume -= p.getVolume())
            {
                p.setRadius(generateRadius());
                Mdouble r2Min = mathsFunc::square(shaftRadius + p.getRadius());
                Mdouble r2Max = mathsFunc::square(casingRadius - p.getRadius());
                Mdouble r2, y, z;
                do
                {
                    y = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                    z = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                    r2 = y * y + z * z;
                } while (r2 > r2Max || r2 < r2Min);
                const Mdouble x = random.getRandomNumber(capLengthIn, casingLength - capLengthOut);
                p.setPosition({x, y, z});
                particleHandler.copyAndAddObject(p);
            }
            logger(INFO, "Added % particles species 1", particleHandler.getSize());
    
            //Uncomment to add a single particle close to the deletion boundary to check its behavior
            //p.setRadius(1.5*particleRadius_);
            //p.setPosition({1.19*casingLength, 0, p.getRadius()-beltHeight});
            //particleHandler.copyAndAddObject(p);
    
            int s1Added = particleHandler.getSize();
    
            for (Mdouble fillVolumeS2 = (1.0 - volFracSpec1_) * initialFillVolume;
                 fillVolumeS2 > 0; fillVolumeS2 -= p.getVolume())
            {
                p.setRadius(generateRadius2());
                Mdouble r2Min = mathsFunc::square(shaftRadius + p.getRadius());
                Mdouble r2Max = mathsFunc::square(casingRadius - p.getRadius());
                Mdouble r2, y, z;
                do
                {
                    y = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                    z = (casingRadius - p.getRadius()) * random.getRandomNumber(-1, 1);
                    r2 = y * y + z * z;
                } while (r2 > r2Max || r2 < r2Min);
                const Mdouble x = random.getRandomNumber(capLengthIn, casingLength - capLengthOut);
                p.setPosition({x, y, z});
                particleHandler.copyAndAddObject(p);
            }
            logger(INFO, "Added % particles species 2", particleHandler.getSize() - s1Added);
        } else {
            CubeInsertionBoundary fillerInsertion;
            //fillerInsertion.setHandler(boundaryHandler);
            BaseParticle particle(speciesHandler.getObject(0));
            Vec3D min = Vec3D(capLengthIn,-casingRadius,-casingRadius);
            Vec3D max = Vec3D(casingLength - capLengthOut,casingRadius,casingRadius);
            fillerInsertion.set(&particle, 100, min, max, Vec3D(0, 0, 0), Vec3D(0, 0, 0), particleRadius_*(1-polydispersity_), particleRadius_*(1+polydispersity_));
            fillerInsertion.setDistribution(CubeInsertionBoundary::Distribution::Normal_1_5);
            if (hasPSD(fillerType)) fillerInsertion.setPSD(getPSD(fillerType));
            fillerInsertion.setInitialVolume(volFracSpec1_ * initialFillVolume);
            fillerInsertion.checkBoundaryBeforeTimeStep(this);
            logger(INFO, "Added % particles species 1", particleHandler.getSize());
    
            CubeInsertionBoundary apiInsertion;
            particle.setSpecies(speciesHandler.getObject(1));
            apiInsertion.set(&particle, 100, min, max, Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1, 1);
            apiInsertion.setDistribution(CubeInsertionBoundary::Distribution::Normal_1_5);
            if (hasPSD(apiType)) apiInsertion.setPSD(getPSD(apiType));
            apiInsertion.setInitialVolume((1 - volFracSpec1_) * initialFillVolume);
            apiInsertion.checkBoundaryBeforeTimeStep(this);
            logger(INFO, "Added % particles species 2", particleHandler.getSize());
        }
        
        //add processor number to file name if necessary
        std::string q = (NUMBER_OF_PROCESSORS > 1)?std::to_string(PROCESSOR_ID):"";
        inflowTracker.open(getName()+".in"+q);
        inflowTracker << "Time\tID\tRadius\tSpecies" << std::endl;
        outflowTracker.open(getName()+".out"+q);
        outflowTracker << "Time\tID\tRadius" << std::endl;
    }

    void actionsOnRestart() override {
        inflowTracker.open(getName()+".in",std::fstream::out | std::fstream::app);
        outflowTracker.open(getName()+".out",std::fstream::out | std::fstream::app);
    }
    
    void actionsAfterSolve () override {
        inflowTracker.close();
        outflowTracker.close();
    }
    
    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        // Functions after this statement only get executed every 100th timestep (if counter==100)
        static unsigned counter = 99;
        if (++counter!=100) return;
        else counter = 0;
    
        inflowTracker << getTime()
            << '\t' << fillerInsertionBoundary->getMassOfParticlesInserted()
            << '\t' << apiInsertionBoundary->getMassOfParticlesInserted()
            << std::endl;
    }

    void checkInteractionWithBoundaries () override {
        for (const auto p: particleHandler) {
            if (p->getPosition().X>1.2*getXMax()) {
                outflowTracker << getTime() << '\t' << p->getId() << '\t' << p->getRadius() << '\n';
            }
        }
        outflowTracker.flush();
        DPMBase::checkInteractionWithBoundaries();
    }
    
    void write(std::ostream& os, bool writeAllParticles) const override
    {
        GCGPeriodic::write(os, writeAllParticles);
        os << "fillRate " << fillRate
        << " bladeNumber " << bladeNumber << '\n';
    }

    void read(std::istream& is)
    {
        GCGPeriodic::read(is);

        std::stringstream line;
        helpers::getLineFromStringStream(is, line);
        std::string dummy;

        line >> dummy >> fillRate >> dummy >> bladeNumber;
        logger(INFO,"% % %",SpeciesType::MPT,fillerType);
        
        logger(INFO,"restarted with fillRate %, bladeNumber % % % %", fillRate, bladeNumber, fillerType, apiType, boundaryHandler.getSize());
    }

    void addCG3D () {
        TimeSmoothedCG<CGCoordinates::XYZ,CGFunctions::Heaviside> cg;
        cg.setXGrid(0, casingLength, 0.01);
        cg.setYGrid(-casingRadius, casingRadius, 0.007);
        cg.setZGrid(-casingRadius, casingRadius, 0.007);
        cg.setWidth(0.005);
        cg.statFile.setSaveCount(200);
        cg.statFile.setName(getName()+".XYZ.stat");
        cg.setTimeStep(0.2);
        cg.setWidthTime(0.2);
        cgHandler.copyAndAddObject(cg);
        
        //if we have mixed inflow, add species-0 statistics to compute the mixing index
        if (volFracSpec1_==1.0) return;
        cg.selectSpecies(0);
        cg.statFile.setName(getName()+".XYZ.S1.stat");
        cgHandler.copyAndAddObject(cg);
    }

    void addCG2D () {
        TimeSmoothedCG<CGCoordinates::XZ,CGFunctions::Heaviside> cg;
        cg.setXGrid(0, casingLength, 0.005);
        cg.setZGrid(-casingRadius, casingRadius, 0.005);
        cg.setWidth(2.0*particleRadius_);
        cg.statFile.setSaveCount(50);
        cg.statFile.setName(getName()+".XZ.stat");
        cg.setTimeStep(0.2);
        cg.setWidthTime(0.2);
        cgHandler.copyAndAddObject(cg);
    
        //if we have mixed inflow, add species-0 statistics to compute the mixing index
        if (volFracSpec1_==1.0) return;
        cg.selectSpecies(0);
        cg.statFile.setName(getName()+".XZ.S1.stat");
        cgHandler.copyAndAddObject(cg);
    }
    
    void addCG1D () {
        TimeSmoothedCG<CGCoordinates::X,CGFunctions::Heaviside> cg;
        cg.setXGrid(0, casingLength, 0.005);
        cg.setWidth(2.0*particleRadius_);
        cg.statFile.setSaveCount(50);
        cg.statFile.setName(getName()+".X.stat");
        cg.setTimeStep(0.2);
        cg.setWidthTime(0.2);
        cgHandler.copyAndAddObject(cg);
    
        //if we have mixed inflow, add species-0 statistics to compute the mixing index
        if (volFracSpec1_==1.0) return;
        cg.selectSpecies(0);
        cg.statFile.setName(getName()+".X.S1.stat");
        cgHandler.copyAndAddObject(cg);
    }
    
    void addCG0D () {
        CG<CGCoordinates::O> cg;
        cg.statFile.setSaveCount(50);
        cg.statFile.setName(getName()+".O.stat");
        cgHandler.copyAndAddObject(cg);
    
        //if we have mixed inflow, add species-0 statistics to compute the mixing index
        if (volFracSpec1_==1.0) return;
        cg.selectSpecies(0);
        cg.statFile.setName(getName()+".O.S1.stat");
        cgHandler.copyAndAddObject(cg);
    }
    
    void setGeometry (Mdouble scaleFactor = 1.0) {
        GCGPeriodic::setGeometry(scaleFactor);
        capLengthOut = 37.5e-3*scaleFactor; //length of caps, i.e. mixer length beyond in and outflow
        capLengthIn = 42.5e-3*scaleFactor; //length of caps, i.e. mixer length beyond in and outflow
        inflowRadius = 0.875*25.4e-3*scaleFactor; //radius of the in- and outflow (diam=1.75'')
        exitHeight = 45e-3*scaleFactor; //height of the in and outflow cylinders
        hopperHeight = 275e-3;
        hopperRadius = 96e-3;
        hopperOffset = 51e-3;
        //determines where the particles flow in
        // (the value in the sqrt is the factor by which the density is increased in the hopper)
        beltHeight = casingRadius+1.5*exitHeight; //radial distance of outflow belt from mixer axis
        screwPitch = 20e-3*scaleFactor; //offset between the screw blades at the in- and outflow
        inflowScrewLength = 2.75*screwPitch; //length of the inflow screw
        outflowScrewLength = 2*screwPitch; //length of the outflow screw
        casingLength = 10e-3*scaleFactor + inflowScrewLength + outflowScrewLength + bladeNumber * bladeDistance;
    }
    
    Mdouble fillRate;
    Mdouble bladeNumber; //< number of blade pairs
    Mdouble initialFillVolume = 0; //determines if particles should be inserted into the mixer at the beginning
    //Geometric parameters
    Mdouble capLengthOut; //length of caps, i.e. mixer length beyond in and outflow
    Mdouble capLengthIn; //length of caps, i.e. mixer length beyond in and outflow
    Mdouble inflowRadius; //radius of the in- and outflow (diam=1.75'')
    Mdouble exitHeight; //height of the in and outflow cylinders
    Mdouble hopperHeight;
    Mdouble hopperRadius;
    Mdouble hopperOffset;
    //determines where the particles flow in
    // (the value in the sqrt is the factor by which the density is increased in the hopper)
    Mdouble beltHeight; //radial distance of outflow belt from mixer axis
    Mdouble screwPitch; //offset between the screw blades at the in- and outflow
    Mdouble inflowScrewLength; //length of the inflow screw
    Mdouble outflowScrewLength; //length of the outflow screw
    Mdouble casingLength; //total length of the casing
    //Output files
    std::ofstream inflowTracker;
    std::ofstream outflowTracker;
    //Pointers to the inflow boundaries needed in actionBeforeTimeStep
    CubeInsertionBoundary* fillerInsertionBoundary = nullptr;
    CubeInsertionBoundary* apiInsertionBoundary = nullptr;
};

#endif //MERCURY_GCG_H
