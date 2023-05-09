#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/WearableTriangleMeshWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include <cstring>
#include <chrono>

class Screen : public Mercury3D
{
public:
    
    void setupInitialConditions() override
    {
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(1);
        species->setCollisionTimeAndRestitutionCoefficient(20.0*getTimeStep(), 0.1, species->getMassFromRadius(0.5 * (minParticleRadius_ + maxParticleRadius_)));
        
        species->setRollingFrictionCoefficient(0.4);
        species->setRollingStiffness(species->getStiffness() * 2.0 / 7.0);
        species->setRollingDissipation(species->getDissipation() * 2.0 / 7.0);
        
        species->setSlidingFrictionCoefficient(0.25);
        species->setSlidingStiffness(species->getStiffness() * 2.0 / 7.0);
        species->setSlidingDissipation(species->getDissipation() * 2.0 / 7.0);
        
        species->setTorsionFrictionCoefficient(0.1);
        species->setTorsionStiffness(species->getStiffness() * 2.0 / 7.0);
        species->setTorsionDissipation(species->getDissipation() * 2.0 / 7.0);
        
        // --- Screen with grizzly bars --------------------------------------------------------------------------------
        for (int i = 0; i < numberOfGrizzlyBars_; i++)
        {
            auto w = createGrizzlyBar(grizzlyBarHeight_, grizzlyBarWidthBack_, grizzlyBarWidthFront_, grizzlyBarLength_,
                                      numSegmentsX_, numSegmentsY_, numSegmentsZ_);
            w.setSpecies(species);
            w.setWearCoefficient(wearCoefficient_);
            w.setHardness(hardness_);
            w.setWearAcceleration(wearAcceleration_);
            
            const Vec3D restPositionGrizzlyBar = screenRestPosition_ + Vec3D(0.0, (i + 0.5) * grizzlyBarSpacing_, 0.0);
            w.setOrientationViaEuler(Vec3D(0.0, screenAngle_, 0.0));
            w.setPosition(restPositionGrizzlyBar);
            w.setPrescribedPosition(getVibratingMotion(vibrationAmplitude_, vibrationFrequency_, screenAngle_, restPositionGrizzlyBar));
            
            wallHandler.copyAndAddObject(w);
        }
    
        // --- Periodic boundaries -------------------------------------------------------------------------------------
        // Set YMax so that domain fits tightly around grizzly bars. Minimum of 1 spacing distance.
        setYMax((numberOfGrizzlyBars_ <= 1 ? 1 : numberOfGrizzlyBars_) * grizzlyBarSpacing_);
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(getXMin(), getYMin(), getZMin()), Vec3D(getXMin(), getYMax(), getZMin()));
        boundaryHandler.copyAndAddObject(b0);
    
        // --- Input screen --------------------------------------------------------------------------------------------
        const Mdouble screenSpacing = screenSpacing_ + 2.0 * grizzlyBarHeight_;
        const Vec3D restPositionInputScreen = screenRestPosition_ + Vec3D(screenSpacing * std::sin(screenAngle_), 0.0, screenSpacing * std::cos(screenAngle_));

        IntersectionOfWalls w1;
        w1.setSpecies(species);
        w1.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, 0.0));
        w1.addObject(Vec3D(-1.0, 0.0, 0.0), Vec3D(0.0, 0.0, 0.0));
        w1.setPosition(restPositionInputScreen);
        w1.setOrientationViaEuler(Vec3D(0.0, screenAngle_, 0.0));
        w1.setPrescribedPosition(getVibratingMotion(vibrationAmplitude_, vibrationFrequency_, screenAngle_, restPositionInputScreen));
        wallHandler.copyAndAddObject(w1);
    
        // --- Particle insertion --------------------------------------------------------------------------------------
        SphericalParticle p0;
        p0.setSpecies(species);
        
        CubeInsertionBoundary ib;
        ib.set(p0, 1, Vec3D(-grizzlyBarLength_, getYMin() + maxParticleRadius_, restPositionInputScreen.Z + grizzlyBarLength_ * std::sin(screenAngle_)),
               Vec3D(-0.5 * grizzlyBarLength_, getYMax() - maxParticleRadius_, getZMax()), Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 0.0, 0.0),
               minParticleRadius_, maxParticleRadius_);
        boundaryHandler.copyAndAddObject(ib);
        
        DeletionBoundary db;
        db.set(Vec3D(0.0, 0.0, -1.0), 0.0);
        boundaryHandler.copyAndAddObject(db);
    }
    
    void actionsOnRestart() override
    {
        for (int i = 0; i < wallHandler.getNumberOfObjects() - 1; i++)
        {
            const Vec3D restPositionGrizzlyBar = screenRestPosition_ + Vec3D(0.0, (i + 0.5) * grizzlyBarSpacing_, 0.0);
            auto w = wallHandler.getObject(i);
            w->setPrescribedPosition(getVibratingMotion(vibrationAmplitude_, vibrationFrequency_, screenAngle_, restPositionGrizzlyBar));
        }
    
        const Mdouble screenSpacing = screenSpacing_ + 2.0 * grizzlyBarHeight_;
        const Vec3D restPositionInputScreen = screenRestPosition_ + Vec3D(screenSpacing * std::sin(screenAngle_), 0.0, screenSpacing * std::cos(screenAngle_));
        auto w = wallHandler.getLastObject();
        w->setPrescribedPosition(getVibratingMotion(vibrationAmplitude_, vibrationFrequency_, screenAngle_, restPositionInputScreen));
    }
    
    void actionsAfterTimeStep() override
    {
    }
    
private:
    std::function<Vec3D(Mdouble)> getVibratingMotion(Mdouble amplitude, Mdouble frequency, Mdouble angle, Vec3D restPosition)
    {
        // Perpendicular sin
//        return [amplitude, frequency, angle, restPosition] (Mdouble time) {
//            Mdouble val = amplitude * std::sin(2.0 * constants::pi * frequency * time);
//            Vec3D dir(std::sin(angle), 0.0, std::cos(angle));
//            return restPosition + val * dir;
//        };

        // Circular
//        return [amplitude, frequency, angle, restPosition] (Mdouble time) {
//            Mdouble t = 2.0 * constants::pi * frequency * time;
//            return restPosition + amplitude * Vec3D(std::cos(-t), 0.0, std::sin(-t));
//        };

        // Horizontal jerking motion
        return [amplitude, frequency, angle, restPosition] (Mdouble time) {
            Mdouble t = 2.0 * constants::pi * frequency * time;
            Mdouble lengthStage2 = 0.2 * constants::pi;
            Mdouble tmod = std::fmod(t, constants::pi + lengthStage2);
            Mdouble val = tmod < constants::pi ? -std::cos(tmod) + 1.0 : -2.0 / lengthStage2 * (tmod - constants::pi) + 2.0;
            return restPosition + amplitude * Vec3D(val, 0.0, 0.0);
        };
    }
    
    WearableTriangleMeshWall createGrizzlyBar(Mdouble height, Mdouble widthBack, Mdouble widthFront, Mdouble length,
                                              int numSegmentsX, int numSegmentsY, int numSegmentsZ)
    {
        WearableTriangleMeshWall w, temp;
        
        // Top
        w.createFourPointMesh(Vec3D(0.0, -widthBack / 2, height), Vec3D(0.0, widthBack / 2, height),
                              Vec3D(length, widthFront / 2, height), Vec3D(length, -widthFront / 2, height),
                              numSegmentsX, numSegmentsY);
    
        // Bottom
        temp.createFourPointMesh(Vec3D(0.0, -widthBack / 2, 0.0), Vec3D(0.0, widthBack / 2, 0.0),
                                 Vec3D(length, widthFront / 2, 0.0), Vec3D(length, -widthFront / 2, 0.0),
                                 numSegmentsX, numSegmentsY);
        w.addToMesh(temp);
    
        // Left
        temp.createFourPointMesh(Vec3D(0.0, -widthBack / 2, 0.0), Vec3D(0.0, -widthBack / 2, height),
                                 Vec3D(length, -widthFront / 2, height), Vec3D(length, -widthFront / 2, 0.0),
                                 numSegmentsX, numSegmentsZ);
        w.addToMesh(temp);
    
        // Right
        temp.createFourPointMesh(Vec3D(0.0, widthBack / 2, 0.0), Vec3D(0.0, widthBack / 2, height),
                                 Vec3D(length, widthFront / 2, height), Vec3D(length, widthFront / 2, 0.0),
                                 numSegmentsX, numSegmentsZ);
        w.addToMesh(temp);
    
        // Back
        temp.createFourPointMesh(Vec3D(0.0, -widthBack / 2, 0.0), Vec3D(0.0, -widthBack / 2, height),
                                 Vec3D(0.0, widthBack / 2, height), Vec3D(0.0, widthBack / 2, 0.0),
                                 numSegmentsY, numSegmentsZ);
        w.addToMesh(temp);
    
        // Front
        temp.createFourPointMesh(Vec3D(length, -widthFront / 2, 0.0), Vec3D(length, -widthFront / 2, height),
                                 Vec3D(length, widthFront / 2, height), Vec3D(length, widthFront / 2, 0.0),
                                 numSegmentsY, numSegmentsZ);
        w.addToMesh(temp);
        
        return w;
    }

public:
    
    bool readNextArgument(int& i, int argc, char* argv[]) override
    {
        // Call base function and if it has read the next argument, return true.
        if (Mercury3D::readNextArgument(i, argc, argv))
            return true;
        
        // When base function has not found the argument, it might be a custom argument defined here.
        
        if (!strcmp(argv[i], "-wearCoefficient"))
        {
            wearCoefficient_ = std::stod(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-hardness"))
        {
            hardness_ = std::stod(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-wearAcceleration"))
        {
            wearAcceleration_ = std::stod(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-omp"))
        {
            setNumberOfOMPThreads(std::stoi(argv[i + 1]));
        }
        else if (!strcmp(argv[i], "-numSegmentsX"))
        {
            numSegmentsX_ = std::stoi(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-numSegmentsY"))
        {
            numSegmentsY_ = std::stoi(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-numSegmentsZ"))
        {
            numSegmentsZ_ = std::stoi(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-particleRadiusRange"))
        {
            Mdouble min = std::stod(argv[i + 1]);
            Mdouble max = std::stod(argv[i + 2]);
            i++;
            setParticleRadiusRange(min, max);
        }
        else
        {
            // Argument passed is not a custom flag
            return false;
        }
        
        // Argument passed has been found as a custom flag
        return true;
    }
    
    void setScreenProperties(Vec3D restPosition, Mdouble angle, Mdouble spacing, int numberOfGrizzlyBars)
    {
        screenRestPosition_ = restPosition;
        screenAngle_ = angle;
        screenSpacing_ = spacing;
        numberOfGrizzlyBars_ = numberOfGrizzlyBars;
    }
    
    void setVibrationProperties(Mdouble amplitude, Mdouble frequency)
    {
        vibrationAmplitude_ = amplitude;
        vibrationFrequency_ = frequency;
    }
    
    void setGrizzlyBarDimensions(Mdouble height, Mdouble widthBack, Mdouble widthFront, Mdouble length, Mdouble spacing)
    {
        grizzlyBarHeight_ = height;
        grizzlyBarWidthBack_ = widthBack;
        grizzlyBarWidthFront_ = widthFront;
        grizzlyBarLength_ = length;
        grizzlyBarSpacing_ = spacing;
    }
    
    void setNumberOfSegments(int numX, int numY, int numZ)
    {
        numSegmentsX_ = numX;
        numSegmentsY_ = numY;
        numSegmentsZ_ = numZ;
    }
    
    void setWearCoefficient(Mdouble wearCoefficient)
    { wearCoefficient_ = wearCoefficient; }
    
    void setHardness(Mdouble hardness)
    { hardness_ = hardness; }
    
    void setWearAcceleration(Mdouble wearAcceleration)
    { wearAcceleration_ = wearAcceleration; }
    
    void setParticleRadiusRange(Mdouble min, Mdouble max)
    {
        logger.assert_always(max >= min, "Maximum particle radius must be greater or equal to minimum particle radius!");
        minParticleRadius_ = min;
        maxParticleRadius_ = max;
    }
    
private:
    Vec3D screenRestPosition_;
    Mdouble screenAngle_, screenSpacing_, vibrationAmplitude_, vibrationFrequency_;
    int numberOfGrizzlyBars_;
    Mdouble grizzlyBarHeight_, grizzlyBarWidthBack_, grizzlyBarWidthFront_, grizzlyBarLength_, grizzlyBarSpacing_;
    int numSegmentsX_, numSegmentsY_, numSegmentsZ_;
    Mdouble wearCoefficient_, hardness_, wearAcceleration_;
    Mdouble minParticleRadius_, maxParticleRadius_;
};

int main(int argc, char** argv)
{
    Screen problem;
    problem.setName("Screen_008");
    
    problem.setTimeStep(1e-5);
    problem.setTimeMax(1.0);
    problem.setSaveCount((1.0 / problem.getTimeStep() / 30)); // 30 fps
    problem.setMin(Vec3D(0.0, 0.0, 0.0));
    problem.setMax(Vec3D(1.0, 1.0, 1.0)); // YMax is overridden in setup to fit tightly around screen, depending on number of grizzly bars used.
    problem.setGravity(Vec3D(0.0, 0.0, -9.81));
    
    problem.setParticlesWriteVTK(true);
    problem.wallHandler.setWriteVTK(true);
    problem.wallHandler.setWriteDetailsVTK(WallHandler::DetailsVTKOptions::BOUNDINGBOX, FileType::ONE_FILE);
    
    // For parallel computing
    problem.setNumberOfOMPThreads(1);
    
    // Custom
    problem.setScreenProperties(Vec3D(0.0, 0.0, 0.5), 15.0 * constants::pi / 180.0, 50.0e-3, 5);
    problem.setVibrationProperties(0.01, 3.0);
    problem.setGrizzlyBarDimensions(12.0e-3, 40.0e-3, 20.0e-3, 750.0e-3, 80.0e-3);
    problem.setNumberOfSegments(180, 5, 3);
    problem.setWearCoefficient(1.0e-6);
    problem.setHardness(1.0);
    problem.setWearAcceleration(1.0);
    problem.setParticleRadiusRange(20.0e-3, 30.0e-3);
    
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();
    
    problem.solve(argc, argv);
    
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Total computation time: % s", elapsed_seconds.count());
    
    return 0;
}