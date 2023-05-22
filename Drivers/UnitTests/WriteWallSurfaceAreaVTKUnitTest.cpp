#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Walls/InfiniteWall.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>

class AreaVTK : public Mercury3D {
public:
    AreaVTK(std::string const& wallSource)
        : wallSource_(wallSource)
    {
        species_ = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

    void setupInitialConditions() override
    {
        setName(wallSource_);

        setTimeStep(1);
        setTimeMax(1);

        removeOldFiles();
        setSaveCount(1);
        setFileType(FileType::NO_FILE);
        setWritePythonFileForVTKVisualisation(false);
        setParticlesWriteVTK(false);
        wallHandler.setWriteVTK(true);
        wallHandler.setWriteWallSurfaceAreaVTK(true);

        setMin({ -1, -1, -1 });
        setMax({ 1, 1, 1 });

        if (wallSource_ == "InfiniteWall") {
            InfiniteWall w0;
            w0.setSpecies(species_);
            w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, 0));
            wallHandler.copyAndAddObject(w0);
        } else {
            wallHandler.readTriangleWall(wallSource_, species_);
        }
    }

    double getTotalsurfaceAreaFromFile(std::string const& fileName)
    {
        FILE* fp = fopen(fileName.c_str(), "r");
        if (fp == NULL) {
            logger(ERROR, "%", strerror(errno));
        }

        auto trimLine = [this](char* line) {
            char* start = line;
            while (isspace(*start)) {
                start++;
            }

            char* end = line + strlen(line) - 1;

            while (isspace(*end)) {
                end--;
            }

            int length = end - start + 1;
            if (length <= 0) {
                logger(ERROR, "Unable to trim `%`", line);
            }

            return std::string(start, length);
        };

        auto getDouble = [this](std::string const& str) {
            char* next;
            double val = strtod(str.c_str(), &next);
            if ((next == str) || (*next != '\0')) {
                logger(ERROR, "Invalid number: `%`", str);
            };
            return val;
        };

        char line[1024];
        bool inSurfaceArea = false;
        double totalSurfaceArea = 0;
        size_t numberOfSurfaces = 0;

        while (fgets(line, sizeof(line), fp)) {
            if (strstr(line, "SurfaceArea") != NULL) {
                inSurfaceArea = true;
            } else if (inSurfaceArea && strstr(line, "</DataArray>") != NULL) {
                inSurfaceArea = false;
                break;
            } else if (inSurfaceArea) {
                numberOfSurfaces++;
                double surfaceArea = getDouble(trimLine(line));
                totalSurfaceArea += surfaceArea;
            } else {
                // do nothing
            }
        }

        if (numberOfSurfaces == 0) {
            logger(ERROR, "No 'SurfaceArea' found in %", fileName);
        }

        if (inSurfaceArea) {
            logger(ERROR, "Corrupted 'SurfaceArea' in %", fileName);
        }

        fclose(fp);

        return totalSurfaceArea;
    }

private:
    LinearViscoelasticSpecies* species_;
    std::string wallSource_;
};

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "C");

    {
        AreaVTK problem("InfiniteWall");
        problem.solve(argc, argv);

        double A0 = problem.getTotalsurfaceAreaFromFile(problem.getName() + "Wall_0.vtu");
        double A1 = problem.getTotalsurfaceAreaFromFile(problem.getName() + "Wall_1.vtu");

        double checkedArea = 4.0;   // calculated from domain size
        helpers::check(A0, A1, 0, "Surface areas of first object are equal");
        helpers::check(A0, checkedArea, 1e-3, "Surface area of first object");
    }

    {
        AreaVTK problem("../SelfTests/Walls/Casing.stl");
        problem.solve(argc, argv);

        double A0 = problem.getTotalsurfaceAreaFromFile(problem.getName() + "Wall_0.vtu");
        double A1 = problem.getTotalsurfaceAreaFromFile(problem.getName() + "Wall_1.vtu");

        double checkedArea = 675075.0;  // checked via MeshLab
        helpers::check(A0, A1, 0, "Surface areas of second object are equal");
        helpers::check(A0, checkedArea, 1e-1, "Surface area of second object");
    }

    return 0;
}
