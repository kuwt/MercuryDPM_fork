/*
    2D implementation of mixing indices on the following publication:
        Shih-Hao Chou, Yue-Lou Song and Shu-San Hsiau (2017):
        A Study of the Mixing Index in Solid Particles,
        KONA Powder and Particle Journal, No. 34, pp. 275–281
        DOI: 10.14356/kona.2017018

    2023, Tóth János, Hungarian University of Agriculture and Life Sciences, Gödöllő
*/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"

#include <Mercury3D.h>

#include <Particles/SphericalParticle.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Walls/InfiniteWall.h>

#pragma GCC diagnostic pop

#include "Grid2D.h"

class MixingIndices : public Mercury3D {
public:
    MixingIndices(int patternCase)
        : patternCase_(patternCase)
    {
        LinearViscoelasticSpecies species;
        species.setDensity(500);
        species.setStiffness(1e2);
        species.setDissipation(10e-3);
        speciesHandler.copyAndAddObject(species);
    }

    void setupInitialConditions() override
    {
        setName("MixingIndices" + std::to_string(patternCase_) + "_");
        setSystemDimensions(2);
        setGravity({ 0, 0, 0 });

        setSaveCount(100);
        setFileType(FileType::NO_FILE);
        setWritePythonFileForVTKVisualisation(false);
        wallHandler.setWriteVTK(true);
        setParticlesWriteVTK(true);

        setTimeMax(1.0);
        setTimeStep(1.0);

        setXMin(0);
        setXMax(12 * 2 * particleRadius);
        setYMin(0);
        setYMax(12 * 2 * particleRadius);

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set({ 1.0, 0.0, 0.0 }, { getXMax(), 0.0, 0.0 });
        wallHandler.copyAndAddObject(w0);
        w0.set({ -1.0, 0.0, 0.0 }, { getXMin(), 0.0, 0.0 });
        wallHandler.copyAndAddObject(w0);
        w0.set({ 0.0, 1.0, 0.0 }, { 0.0, getYMax(), 0.0 });
        wallHandler.copyAndAddObject(w0);
        w0.set({ 0.0, -1.0, 0.0 }, { 0.0, getYMin(), 0.0 });
        wallHandler.copyAndAddObject(w0);

        insertParticles();

        logger(INFO, "Total particles: %", particleHandler.getSize());
    }

    void insertParticles()
    {
        const double startX = getXMin() + particleRadius;
        const double startY = getYMin() + particleRadius;

        const double endX = getXMax() - particleRadius;

        const size_t maxParticles1D = static_cast<size_t>((endX - startX) / (particleRadius * 2.0)) + 1;
        logger.assert_always(maxParticles1D == 12, "wrong geometry");

        auto case1 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (j > 5) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case2 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (i > 5) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        // NOTE: not the same, but good enough
        auto case3 = [this, maxParticles1D]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (i > j) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case5 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (j < 4) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else if (j >= 4 && j < 8) {
                size_t index = i;
                index += j % 2;
                if (index % 2 == 0) {
                    particleHandler.getLastObject()->setId(TYPE_A + n);
                } else {
                    particleHandler.getLastObject()->setId(TYPE_B + n);
                }
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case6 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (j < 2) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else if (j >= 2 && j < 10) {
                size_t index = i;
                index += j % 2;
                if (index % 2 == 0) {
                    particleHandler.getLastObject()->setId(TYPE_A + n);
                } else {
                    particleHandler.getLastObject()->setId(TYPE_B + n);
                }
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case7 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            size_t index = i;
            index += j % 2;
            if (index % 2 == 0) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case8 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (j % 2 == 0) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case9 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (i % 2 == 0) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        auto case11 = [this]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j, size_t n) {
            if (i == 0 || i == 11 || j == 0 || j == 11) {
                particleHandler.getLastObject()->setId(TYPE_A + n);
            } else {
                particleHandler.getLastObject()->setId(TYPE_B + n);
            }
        };

        size_t n = 0;
        for (size_t i = 0; i < maxParticles1D; ++i) {
            for (size_t j = 0; j < maxParticles1D; ++j) {
                double x = startX + i * particleRadius * 2.0;
                double y = startY + j * particleRadius * 2.0;

                SphericalParticle p;
                p.setSpecies(speciesHandler.getObject(0));
                p.setRadius(particleRadius);
                p.setPosition({ x, y, 0.0 });
                particleHandler.copyAndAddObject(p);

                switch (patternCase_) {
                case 1:
                    case1(i, j, n);
                    break;
                case 2:
                    case2(i, j, n);
                    break;
                case 3:
                    case3(i, j, n);
                    break;
                case 5:
                    case5(i, j, n);
                    break;
                case 6:
                    case6(i, j, n);
                    break;
                case 7:
                    case7(i, j, n);
                    break;
                case 8:
                    case8(i, j, n);
                    break;
                case 9:
                    case9(i, j, n);
                    break;
                case 11:
                    case11(i, j, n);
                    break;
                default:
                    logger(ERROR, "Invalid pattern case: %", patternCase_);
                }

                ++n;
            }
        }
    }

    double getLaceyIndex(double gridSize)
    {
        Grid2D g(getXMin(), getXMax(), getYMin(), getYMax(), gridSize);

        for (auto& p : particleHandler) {
            auto pos = p->getPosition();
            g.addParticle(pos.getX(), pos.getY(), p->getId());
        }

        logger(INFO, "Total: %", g.getNumberOfParticles());

        return g.getLaceyIndex();
    }

    double getKramerIndex(double gridSize)
    {
        Grid2D g(getXMin(), getXMax(), getYMin(), getYMax(), gridSize);

        for (auto& p : particleHandler) {
            auto pos = p->getPosition();
            g.addParticle(pos.getX(), pos.getY(), p->getId());
        }

        logger(INFO, "Total: %", g.getNumberOfParticles());

        return g.getKramerIndex();
    }

    static constexpr double particleRadius = (1.0 / 10.0) / 2.0;

private:
    int patternCase_ { 1 };
};

int main(int argc, char** argv)
{
    std::vector<std::tuple<int, double, double>> patternMixingIndices;

    const double r = MixingIndices::particleRadius;

    double gridSize = helpers::readFromCommandLine(argc, argv, "-grid", r * 4.0);
    helpers::removeFromCommandline(argc, argv, "-grid", 1);

    for (auto const& pattern : { 1, 2, 3, 5, 6, 7, 8, 9, 11 }) {
        MixingIndices problem(pattern);
        problem.solve(argc, argv);
        const double lacy = problem.getLaceyIndex(gridSize);
        const double kramer = problem.getKramerIndex(gridSize);
        patternMixingIndices.emplace_back(std::make_tuple(pattern, lacy, kramer));
    }

    for (auto const& [pattern, lacey, kramer] : patternMixingIndices) {
        logger(INFO, "Pattern %: %, %", pattern, lacey, kramer);
    }

    return 0;
}
