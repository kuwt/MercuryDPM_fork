/*
    Cell/Grid structures for the standard deviation based mixing indices.

    2023, Tóth János, Hungarian University of Agriculture and Life Sciences, Gödöllő
*/

#pragma once

#include <Mercury3D.h>

static constexpr double TYPE_A = 100'000;
static constexpr double TYPE_B = 200'000;

class Cell2D {
public:
    Cell2D(double x, double y, double size)
        : x_(x)
        , y_(y)
        , size_(size) {};

    bool contains(double x, double y)
    {
        return x >= x_ && x < x_ + size_ && y >= y_ && y < y_ + size_;
    }

    void addType(unsigned int id)
    {
        if (id >= TYPE_A && id < TYPE_B) {
            numberOfA_++;
        } else if (id >= TYPE_B) {
            numberOfB_++;
        }
    }

    size_t getTotalNumber() const
    {
        return numberOfA_ + numberOfB_;
    }

    size_t getNumberOfA() const
    {
        return numberOfA_;
    };

    size_t getNumberOfB() const
    {
        return numberOfB_;
    };

    double getFractionOfA() const
    {
        return (double)(numberOfA_) / (double)(numberOfA_ + numberOfB_);
    }

    double getX() const
    {
        return x_;
    }

    double getY() const
    {
        return y_;
    }

    double getSize() const
    {
        return size_;
    }

private:
    double x_;
    double y_;
    double size_;

    size_t numberOfA_ { 0 };
    size_t numberOfB_ { 0 };
};

class Grid2D {
public:
    Grid2D(double xMin, double xMax, double yMin, double yMax, double size)
    {
        for (double x = xMin; x <= xMax + size; x += size) {
            for (double y = yMin; y <= yMax + size; y += size) {
                cells_.emplace_back(Cell2D(x, y, size));
            }
        }

        logger(INFO, "Grid: [%, %]x[%, %](%)", xMin, xMax, yMin, yMax, size);
    }

    void addParticle(double x, double y, unsigned int id)
    {
        for (auto& c : cells_) {
            if (c.contains(x, y)) {
                c.addType(id);
                return;
            }
        }
    }

    size_t getNumberOfParticles() const
    {
        size_t total = 0;
        for (auto const& c : cells_) {
            total += c.getTotalNumber();
        }

        return total;
    }

    double getLaceyIndex() const
    {
        double totalCells = 0;
        double totalParticles = 0.0;
        double totalA = 0.0;
        for (auto const& c : cells_) {
            if (c.getTotalNumber() > 0) {
                totalCells += 1.0;
                totalParticles += c.getTotalNumber();
                totalA += c.getNumberOfA();
            }
        }

        if (totalCells < 1.0) {
            logger(ERROR, "invalid mixing");
        }

        const double N = totalCells;
        const double n = totalParticles / totalCells;
        const double phi_m = totalA / totalParticles;

        double sigma2 = 0;
        for (auto const& c : cells_) {
            if (c.getTotalNumber() > 0) {
                const double phi_i = c.getFractionOfA();
                sigma2 += ((phi_i - phi_m) * (phi_i - phi_m)) / (N - 1);
            }
        }

        const double sigma_02 = phi_m * (1 - phi_m);
        const double sigma_r2 = (phi_m * (1 - phi_m)) / n;

        logger(INFO, "  N = %", N);
        logger(INFO, "  n = %", n);
        logger(INFO, "  phi_m = %", phi_m);
        logger(INFO, "  sigma2 = %", sigma2);
        logger(INFO, "  sigma_02 = %", sigma_02);
        logger(INFO, "  sigma_r2 = %", sigma_r2);

        const double laceyIndex = (sigma_02 - sigma2) / (sigma_02 - sigma_r2);

        if (laceyIndex > 1.0) {
            logger(WARN, "Lacey index is greater than 1.0, setting value to 1.0");
            return 1.0;
        } else if (laceyIndex < 0.0) {
            logger(WARN, "Lacey index is less than 0.0, setting value to 0.0");
            return 0.0;
        } else {
            return laceyIndex;
        }
    }

    double getKramerIndex() const
    {
        double totalCells = 0;
        double totalParticles = 0.0;
        double totalA = 0.0;
        for (auto const& c : cells_) {
            if (c.getTotalNumber() > 0) {
                totalCells += 1.0;
                totalParticles += c.getTotalNumber();
                totalA += c.getNumberOfA();
            }
        }

        if (totalCells < 1.0) {
            logger(ERROR, "invalid mixing");
        }

        const double N = totalCells;
        const double n = totalParticles / totalCells;
        const double phi_m = totalA / totalParticles;

        double sigma2 = 0;
        for (auto const& c : cells_) {
            if (c.getTotalNumber() > 0) {
                const double phi_i = c.getFractionOfA();
                sigma2 += ((phi_i - phi_m) * (phi_i - phi_m)) / (N - 1);
            }
        }

        const double sigma_02 = phi_m * (1 - phi_m);
        const double sigma_r2 = (phi_m * (1 - phi_m)) / n;

        const double sigma = sqrt(sigma2);
        const double sigma_0 = sqrt(sigma_02);
        const double sigma_r = sqrt(sigma_r2);

        logger(INFO, "  N = %", N);
        logger(INFO, "  n = %", n);
        logger(INFO, "  phi_m = %", phi_m);
        logger(INFO, "  sigma = %", sigma);
        logger(INFO, "  sigma_0 = %", sigma_0);
        logger(INFO, "  sigma_r = %", sigma_r);

        const double kramerIndex = (sigma_0 - sigma) / (sigma_0 - sigma_r);

        if (kramerIndex > 1.0) {
            logger(WARN, "Kramer index is greater than 1.0, setting value to 1.0");
            return 1.0;
        } else if (kramerIndex < 0.0) {
            logger(WARN, "Kramer index is less than 0.0, setting value to 0.0");
            return 0.0;
        } else {
            return kramerIndex;
        }
    }

private:
    std::vector<Cell2D> cells_;
};
