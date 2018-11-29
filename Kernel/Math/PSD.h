#ifndef MERCURY_PSD_H
#define MERCURY_PSD_H

#include <fstream>

/**
 * Stores a radius and a cumulative number density:
 * To be used as a vector, std::vector<PSD> psd
 * For each particle, there is the probability p_i that its radius r is less than r_i.
 *     p(r<r_i)=n_i.
 */
struct PSD
{
   
    /**
     * required to use std::lower_bound for finding when the probability is higher than a certain value
     *   std::vector<PSD> psd;
     *   double probability;
     *   std::lower_bound(psd.begin(),psd.end(),probability);
     * @param probability
     * @return
     */
    bool operator<(const double probability) const
    {
        return this->probability < probability;
    }
    
    /*!
     * \brief Writes to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const PSD& psd) {
        os << psd.radius << ' ' << psd.probability;
        return os;
    }
    
    /*!
     * \brief Reads from input stream
     */
    friend std::istream& operator>>(std::istream& is, PSD& psd) {
        is >> psd.radius >> psd.probability;
        return is;
    }
    
    double radius;
    double probability;
};



#endif //MERCURY_PSD_H
