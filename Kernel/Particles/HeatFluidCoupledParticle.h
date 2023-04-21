//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef HeatFluidCoupledParticle_H
#define HeatFluidCoupledParticle_H

#include "ThermalParticle.h"

/*!
 * \class HeatFluidCoupledParticle
 * \brief Class that implements particles which store both temperature/heat capacity and liquid content
 *        which is adapted for the CFD-DEM studies.
 * \details In some CFD-DEM studies, a drying process needs to be simulated.
 *          Therefore, we selected the drying model of Azmir et al. (2018), which considers heat and mass transfer.
 *          Their model considers convection, conduction, radiation, and evaporation.
 *          Conduction and convection had already been implemented in ThermalParticle class.
 *          This class implements evaporation model.
 */
class HeatFluidCoupledParticle final : public ThermalParticle
{
public:
    /*!
     * \breif HeatFluidCoupledParticle constructor creates a HeatFluidCoupledParticle at (0,0,0) with radius,
     *        mass and inertia equal to 1.
     */
    HeatFluidCoupledParticle()
    {
        liquidVolume_ = 0;
    }

    /*!
     * \brief HeatFluidCoupledParticle copy constructor, which accepts as input a reference to a HeatFluidCoupledParticle.
     *        It creates a copy of this HeatFluidCoupledParticle and all it's information.
     *        Usually it is better to use the copy() function for polymorphism.
     * \details Constructor that copies most of the properties of the given particle.
     *          Please note that not everything is copied, for example the position
     *          in the HGrid is not determined yet by the end of this constructor.
     *          It also does not copy the interactions and the pointer to the handler
     *          that handles this particle. Use with care.
     * \param[in,out] p  Reference to the HeatFluidCoupledParticle this one should become a copy of.
     */
    HeatFluidCoupledParticle(const HeatFluidCoupledParticle& p)
    {
        liquidVolume_ = p.liquidVolume_;
    }

    /*!
     * \brief HeatFluidCoupledParticle destructor, needs to be implemented and checked if it removes tangential spring information.
     * \details Destructor that asks the ParticleHandler to check if this was the smallest or largest particle and adjust itself accordingly.
     */
    ~HeatFluidCoupledParticle() override
    = default;

    /*!
     * \brief HeatFluidCoupledParticle copy method. Use copy constructor of this HeatFluidCoupledParticle to create a copy on the heap,
     *        useful for polymorphism.
     *\return pointer to the particle's copy.
     */
    HeatFluidCoupledParticle* copy() const override
    {
        return new HeatFluidCoupledParticle(*this);
    }

    /*!
     * \brief HeatFluidCoupledParticle write function, writes HeatFluidCoupledParticle information to the given output-stream,
     *        for example a restart-file.
     * \details HeatFluidCoupledParticle print method, which accepts an os std::ostream as input.
     *          It prints human readable HeatFluidCoupledParticle information to the std::ostream.
     * \param[in,out] os    stream to which the info is written.
     */
    void write(std::ostream& os) const override
    {
        ThermalParticle::write(os);
        os << " liquidVolume " << liquidVolume_;
    }

    /*!
     * \breif Returns the name of the object; in this case "HeatFluidCoupledParticle".
     * \return The object name.
     */
    std::string getName() const override
    {
        return "HeatFluidCoupledParticle";
    }

    /*!
     * \brief HeatFluidCoupledParticle read function, reads in the information for this HeatFluidCoupledParticle from the given input-stream,
     *        for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Returns the volume of the Liquid.
     * \return The actual volume of the liquid.
     */
    Mdouble getLiquidVolume() const
    {
        return liquidVolume_;
    }

    /*!
     * \brief Sets the volume of the Liquid.
     */
    void setLiquidVolume(Mdouble liquidVolume)
    {
        liquidVolume_ = liquidVolume;
    }

    /*!
     * \brief Adds the volume of the Liquid.
     */
    void addLiquidVolume(Mdouble liquidVolume)
    {
        liquidVolume_ += liquidVolume;
    }

    unsigned getNumberOfFieldsVTK() const override
    {
        return 4;
    }

    std::string getTypeVTK(unsigned i) const override
    {
        return "Float32";
    }

    std::string getNameVTK(unsigned i) const override;

    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    bool isSphericalParticle() const override {return true;}

    void actionsAfterTimeStep() override;

private:

    /// f1 is used in Runge–Kutta method.
    double f1(double liquidVolume,double temperature);

    /// f2 is used in Runge–Kutta method.
    double f2(double liquidVolume,double temperature);

    //Volume of the liquid
    Mdouble liquidVolume_;

};

#endif