/*
 * RHQInsertionV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */
 
/*  END LEGAL */

#ifndef Hadrons_MRHQ_RHQInsertionV_hpp_
#define Hadrons_MRHQ_RHQInsertionV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            RHQInsertionV                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MRHQ)

class RHQInsertionVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionVPar,
                                    std::string, q,
                                    std::string, gauge);
};

template <typename FImpl, typename GImpl>
class TRHQInsertionV: public Module<RHQInsertionVPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,)
    SINK_TYPE_ALIASES();
public:
    // constructor
    TRHQInsertionV(const std::string name);
    // destructor
    virtual ~TRHQInsertionV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionV, ARG(TRHQInsertionV<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                    TRHQInsertionV implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionV<FImpl, GImpl>::TRHQInsertionV(const std::string name)
: Module<RHQInsertionVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionV<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionV<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionV<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());//, 1, env().getDim(Tp));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionV<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Applying Improvement term V to'" << par().q
                 << std::endl;

    auto &field    = envGet(PropagatorField, par().q);
    const auto &gaugefield  = envGet(GaugeField, par().gauge);
    const auto gauge_t = peekLorentz(gaugefield, 3);

    Gamma gt(Gamma::Algebra::GammaT);
    PropagatorField  insertion = gt * (GImpl::CovShiftForward( gauge_t,3,field) - GImpl::CovShiftBackward( gauge_t,3,field));
     
    auto &out  = envGet(PropagatorField, getName());
    out = insertion;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionV_hpp_