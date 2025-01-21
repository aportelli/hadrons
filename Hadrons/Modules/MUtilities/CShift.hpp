/*
 * CShift.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Matthew Black <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MUtilities_CShift_hpp_
#define Hadrons_MUtilities_CShift_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Field CShift                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class CShiftPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CShiftPar,
                                    std::string, field,
                                    std::string, shift);
};

template <typename Field>
class TCShift: public Module<CShiftPar>
{
public:
    // constructor
    TCShift(const std::string name);
    // destructor
    virtual ~TCShift(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CShiftGaugeField, TCShift<GIMPL::GaugeField>, MUtilities);
MODULE_REGISTER_TMP(CShiftPropagatorField, TCShift<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                     TCShift implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TCShift<Field>::TCShift(const std::string name)
: Module<CShiftPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TCShift<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename Field>
std::vector<std::string> TCShift<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TCShift<Field>::setup(void)
{
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TCShift<Field>::execute(void)
{
    Coordinate xvec = strToVec<int>(par().shift);
    LOG(Message) << "Shifting '" << par().field << "' by vector [" 
                 << xvec[0] << " " << xvec[1] << " "
                 << xvec[2] << " " << xvec[3] << "]." << std::endl;

    auto &U   = envGet(Field, par().field);
    auto &Utr = envGet(Field, getName());
    Utr = U;
    
    for (int dir = 0; dir < Nd; dir++) {
        int coord = xvec[dir];
        if (coord != 0) Utr = Cshift(Utr,dir,coord);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_CShift_hpp_
