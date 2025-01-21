#include <Hadrons/Modules/MContraction/QEDBurgerShort.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

template class Grid::Hadrons::MContraction::TQEDBurgerShort<FIMPL, FIMPL::PropagatorField, vComplex>;
template class Grid::Hadrons::MContraction::TQEDBurgerShort<FIMPL, FIMPL::FermionField,    vComplex>;
