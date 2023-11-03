#include <iostream>
#include <cstdlib>

#include "linalg.hpp"
#include "io.hpp"

#include "cfield.h"
#include "globalfunctions.h"

#include "corbifold.h"

#define CHECKERROR true

using std::cout;
using std::endl;
using std::exit;
using std::vector;



/* ########################################################################################
######   CField()                                                                    ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CField object. No content is specified.           ######
######################################################################################## */
CField::CField()
{
  this->internalIndex.assign(3,0);
  this->Labels.push_back("");
  this->Numbers.push_back(1);

  this->BmLCharge = 0.0;
  this->OriginStandardConfig = -1;
}



/* ########################################################################################
######   ~CField()                                                                   ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CField object.                                     ######
######################################################################################## */
CField::~CField()
{
}



/* ########################################################################################
######   SetInternalIndex(const unsigned &i, const unsigned &j, ...)                 ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i : index of this fields' sector in the member variable "Sectors" in the ######
######          class "COrbifold"                                                    ######
######   2) j : index of this fields' fixed brane in the member variable             ######
######          "FixedBranes" in the class "Sector"                                  ######
######   3) k : index of this fields' state in the member variable "InvariantStates" ######
######          in the class "FixedBrane"                                            ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Set the internal index of this field (specifies the sector, the fixed brane ######
######   and the state of this field).                                               ######
######################################################################################## */
void CField::SetInternalIndex(const unsigned &i, const unsigned &j, const unsigned &k)
{
  this->internalIndex[0] = i;
  this->internalIndex[1] = j;
  this->internalIndex[2] = k;
}



/* ########################################################################################
######   GetDiscreteCharge(const SDiscreteSymmetry &DiscreteSymmetry, ...) const     ######
######                                                                               ######
######   Version: 23.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) DiscreteSymmetry : a SDiscreteSymmetry object specifying the discrete    ######
######                         symmetry                                              ######
######   output:                                                                     ######
######   return value        : this fields' discrete charge with respect to          ######
######                         "DiscreteSymmetry"                                    ######
###########################################################################################
######   description:                                                                ######
######   Compute the discrete charge of this field with respect to the symmetry      ######
######   specified by "DiscreteSymmetry".                                            ######
######################################################################################## */
rational<int> CField::GetDiscreteCharge(const SDiscreteSymmetry &DiscreteSymmetry) const
{
  
  if (DiscreteSymmetry.ChargeOperator.size() == 4)
  {
    rational<int> DiscreteCharge(0);

    for (unsigned k = 1; k < 4; ++k)
      DiscreteCharge += (D2Rat(this->q_sh[k] + this->OsciContribution[k]) * DiscreteSymmetry.ChargeOperator[k]);
    //DiscreteCharge += (D2Rat(this->q_sh[k]) * DiscreteSymmetry.ChargeOperator[k]);

    return DiscreteCharge;
  }
  /*else
  if (DiscreteSymmetry.ChargeOperator.size() == 4)
  {
    rational<int> DiscreteCharge(0);

    for (unsigned k = 1; k < 4; ++k)
      DiscreteCharge += (D2Rat(this->q_sh[k] + this->OsciContribution[k]) * DiscreteSymmetry.ChargeOperator[k]);

    return DiscreteCharge;
}*/

  return this->SGElement * DiscreteSymmetry.ChargeOperator;
}



/* ########################################################################################
######   GetModularWeight(const SModularSymmetry &ModularSymmetry, ...) const        ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ModularSymmetry : a SModularSymmetry object specifying the modular       ######
######                        symmetry                                               ######
######   2) Sectors         : the set of all sectors of the corresponding orbifold   ######
######                        needed because the oscillators are stored there        ######
######   output:                                                                     ######
######   return value       : this fields' transformation property with respect to   ######
######                        "ModularSymmetry"                                      ######
###########################################################################################
######   description:                                                                ######
######   Compute the transformation property of this field with respect to the       ######
######   modular symmetry specified by "ModularSymmetry".                            ######
######################################################################################## */
rational<int> CField::GetModularWeight(const SModularSymmetry &ModularSymmetry, const vector<CSector> &Sectors) const
{
  // begin: contribution of moduli to grav-grav-modular anomaly
  // note: this is experimental! 
  if ((this->Multiplet == LCModulus) || (this->Multiplet == RCModulus))
  {
    const size_t s1 = this->GetNumberOfOscillators(Sectors);
    if (this->GetNumberOfOscillators(Sectors) != 1)
    {
      cout << "  Warning in rational<int> CField::GetModularWeight(...) const : Modulus is not excited by one oscillator. Return 0." << endl;
      return 0;
    }
    const CModedOscillator &Osci = this->GetOscillator(0, Sectors);
    
    // T_i = T_ii or U_m modulus
    if (D2Rat(this->q_sh[Osci.GetIndex()]) != 0)
    {
      // T_i = T_ii modulus
      if (Osci.GetComplex())
      {
        if (Osci.GetIndex() == ModularSymmetry.Index)
          return rational<int>(1,1);
        else
          return rational<int>(-1,1);
      }
      // U_m modulus
      else
        return rational<int>(1,1);
    }
    // T_ij modulus
    else
    {
      if (((this->Multiplet == LCModulus) && !Osci.GetComplex()) || ((this->Multiplet == RCModulus) && Osci.GetComplex()))
      {
        cout << "  Warning in rational<int> CField::GetModularWeight(...) const : U_ij Modulus does not exist. Return 0." << endl;
        return 0;
      }
      for (unsigned i = 0; i < 4; ++i)
      {
        if (D2Rat(this->q_sh[i]) != 0)
        {
          if ((ModularSymmetry.Index == i) || (ModularSymmetry.Index == Osci.GetIndex()))
            return rational<int>(-1,1);
          else
            return rational<int>(1,1);
        }
      }
    }
    return 0;
  }
  // end: contribution of moduli to grav-grav-modular anomaly 

  // q_sh
  rational<int> q_sh = D2Rat(this->q_sh[ModularSymmetry.Index]);
  if (q_sh.denominator() == 1)
    return q_sh;

  rational<int> ModularWeight = ModularSymmetry.Const + (q_sh * ModularSymmetry.ChargeOperator[0]);

  // N^j and N^\bar{j}
  const size_t s1 = this->GetNumberOfOscillators(Sectors);
  for (unsigned i = 0; i < s1; ++i)
  {
    const CModedOscillator &Osci = this->GetOscillator(i, Sectors);
    if (Osci.GetIndex() == ModularSymmetry.Index)
    {
      if (Osci.GetComplex())
        ModularWeight += ModularSymmetry.ChargeOperator[2];
      else
        ModularWeight += ModularSymmetry.ChargeOperator[1];
    }
  }
  return ModularWeight;
}



/* ########################################################################################
######   GetNumberOfOscillators(const vector<CSector> &Sectors) const                ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors   : the set of all sectors of the corresponding orbifold         ######
######                  needed because the oscillators are stored there              ######
######   output:                                                                     ######
######   return value : the number of oscillators acting on this field               ######
###########################################################################################
######   description:                                                                ######
######   Returns the number of oscillators acting on this field.                     ######
######################################################################################## */
size_t CField::GetNumberOfOscillators(const vector<CSector> &Sectors) const
{
  const CFixedBrane &FixedBrane = Sectors[this->internalIndex[0]].GetFixedBrane(this->internalIndex[1]);

  return FixedBrane.GetMasslessLeftMover(FixedBrane.GetInvariantState(this->internalIndex[2]).GetLeftMover().GetIndex()).Excitation.OscillatorIndices.size();
}



/* ########################################################################################
######   &GetOscillator(const unsigned &i, const vector<CSector> &Sectors) const     ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify this fields' i-th oscillator stored in "Sectors"     ######
######   2) Sectors   : the set of all sectors of the corresponding orbifold         ######
######                  needed because the oscillators are stored there              ######
######   output:                                                                     ######
######   return value : constant reference to the i-th oscillator of this field      ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' i-th oscillator.                     ######
######################################################################################## */
const CModedOscillator &CField::GetOscillator(const unsigned &i, const vector<CSector> &Sectors) const
{
  const CSector     &Sector     = Sectors[this->internalIndex[0]];
  const CFixedBrane &FixedBrane = Sector.GetFixedBrane(this->internalIndex[1]);
  const CState      &State      = FixedBrane.GetInvariantState(this->internalIndex[2]);

  const vector<unsigned> &OscillatorIndices = FixedBrane.GetMasslessLeftMover(State.GetLeftMover().GetIndex()).Excitation.OscillatorIndices;

  #ifdef CHECKERROR
  if (i >= OscillatorIndices.size())
  {
    cout << "\n  Warning in const CModedOscillator &CField::GetOscillator(...) const : Index i out of range. Set i = 0." << endl;
    return Sector.GetLM_Oscillator(OscillatorIndices[0]);
  }
  #endif

  return Sector.GetLM_Oscillator(OscillatorIndices[i]);
}



/* ########################################################################################
######   &GetLMWeight(const unsigned &i, const vector<CSector> &Sectors) const       ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify this fields' i-th massless p_sh stored in "Sectors"  ######
######   2) Sectors   : the set of all sectors of the corresponding orbifold         ######
######                  needed because the massless left-movers are stored there     ######
######   output:                                                                     ######
######   return value : constant reference to the i-th momentum p_sh of this field   ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' i-th massless left-moving momentum   ######
######   p_sh.                                                                       ######
######################################################################################## */
const CVector &CField::GetLMWeight(const unsigned &i, const vector<CSector> &Sectors) const
{
  const CFixedBrane      &FixedBrane    = Sectors[this->internalIndex[0]].GetFixedBrane(this->internalIndex[1]);
  const CState           &State         = FixedBrane.GetInvariantState(this->internalIndex[2]);
  const CHalfState       &LeftMover     = State.GetLeftMover();

  #ifdef CHECKERROR
  if (i >= this->WeightIndices.size())
  {
    cout << "\n  Warning in const CVector &CField::GetLMWeight(...) const : Index i out of range. Set i = 0." << endl;
    return FixedBrane.GetMasslessLeftMover(LeftMover.GetIndex()).Weights[this->WeightIndices[0]];
  }
  #endif

  return FixedBrane.GetMasslessLeftMover(LeftMover.GetIndex()).Weights[this->WeightIndices[i]];
}



/* ########################################################################################
######   GetNumberOfRMWeights(const vector<CSector> &Sectors) const                  ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors   : the set of all sectors of the corresponding orbifold needed  ######
######                  because the right-moving momenta q_sh are stored there       ######
######   output:                                                                     ######
######   return value : the number of this fields' massless right-moving momenta q_sh######
###########################################################################################
######   description:                                                                ######
######   Returns the number of this fields' massless right-moving momenta q_sh.      ######
######################################################################################## */
size_t CField::GetNumberOfRMWeights(const vector<CSector> &Sectors) const
{
  return Sectors[this->internalIndex[0]].GetFixedBrane(this->internalIndex[1]).GetInvariantState(this->internalIndex[2]).GetRightMover().Weights.size();
}



/* ########################################################################################
######   &GetRMWeight(const unsigned &i, const vector<CSector> &Sectors) const       ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify this fields' i-th massless q_sh stored in "Sectors"  ######
######   2) Sectors   : the set of all sectors of the corresponding orbifold         ######
######                  needed because the massless right-movers are stored there    ######
######   output:                                                                     ######
######   return value : constant reference to the i-th momentum q_sh of this field   ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' i-th massless right-moving momentum  ######
######   q_sh.                                                                       ######
######################################################################################## */
const CVector &CField::GetRMWeight(const unsigned &i, const vector<CSector> &Sectors) const
{
  const CSector    &Sector     = Sectors[this->internalIndex[0]];
  const CHalfState &RightMover = Sector.GetFixedBrane(this->internalIndex[1]).GetInvariantState(this->internalIndex[2]).GetRightMover();

  #ifdef CHECKERROR
  if (i >= RightMover.Weights.size())
  {
    cout << "\n  Warning in const CVector &CField::GetRMWeight(...) const : Index i out of range. Set i = 0." << endl;
    return Sector.GetMasslessRightMover(RightMover.GetIndex()).Weights[RightMover.Weights[0]];
  }
  #endif

  return Sector.GetMasslessRightMover(RightMover.GetIndex()).Weights[RightMover.Weights[i]];
}



/* ########################################################################################
######   &GetSector(const vector<CSector> &Sectors) const                            ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors   : the set of all sectors of the corresponding orbifold         ######
######   output:                                                                     ######
######   return value : constant reference to this fields' sector                    ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' sector.                              ######
######################################################################################## */
const CSector &CField::GetSector(const vector<CSector> &Sectors) const
{
  return Sectors[this->internalIndex[0]];
}



/* ########################################################################################
######   &GetFixedBrane(const vector<CSector> &Sectors) const                        ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors   : the set of all sectors of the corresponding orbifold         ######
######   output:                                                                     ######
######   return value : constant reference to this fields' fixed point / fixed brane ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' fixed point / fixed brane.           ######
######################################################################################## */
const CFixedBrane &CField::GetFixedBrane(const vector<CSector> &Sectors) const
{
  return Sectors[this->internalIndex[0]].GetFixedBrane(this->internalIndex[1]);
}



/* ########################################################################################
######   &GetState(const vector<CSector> &Sectors) const                             ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors   : the set of all sectors of the corresponding orbifold         ######
######   output:                                                                     ######
######   return value : constant reference to this fields' state                     ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to this fields' state.                               ######
######################################################################################## */
const CState &CField::GetState(const vector<CSector> &Sectors) const
{
  return Sectors[this->internalIndex[0]].GetFixedBrane(this->internalIndex[1]).GetInvariantState(this->internalIndex[2]);
}
  
