#include "cyukawacouplings.h"
#include "corbifold.h"
#include "globalfunctions.h"
#include "cgaugeinvariance.h"
#include "cprint.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;
using std::vector;

bool HWcriterion(const vector<int> &HW1, const vector<int> &HW2);

CVector Convert_complexVector2CVector(const complexVector &input)
{
  const size_t size = input.size();
  CVector result(2*size);
  for (unsigned i = 0; i < size; ++i)
  {
    result[2*i]      = input[i].real();
    result[2*i + 1 ] = input[i].imag();
  }
  return result;
}

/* ########################################################################################
######   CYukawaCouplings::CYukawaCouplings(...)                                     ######
######                                                                               ######
######   Version: 25.5.2009                                                          ######
######################################################################################## */
CYukawaCouplings::CYukawaCouplings()
{
  this->tmp_file = NULL;
}


/* ########################################################################################
######   CYukawaCouplings::~CYukawaCouplings()                                       ######
######################################################################################## */
CYukawaCouplings::~CYukawaCouplings()
{
}



/* ########################################################################################
######   void CYukawaCouplings::CreatePreCouplings(...)                              ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   output:                                                                     ######
######################################################################################## */
void CYukawaCouplings::CreatePreCouplings(const COrbifold &Orbifold, const SConfig &Vacuum, const vector<unsigned> &MaxDigits, const vector<unsigned> &vec_adjoining_Positions, const vector<vector<unsigned> > &PoolsOfFields, unsigned PoolIndex, const vector<unsigned> &NumbersOfFieldsFromPool)
{
  const double   prec      = 0.0001;
  const unsigned Order     = MaxDigits.size();
  const unsigned OrderMOne = Order - 1;

  const vector<SDiscreteSymmetry> &DiscreteRSymmetries    = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteRSymmetries;
  const vector<SDiscreteSymmetry> &DiscreteNonRSymmetries = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteNonRSymmetries;

  const size_t number_of_R_symmetries    = DiscreteRSymmetries.size();
  const size_t number_of_nonR_symmetries = DiscreteNonRSymmetries.size();

  unsigned i = 1;
  unsigned j = 1;
  bool change = false;

  unsigned number_of_X = 0;
  double Sum  = 0.0; // for U(1) charges
  rational<int> RSum(0); // for R-charges
  unsigned tmp_index = 0;
  bool coupling_allowed = true;

  const vector<CField> &Fields = Vacuum.Fields;

  vector<unsigned> FieldCoupling(Order, 0);
  vector<unsigned> currentNumber(Order, 0);

  const size_t s1 = PoolsOfFields.size();
  while (true)
  {
    // begin: create the vector of field indices ("FieldCoupling") from "currentNumber"
    tmp_index = 0;
    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &CurrentPoolOfFields = PoolsOfFields[i];

      for (j = 0; j < NumbersOfFieldsFromPool[i]; ++j)
      {
        FieldCoupling[tmp_index] = CurrentPoolOfFields[currentNumber[tmp_index]];
        ++tmp_index;
      }
    }
      
    // end: create the vector of field indices ("FieldCoupling") from "currentNumber"
    coupling_allowed = true;

    // begin: U1 gauge invariance
    number_of_X = Fields[0].U1Charges.size(); // for the number of U(1)s and the number of gammas

    for (i = 0; i < number_of_X; ++i)
    {
      Sum = 0.0;
      for (j = 0; j < Order; ++j)
        Sum += Fields[FieldCoupling[j]].U1Charges[i];

      if (fabs(Sum) > prec)
      {
        coupling_allowed = false;
        break;
      }
    }
    // end: U1 gauge invariance

    if (coupling_allowed)
    {
      // begin: G A M M A  selection rule
      number_of_X = Fields[0].gamma_phases.size();

      /*for (i = 0; i < number_of_X; ++i)
      {
        Sum = 0.0;
        for (j = 0; j < Order; ++j)
          Sum += Fields[FieldCoupling[j]].gamma_phases[i];

        if (!is_integer(Sum))
        {
          coupling_allowed = false;
          break;
        }
    }*/
      // end: G A M M A  selection rule

      if (coupling_allowed)
      {
        // begin: R-charge selection rule
        for (i = 0; i < number_of_R_symmetries; ++i)
        {
          const SDiscreteSymmetry &DiscreteRSymmetry = DiscreteRSymmetries[i];

          RSum = 0;
          for (j = 0; j < Order; ++j)
            RSum += Fields[FieldCoupling[j]].GetDiscreteCharge(DiscreteRSymmetry);

          if (RSum != DiscreteRSymmetry.SuperpotentialCharge)
          {
            RSum -= DiscreteRSymmetry.SuperpotentialCharge;
            if (RSum.denominator() != 1)
            {
              coupling_allowed = false;
              break;
            }
            else
            {
              if (RSum.numerator() % (int)(DiscreteRSymmetry.Order) != 0)  
              {
                coupling_allowed = false;
                break;
              }
            }
          }
        }
        // end: R-charge selection rule

        if (coupling_allowed)
        {
          // begin: space-group selection rule
          for (i = 0; i < number_of_nonR_symmetries; ++i)
          {
            const SDiscreteSymmetry &DiscreteNonRSymmetry = DiscreteNonRSymmetries[i];

            RSum = 0;
            for (j = 0; j < Order; ++j)
              RSum += Fields[FieldCoupling[j]].GetDiscreteCharge(DiscreteNonRSymmetry);

            if (RSum.denominator() != 1)
            {
              coupling_allowed = false;
              break;
            }
            else
            {
              if (RSum.numerator() % DiscreteNonRSymmetry.Order != 0)
              {
                coupling_allowed = false;
                break;
              }
            }
          }
          // end: space-group selection rule

          // N O N - A B E L I A N  G A U G E Invariance selection rule
          if (coupling_allowed && this->GI.CheckInvariance(Orbifold, Vacuum, FieldCoupling))
          {
            // begin: print allowed pre-coupling
            (*this->tmp_file) << PoolIndex << "  ";
            for (i = 0; i < Order; ++i)
              (*this->tmp_file) << currentNumber[i] << " ";
            (*this->tmp_file) << "\n";
            // end: print allowed pre-coupling
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // COUNTER
    //////////////////////////////////////////////////////////////////////////////////////////
    // begin: increase the backmost entry if possible
    change = false;
    for (i = OrderMOne; i > 0; --i)
    {
      // neighbour entries are correlated
      if (vec_adjoining_Positions[i] == vec_adjoining_Positions[i-1])
      {
        if ((currentNumber[i] < MaxDigits[i]) && (currentNumber[i] < currentNumber[i-1]))
        {
          ++currentNumber[i];
          for (j = i+1; j < Order; ++j)
            currentNumber[j] = 0;
          change = true;
          break;
        }
      }
      else
      // neighbour entries are uncorrelated
      {
        if (currentNumber[i] < MaxDigits[i])
        {
          ++currentNumber[i];
          for (j = i+1; j < Order; ++j)
            currentNumber[j] = 0;
          change = true;
          break;
        }
      }
    }
    // end: increase the backmost entry if possible

    // begin: increase the first entry if possible
    if (!change)
    {
      // STOP
      if (currentNumber[0] == MaxDigits[0])
        return;
      ++currentNumber[0];
      for (j = 1; j < Order; ++j)
        currentNumber[j] = 0;
    }
    // end: increase the first entry if possible
    //////////////////////////////////////////////////////////////////////////////////////////
    // COUNTER
    //////////////////////////////////////////////////////////////////////////////////////////
  }
}



bool length_criterion(const CVector &Vector1, const CVector &Vector2)
{
  if (Vector1.size() != Vector2.size())
  {
    cout << "\n  Warning in bool length_criterion(...) : Return false." << endl;
    return false;
  }
  const double length1 = Vector1 * Vector1;
  const double length2 = Vector2 * Vector2;

  if (length1 < length2 - 0.0001)
    return true;
  else
  {
    if (fabs(length1 - length2) < 0.0001)
    {
      const size_t s = Vector1.size();
      for (unsigned i = 0; i < s; ++i)
      {
        if (Vector1[i] < Vector2[i] - 0.0001)
          return true;
        else
        if (Vector1[i] > Vector2[i] + 0.0001)
          return false;
      }
    }
    else
    if (length1 > length2 + 0.0001)
      return false;
  }
  return true;
}



/* ########################################################################################
######   bool CYukawaCouplings::AutoCreateMassMatrix(...)                            ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CYukawaCouplings::AutoCreateMassMatrix(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &Labels, unsigned maxOrderOfSinglets, const string &SingletLabel, const vector<string> *AllowedFields)
{
  const size_t l1 = Labels.size();
  if (l1 < 2)
  {
    cout << "\n  Warning in bool CYukawaCouplings::AutoCreateMassMatrix(...): \"Labels\" needs to contain at least two field labels. Return false." << endl;
    return false;
  }

  unsigned tmpcounter = 0;
  
  // number of singlets "SingletLabel"
  const int number_of_Singlets = maxOrderOfSinglets - l1 + 2;
  
  // begin: create upto two singlets for all entries
  vector<string> UseLabels = Labels;

  // no additional singlet
  if (l1 != 2)
  {
    tmpcounter = this->AddCoupling(Orbifold, Vacuum, UseLabels, AllowedFields);

    if (tmpcounter != 0)
    {
      for (unsigned j = 0; j < UseLabels.size(); ++j)
        cout << UseLabels[j] << " ";
      cout << " -> " << tmpcounter << endl;
    }
  }
  
  if (number_of_Singlets == 0)
    return true;

  // begin: add one singlet "SingletLabel"
  UseLabels.push_back(SingletLabel);
  tmpcounter = this->AddCoupling(Orbifold, Vacuum, UseLabels, AllowedFields);

  if (tmpcounter != 0)
  {
    for (unsigned j = 0; j < UseLabels.size(); ++j)
      cout << UseLabels[j] << " ";
    cout << " -> " << tmpcounter << endl;
  }
  // end: add one singlet "SingletLabel"

  // next, add 2 or more singlets "SingletLabel"
  unsigned Singlet_counter = 2;
      
  // create the mass matrix
  CMassMatrix Test_MassMatrix(Vacuum, UseLabels[0], UseLabels[1], true);

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  const vector<unsigned> &FieldIndices_Row    = Test_MassMatrix.FieldIndices_Row;
  const vector<unsigned> &FieldIndices_Column = Test_MassMatrix.FieldIndices_Column;

  const size_t number_of_columns = FieldIndices_Column.size();
  const size_t number_of_rows    = FieldIndices_Row.size();

  const string Label_Row    = Test_MassMatrix.Label_Row;
  const string Label_Column = Test_MassMatrix.Label_Column;

  string tmp_Label1 = "";
  string tmp_Label2 = "";
  for (i = 0; i < number_of_rows; ++i)
  {
    std::ostringstream os1;
    os1 << Vacuum.Fields[FieldIndices_Row[i]].Numbers[Vacuum.use_Labels];
    tmp_Label1 = Label_Row + "_" + os1.str();

    for (j = 0; j < number_of_columns; ++j)
    {
      std::ostringstream os2;
      os2 << Vacuum.Fields[FieldIndices_Column[j]].Numbers[Vacuum.use_Labels];
      tmp_Label2 = Label_Column + "_" + os2.str();

      Singlet_counter = 2;
      while (Test_MassMatrix.GetMatrix()[i][j].size() == 0)
      {
        if (Singlet_counter > number_of_Singlets)
          break;

        UseLabels.clear();
        UseLabels.push_back(tmp_Label1);
        UseLabels.push_back(tmp_Label2);

        for (k = 0; k < Singlet_counter; ++k)
          UseLabels.push_back(SingletLabel);

        tmpcounter = this->AddCoupling(Orbifold, Vacuum, UseLabels, AllowedFields);

        if (tmpcounter != 0)
        {
          for (unsigned j = 0; j < UseLabels.size(); ++j)
            cout << UseLabels[j] << " ";
          cout << " -> " << tmpcounter << endl;
        }
        Test_MassMatrix.Update(Vacuum);

        ++Singlet_counter;
      }
    }
  }
  return true;
}



/* ########################################################################################
######   bool CYukawaCouplings::AutoCreateFullMassMatrix(...)                        ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CYukawaCouplings::AutoCreateFullMassMatrix(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &Labels, unsigned maxOrderOfSinglets, const string &SingletLabel, const vector<string> *AllowedFields)
{
  const size_t l1 = Labels.size();
  if (l1 < 2)
  {
    cout << "\n  Warning in bool CYukawaCouplings::AutoCreateFullMassMatrix(...): \"Labels\" needs to contain at least two field labels. Return false." << endl;
    return false;
  }

  unsigned tmpcounter = 0;
  
  // number of singlets "SingletLabel"
  const int number_of_Singlets = maxOrderOfSinglets - l1 + 2;
  
  vector<string> UseLabels = Labels;

  // no additional singlet
  if (l1 != 2)
  {
    tmpcounter = this->AddCoupling(Orbifold, Vacuum, UseLabels, AllowedFields);

    if (tmpcounter != 0)
    {
      for (unsigned j = 0; j < UseLabels.size(); ++j)
        cout << UseLabels[j] << " ";
      cout << " -> " << tmpcounter << endl;
    }
  }
  
  if (number_of_Singlets == 0)
    return true;

  for (unsigned i = 0; i < number_of_Singlets; ++i)
  {
    UseLabels.push_back(SingletLabel);
    tmpcounter = this->AddCoupling(Orbifold, Vacuum, UseLabels, AllowedFields);
    
    if (tmpcounter != 0)
    {
      for (unsigned j = 0; j < UseLabels.size(); ++j)
        cout << UseLabels[j] << " ";
      cout << " -> " << tmpcounter << endl;
    }
  }
  return true;
}



bool HWcriterion(const vector<int> &HW1, const vector<int> &HW2)
{
  const size_t s1 = HW1.size();
  if (s1 != HW2.size())
  {
    cout << "\n  Warning in bool HWcriterion(...) : Return false." << endl;
    return false;
  }
  for (unsigned i = 0; i < s1; ++i)
  {
    if (HW1[i] > HW2[i])
      return true;
    if (HW1[i] < HW2[i])
      return false;
  }
  return false;
}



/* ########################################################################################
######   unsigned CYukawaCouplings::AddCoupling(...)                                 ######
######                                                                               ######
######   Version: 24.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
unsigned CYukawaCouplings::AddCoupling(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &FieldLabels, const vector<string> *AllowedFields)
{
  const unsigned Order = FieldLabels.size();
  if (Order < 3)
    return 0;

  const vector<CField> &Fields = Vacuum.Fields;
  const size_t f1 = Fields.size();
  if (f1 == 0)
    return 0;

  const COrbifoldGroup  &OrbifoldGroup  = Orbifold.OrbifoldGroup;
  const CSpaceGroup     &SpaceGroup     = OrbifoldGroup.GetSpaceGroup();
  
  const vector<SModularSymmetry> &ModularSymmetries = SpaceGroup.ModularSymmetries;
  const size_t number_of_modularsymmetries = ModularSymmetries.size();
  
  rational<int> exponent = -1;

  const unsigned M = SpaceGroup.GetM();
  const unsigned N = SpaceGroup.GetN();

  // Set the precision
  const double prec = 0.0001;

  // begin: ERROR
  if (Vacuum.use_Labels >= Fields[0].Labels.size())
  {
    cout << "\n  Warning in unsigned CYukawaCouplings::AddCoupling(...): \"use_Labels\" out of range. Return 0." << endl;
    return 0;
  }
  // end: ERROR

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  string tmp_Label1 = "";    // label with number,    e.g.  bL_13
  string tmp_Label2 = "";    // label without number, e.g.  bL

  size_t s1 = 0;
  size_t s2 = 0;
  bool found = false;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: how many different fields?

  // "PreNumbersOfFields[i]" countes the number the field "PreDifferentFields[i]" appears in the coupling
  vector<string>   PreDifferentFields;
  vector<unsigned> PreNumbersOfFields;

  // run through the "FieldLabels"
  for (i = 0; i < Order; ++i)
  {
    const string &FieldLabel = FieldLabels[i];

    found = false;

    s1 = PreDifferentFields.size();
    for (j = 0; j < s1; ++j)
    {
      if (PreDifferentFields[j] == FieldLabel)
      {
        ++PreNumbersOfFields[j];
        found = true;
        break;
      }
    }
    // if "FieldLabel" is unknown, add the new label with corresponding counter 1
    if (!found)
    {
      PreDifferentFields.push_back(FieldLabel);
      PreNumbersOfFields.push_back(1);
    }
  }
  // end: how many different fields?
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  //cout << "PreDifferentFields\n" << PreDifferentFields << endl;
  //cout << "PreNumbersOfFields\n" << PreNumbersOfFields << endl;
  
  const vector<CSector> &Spectrum = Orbifold.GetSectors();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: find point group invariant couplings
  vector<vector<unsigned> > PG_Allowed;
  {
    unsigned Sum_k = 0;
    unsigned Sum_l = 0;

    vector<unsigned> currentNumber(Order, 0);
    vector<unsigned> MaxDigits(Order, Spectrum.size());
    vector<unsigned> tmp_vec_adjoining_Positions;

    // insert "PreNumbersOfFields[i]"-times the element "i"
    for (i = 0; i < PreDifferentFields.size(); ++i)
      tmp_vec_adjoining_Positions.insert(tmp_vec_adjoining_Positions.end(), PreNumbersOfFields[i], i);

    vector<vector<unsigned> > vec_adjoining_Positions(1, tmp_vec_adjoining_Positions);
    vector<vector<unsigned> > SectorCouplings;

    RecursiveCounting(currentNumber, 0, MaxDigits[0], MaxDigits, vec_adjoining_Positions, SectorCouplings);
    s1 = SectorCouplings.size();

    Sum_k = 0;
    Sum_l = 0;

    unsigned tmp_Sector = 0;
    vector<unsigned> FieldIndices;
    
    bool right_chiral_Sector    = false;
    bool allowed_from_FixedTori = false;

    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &SectorCoupling = SectorCouplings[i];

      Sum_k = 0;
      Sum_l = 0;
      right_chiral_Sector = false;

      for (j = 0; j < Order; ++j)
      {
        tmp_Sector = SectorCoupling[j];

        FieldIndices.clear();
        Spectrum[tmp_Sector].GetFieldIndices(Fields, LeftChiral, FieldIndices);
        
        if (FieldIndices.size() == 0)
        {
          right_chiral_Sector = true;
          break;
        }
        Sum_k += Spectrum[tmp_Sector].Get_m();		//hacking here!!! not optimized for ZMxZNxZK
        Sum_l += Spectrum[tmp_Sector].Get_n();
      }

      // if the coupling is allowed from point group selection rule
      if (!right_chiral_Sector && (Sum_k % M == 0) && (Sum_l % N == 0))
      {
        allowed_from_FixedTori = true;

        // begin: check fixed tori
        vector<bool> Tori_Fixed(3, false);

        bool UntwistedSector_Found = false;
        for (j = 0; j < Order; ++j)
        {
          tmp_Sector = SectorCoupling[j];

          if (tmp_Sector == 0)
          {
            UntwistedSector_Found = true;
            break;
          }
          else
          {
            if (j == 0)
              Tori_Fixed = this->YC_WhereAreFixedTori[tmp_Sector];
            else
            {
              for (k = 0; k < 3; ++k)
              {
                if (Tori_Fixed[k] && !this->YC_WhereAreFixedTori[tmp_Sector][k])
                  Tori_Fixed[k] = false;
              }
            }
          }
        }
        if (!UntwistedSector_Found)
        {
          // if the all fields involved in the coupling have a fixed tori in the same plane and
          // no untwistwed field is involved, the coupling is forbidden
          for (j = 0; allowed_from_FixedTori && (j < 3); ++j)
          {
            if (Tori_Fixed[j])
              allowed_from_FixedTori = false;
          }
        }
        // end: check fixed tori
        if (allowed_from_FixedTori)
          PG_Allowed.push_back(SectorCoupling);
      }
    }
  }
  // end: find point group invariant couplings
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //cout << "PG_Allowed\n" << PG_Allowed << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: again, how many different fields?

  // "AllNumbersOfFields[i][j]" countes the number the field
  // "AllDifferentFields[i][j]" appears in the coupling,
  // seperated by the point group invariant coupling "PG_Allowed[i]"
  vector<vector<string> >   AllDifferentFields;
  vector<vector<unsigned> > AllNumbersOfFields;
  vector<vector<unsigned> > AllPGSectors;
  {
    unsigned counter = 0;
    vector<string>   tmp_AllDifferentFields;
    vector<unsigned> tmp_AllNumbersOfFields;
    vector<unsigned> tmp_AllPGSectors;

    for (i = 0; i < PG_Allowed.size(); ++i)
    {
      const vector<unsigned> &Sectors = PG_Allowed[i];

      tmp_AllDifferentFields.clear();
      tmp_AllNumbersOfFields.clear();
      tmp_AllPGSectors.clear();

      l = 0;
      for (j = 0; j < PreDifferentFields.size(); ++j)
      {
        counter = 1;
        if (PreNumbersOfFields[j] != 1)
        {
          for (k = 0; k < PreNumbersOfFields[j]-1; ++k)
          {
            if (Sectors[l] == Sectors[l+1])
              ++counter;
            else
            {
              tmp_AllDifferentFields.push_back(PreDifferentFields[j]);
              tmp_AllNumbersOfFields.push_back(counter);
              tmp_AllPGSectors.push_back(Sectors[l]);
              counter = 1;
            }
            ++l;
          }
        }

        tmp_AllDifferentFields.push_back(PreDifferentFields[j]);
        tmp_AllNumbersOfFields.push_back(counter);
        tmp_AllPGSectors.push_back(Sectors[l]);
        ++l;
      }
      AllDifferentFields.push_back(tmp_AllDifferentFields);
      AllNumbersOfFields.push_back(tmp_AllNumbersOfFields);
      AllPGSectors.push_back(tmp_AllPGSectors);
    }
  }
  // end: again, how many different fields?
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /*cout << "AllPGSectors " << AllPGSectors << endl;
  cout << "AllNumbersOfFields " << AllNumbersOfFields << endl;
  for (unsigned i = 0; i < AllDifferentFields.size(); ++i)
  {
    for (unsigned j = 0; j < AllDifferentFields[i].size(); ++j)
       cout << AllDifferentFields[i][j] << " ";
    cout << endl;
  }*/
    
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: create the sets of fields that may couple
  //        (save the result to "PG_PoolsOfFields")
  vector<vector<vector<unsigned> > > PG_PoolsOfFields;
  {
    vector<string>::const_iterator pos1;

    bool FieldAllowed       = true;
    bool Look4AllowedFields = false;

    if (AllowedFields != NULL)
      Look4AllowedFields = true;

    vector<unsigned> tmp_Pool;
    vector<vector<unsigned> > tmp_PG_PoolsOfFields;

    s1 = AllPGSectors.size();
    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &Sectors         = AllPGSectors[i];
      const vector<string>   &DifferentFields = AllDifferentFields[i];

      tmp_PG_PoolsOfFields.clear();

      s2 = DifferentFields.size();
      for (j = 0; j < s2; ++j)
      {
        string   FieldLabel = DifferentFields[j];
        unsigned Sector     = Sectors[j];

        tmp_Pool.clear();
        // run through all fields
        for (k = 0; k < f1; ++k)
        {
          const CField &tmp_Field = Fields[k];
          
          // is this field is left-chiral and from the correct sector?
          if ((tmp_Field.Multiplet == LeftChiral) && (tmp_Field.GetInternalIndex()[0] == Sector))
          {
            std::ostringstream os;
        
            tmp_Label1 = tmp_Field.Labels[Vacuum.use_Labels];
            os << tmp_Field.Numbers[Vacuum.use_Labels];
        
            tmp_Label2 = tmp_Label1;
            tmp_Label1 += "_";
            tmp_Label1 += os.str();

            if ((FieldLabel == tmp_Label2) || (FieldLabel == tmp_Label1) || (FieldLabel == "*"))
            {
              if (Look4AllowedFields)
              {
                FieldAllowed = true;
    
                pos1 = find(AllowedFields->begin(), AllowedFields->end(), tmp_Label1);
                if (pos1 == AllowedFields->end())
                {
                  pos1 = find(AllowedFields->begin(), AllowedFields->end(), tmp_Label2);
                  if (pos1 == AllowedFields->end())
                    FieldAllowed = false;
                }
                if (FieldAllowed)
                  tmp_Pool.push_back(k);
              }
              else
                tmp_Pool.push_back(k);
            }
          }
        }
        if (tmp_Pool.size() != 0)
          tmp_PG_PoolsOfFields.push_back(tmp_Pool);
      }
      PG_PoolsOfFields.push_back(tmp_PG_PoolsOfFields);
    }
  }
  // end: create the sets of fields that may couple
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: create all combinations of fields
  //        (save the result to the file "tmp_combinations")
  std::ostringstream osID;
  osID << getpid();
  const string procID = "ID" + osID.str();
  string combi_filename = "tmp_combinations_" + procID + ".txt";

  {
    std::ofstream combi;
    combi.open(combi_filename.data(), ofstream::out);

    this->tmp_file = &combi;

    vector<unsigned> MaxDigits;
    vector<unsigned> vec_adjoining_Positions;

    s1 = PG_PoolsOfFields.size();
    for (i = 0; i < s1; ++i)
    {
      MaxDigits.clear();
      vec_adjoining_Positions.clear();

      const vector<vector<unsigned> > &PoolsOfFields           = PG_PoolsOfFields[i];
      const vector<unsigned>          &NumbersOfFieldsFromPool = AllNumbersOfFields[i];
      
      s2 = PoolsOfFields.size();

      // only if no pool is empty
      if (s2 == NumbersOfFieldsFromPool.size())
      {
        // begin: Recursive Counting
        for (j = 0; j < s2; ++j)
        {
          const vector<unsigned> &PoolOfFields     = PoolsOfFields[j];
          const size_t           PoolOfFields_size = PoolOfFields.size();
          const unsigned         NumberOfFields    = NumbersOfFieldsFromPool[j];

          for (k = 0; k < NumberOfFields; ++k)
          {
            MaxDigits.push_back(PoolOfFields_size-1);
            vec_adjoining_Positions.push_back(j);
          }
        }
        this->CreatePreCouplings(Orbifold, Vacuum, MaxDigits, vec_adjoining_Positions, PoolsOfFields, i, NumbersOfFieldsFromPool);
      }
      // end: Recursive Counting
    }
    combi.close();
  }
  // end: create all combinations of fields
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  unsigned ti = 0;
  unsigned tj = 0;
  size_t  ts1 = 0;

  vector<vector<unsigned> > vec_adjoining_Positions;
  vector<unsigned>          tmp_Combination(Order, 0);
  vector<unsigned>          Number_Of_Weights(Order, 0);

  CVector Sum_Weights(16);
  const CVector NullVector16(16);

  vector<vector<bool> >     gauge_invariants;
  vector<bool>              tmp_gauge_invariant;

  vector<intMatrix>         all_Weights_DL;
  vector<vector<double> >   Weights_CW_Copy;
  
  // begin: define some constants and variables
  unsigned counter = 0;

  unsigned tk = 0;
  size_t   ts2 = 0;
  size_t   ts3 = 0;
  size_t   ts4 = 0;

  bool SG_coupling_allowed = false;

  CSpaceGroupElement Product;
  CSpaceGroupElement result;

  vector<vector<unsigned> > AllNumbers;
  vector<vector<double> >   double_Basis_of_Lattice;
  vector<vector<double> >   tmp_Basis_of_Lattice;

  vector<unsigned> Number_Of_FixedBranes(Order, 0);
  vector<unsigned> Number_Of_ConjugacyClasses(Order, 1);
  vector<unsigned> tmp_Coupling(Order, 0);

  vector<CVector> pre_Basis_of_Lattice;
  vector<CVector> new_BasisVectors;
  vector<CVector> Basis_of_Lattice;
  vector<CVector> Basis_of_DualLattice;

  unsigned tmp_Sector      = 0;
  unsigned Dim_DualLattice = 0;
  double   lambda_ti       = 0.0;
  bool     Use_DualBasis   = false;

  const CVector NullVector4(4);
  const CVector NullVector6(6);
  
  CVector new_BasisVector(6);
  CVector Test_Product_e;
  complexVector Product_e_complex;
  CVector Product_e;

  vector<unsigned> Number_Of_States(Order, 0);
  vector<unsigned> Number_Of_Reps(Order, 0);

  vector<vector<vector<unsigned> > > tmp1_AllowedCouplingsCState;
  vector<vector<unsigned> >          tmp2_AllowedCouplingsCState;
  vector<vector<unsigned> >          CouplingsCState;

  YukawaCoupling tmp_YukawaCoupling;
  // end: define some constants and variables

  vector<unsigned> ConjugacyClass(Order, 0);

  vector<unsigned> tmp_PossibleFieldCoupling;
  vector<unsigned> PossibleFieldCoupling(Order, 0);
  unsigned tmp_index = 0;

  unsigned temp = 0;

  srand ( time(NULL) );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: Do the newly created couplings have to be compared to the old ones?
  bool Check_if_known = true;
  {
    const size_t t1 = Vacuum.FieldCouplings.size();
    if (t1 == 0)
      Check_if_known = false;
    else
    {
      bool CurrentOrderNotKnown = true;

      // run through all couplings that have been created before
      for (i = 0; i < t1; ++i)
      {
        if (Vacuum.FieldCouplings[i].FieldIndices.size() == Order)
        {
          CurrentOrderNotKnown = false;
          break;
        }
      }
      if (CurrentOrderNotKnown)
        Check_if_known = false;
    }
  }
  // end: Do the newly created couplings have to be compared to the old ones?
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  const CSpaceGroupElement Identity;
    
  std::ifstream in_combi;
  in_combi.open(combi_filename.data(), ifstream::in);

  //vector<YukawaCoupling> TestFieldCouplings;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check the couplings
  string tmp_string = "";
  while (getline(in_combi, tmp_string) )
  {
    tmp_PossibleFieldCoupling.clear();
    temp = 0;

    std::istringstream line(tmp_string);
    while (line >> temp) tmp_PossibleFieldCoupling.push_back(temp);

    tmp_index = 1;
    for (i = 0; i < PG_PoolsOfFields[tmp_PossibleFieldCoupling[0]].size(); ++i)
    {
      const vector<unsigned> &PoolOfFields     = PG_PoolsOfFields[tmp_PossibleFieldCoupling[0]][i];
      const unsigned         NumberOfFields    = AllNumbersOfFields[tmp_PossibleFieldCoupling[0]][i];

      for (j = 0; j < NumberOfFields; ++j)
      {
        PossibleFieldCoupling[tmp_index-1] = PoolOfFields[tmp_PossibleFieldCoupling[tmp_index]];
        ++tmp_index;
      }
    }

    // space group selection rule

    // if the coupling from the twisted sectors, specified in "CouplingCSector", is known,
    // do not need to create the "Basis_of_DualLattice" and the "Basis_of_Lattice"

    vector<unsigned> CouplingCSector(Order,0);
    for (ti = 0; ti < Order; ++ti)
      CouplingCSector[ti] = Fields[PossibleFieldCoupling[ti]].GetInternalIndex()[0];

    stable_sort(CouplingCSector.begin(), CouplingCSector.end());
    
    vector<vector<unsigned> >::iterator pos = find(this->KnownCouplingsCSector.begin(), this->KnownCouplingsCSector.end(), CouplingCSector);
    
    if (pos != this->KnownCouplingsCSector.end())
    {
      const unsigned pos_index = distance(this->KnownCouplingsCSector.begin(), pos);

      Basis_of_DualLattice.clear();
      Basis_of_Lattice.clear();

      Basis_of_DualLattice = this->KnownBasis_of_DualLattice[pos_index];
      Basis_of_Lattice     = this->KnownBasis_of_Lattice[pos_index];
      Dim_DualLattice      = this->KnownDim_DualLattice[pos_index];
      Use_DualBasis        = this->KnownUse_DualBasis[pos_index];
    }
    else
    {
      this->KnownCouplingsCSector.push_back(CouplingCSector);

      // begin: create a pre basis of the lattice T^6
      Basis_of_DualLattice.clear();
      Basis_of_Lattice.clear();
      Dim_DualLattice = 0;
      Use_DualBasis   = false;

      pre_Basis_of_Lattice.clear();

      bool Sector_known      = false;
      bool BasisVector_known = false;

      for (ti = 0; ti < Order; ++ti)
      {
        tmp_Sector = Fields[PossibleFieldCoupling[ti]].GetInternalIndex()[0];

        Sector_known = false;
        for (tj = 0; !Sector_known && (tj < ti); ++tj)
        {
          if (tmp_Sector == Fields[PossibleFieldCoupling[tj]].GetInternalIndex()[0])
            Sector_known = true;
        }

        if (!Sector_known)
        {
          const CSector          &Sector    = Spectrum[tmp_Sector];
          const vector<CVector>  &tmp_Basis = this->YC_Basis_of_Rotated_Lattices_For_kl_Twisted_Sectors[Sector.Get_m()][Sector.Get_n()];		//hacking here!!! not optimized for ZMxZNxZK

          ts1 = tmp_Basis.size();
          for (tj = 0; tj < ts1; ++tj)
          {
            const CVector &tmp_BasisVector = tmp_Basis[tj];
            BasisVector_known = false;

            ts2 = pre_Basis_of_Lattice.size();
            for (tk = 0; !BasisVector_known && (tk < ts2); ++tk)
            {
              const CVector &pre_BasisVector = pre_Basis_of_Lattice[tk];

              if ((tmp_BasisVector == pre_BasisVector) || (tmp_BasisVector * (-1) == pre_BasisVector))
                BasisVector_known = true;
            }
            if (!BasisVector_known)
              pre_Basis_of_Lattice.push_back(tmp_BasisVector);
          }
        }
      }
      ts1 = pre_Basis_of_Lattice.size();

      // begin: create new pre-basis vectors by linear combination
      if (ts1 != 0)
      {
        // begin: find the maximal length
        stable_sort(pre_Basis_of_Lattice.begin(), pre_Basis_of_Lattice.end(), length_criterion);
        const CVector &longest_BasisVector = pre_Basis_of_Lattice[ts1-1];
        const double  length               = longest_BasisVector.GetSqrTo(6);
        // end: find the maximal length

        vector<unsigned> currentNumber(ts1,0);
        vector<unsigned> MaxDigits(ts1,3);
        AllNumbers.clear();

        RecursiveCounting(currentNumber, 0, MaxDigits[0], MaxDigits, AllNumbers);

        new_BasisVectors.clear();

        ts2 = AllNumbers.size();
        for (ti = 0; ti < ts2; ++ti)
        {
          new_BasisVector = NullVector6;

          const vector<unsigned> &Numbers = AllNumbers[ti];
          for (tj = 0; tj < ts1; ++tj)
          {
            const unsigned Number = Numbers[tj];
            if (Number == 0)
              new_BasisVector = new_BasisVector - pre_Basis_of_Lattice[tj];
            else
            if (Number == 2)
              new_BasisVector = new_BasisVector + pre_Basis_of_Lattice[tj];
          }

          BasisVector_known = false;
          // first check the length
          if (new_BasisVector.GetSqrTo(6) > length)
            BasisVector_known = true;

          ts3 = new_BasisVectors.size();
          // then check whether the new pre-basis vector is known
          for (tj = 0; !BasisVector_known && (tj < ts3); ++tj)
          {
            const CVector &pre_BasisVector = new_BasisVectors[tj];
            if ((new_BasisVector == pre_BasisVector) || (new_BasisVector * (-1) == pre_BasisVector))
              BasisVector_known = true;
          }

          for (tj = 0; !BasisVector_known && (tj < ts1); ++tj)
          {
            const CVector &pre_BasisVector = pre_Basis_of_Lattice[tj];
            if ((new_BasisVector == pre_BasisVector) || (new_BasisVector * (-1) == pre_BasisVector))
              BasisVector_known = true;
          }

          if (!BasisVector_known)
            new_BasisVectors.push_back(new_BasisVector);
        }

        pre_Basis_of_Lattice.insert(pre_Basis_of_Lattice.end(), new_BasisVectors.begin(), new_BasisVectors.end());
      }
      // end: create new pre-basis vectors by linear combination
  
      stable_sort(pre_Basis_of_Lattice.begin(), pre_Basis_of_Lattice.end(), length_criterion);
      // end: create a pre basis of the lattice T^6

      // begin: create a basis of the lattice T^6
      double_Basis_of_Lattice.clear();
      tmp_Basis_of_Lattice.clear();

      ts2 = pre_Basis_of_Lattice.size();
      bool continue_loop = true;
      
      for (ti = 0; continue_loop && (ti < ts2); ++ti)
      {
        const CVector        &tmp_Basis   = pre_Basis_of_Lattice[ti];
        const vector<double> double_Basis = tmp_Basis;

        if (double_Basis_of_Lattice.size() != 0)
          tmp_Basis_of_Lattice = findBasis<double>(double_Basis_of_Lattice);

        ts3 = tmp_Basis_of_Lattice.size();

        double_Basis_of_Lattice.push_back(double_Basis);

        tmp_Basis_of_Lattice = findBasis<double>(double_Basis_of_Lattice);
        ts4 = tmp_Basis_of_Lattice.size();

        if (ts4 == ts3 + 1)
          Basis_of_Lattice.push_back(tmp_Basis);
        else
          double_Basis_of_Lattice.pop_back();

        if (Basis_of_Lattice.size() == 6)
          continue_loop = false;
      }
      // end: create a basis of the lattice T^6

      // begin: create a basis of the D U A L  L A T T I C E
      Use_DualBasis   = false;
      Dim_DualLattice = Basis_of_Lattice.size();
      Basis_of_DualLattice.resize(Dim_DualLattice);

      if (Dim_DualLattice != 0)
      {
        Use_DualBasis = true;

        // begin: compute the metric of the rotated six-torus

        // Metric[l][m] = Basis_of_Lattice[l] * Basis_of_Lattice[m]
        intMatrix Metric(Dim_DualLattice);

        for (unsigned l = 0; l < Dim_DualLattice; ++l)
        {
          const CVector &Basis_Vector = Basis_of_Lattice[l];

          intVector row(Dim_DualLattice,0);
          for (unsigned m = 0; m < Dim_DualLattice; ++m)
          {
            const double tmp1 = Basis_Vector * Basis_of_Lattice[m];
            const double tmp2 = round_double_to_int(tmp1);

            if (fabs(tmp1 - tmp2) > prec)
            {
              cout << "\n  Warning in unsigned CYukawaCouplings::AddCoupling(...) : Can not compute the metric of the rotated lattice. Return 0." << endl;
              return 0;
            }
            row[m] = (int)tmp2;
          }
          Metric[l] = row;
        }
        // end: compute the metric of the rotated six-torus

        rationalMatrix Inverse_Metric = inverse(Metric);

        // begin: compute the dual basis of the rotated six-torus
        for (unsigned l = 0; l < Dim_DualLattice; ++l)
        {
          const rationalVector &rat_row = Inverse_Metric[l];

          CVector BasisVector(6);
          for (unsigned m = 0; m < Dim_DualLattice; ++m)
          {
            const rational<int> &rat_g_lm = rat_row[m];
            const double        g_lm      = ((double)rat_g_lm.numerator())/((double)rat_g_lm.denominator());
            const CVector       &e_m      = Basis_of_Lattice[m];
            BasisVector = BasisVector + (g_lm * e_m);
          }
          Basis_of_DualLattice[l] = BasisVector;
        }
        // end: compute the dual basis on the rotated six-torus

        // begin: ERROR check
        doubleMatrix tmp_Basis_of_Lattice(Dim_DualLattice);
        doubleMatrix tmp_Basis_of_DualLattice(Dim_DualLattice);

        for (unsigned l = 0; l < Dim_DualLattice; ++l)
        {
          tmp_Basis_of_Lattice[l]     = Basis_of_Lattice[l];
          tmp_Basis_of_DualLattice[l] = Basis_of_DualLattice[l];
        }
        doubleMatrix unity = tmp_Basis_of_DualLattice * transpose(tmp_Basis_of_Lattice);

        bool error = false;

        const size_t us = unity.size();

        if (us == 0)
          error = true;

        for (unsigned l = 0; !error && (l < us); ++l)
        {
          const doubleVector &tmpVector = unity[l];

          const size_t vs = tmpVector.size();

          if (vs == 0)
            error = true;

          for (unsigned m = 0; !error && (m < vs); ++m)
          {
            const double tmp = tmpVector[m];
            if (l == m)
            {
              if (fabs(tmp - 1.0) > prec)
                error = true;
            }
            else
            {
              if (fabs(tmp) > prec)
                error = true;
            }
          }
        }
        if (error)
        {
          cout << "\n  Warning in unsigned CYukawaCouplings::AddCoupling(...) : Check the dual basis. Return 0." << endl;
          return 0;
        }

        // check whether the vectors of Basis_of_Lattice are the smallest basis of "pre_Basis_of_Lattice"
        for (unsigned l = 0; l < ts1; ++l)
        {
          const CVector &pre_BasisVector = pre_Basis_of_Lattice[l];

          for (unsigned m = 0; m < Dim_DualLattice; ++m)
          {
            const double lambda = pre_BasisVector * Basis_of_DualLattice[m];

            // is lambda integer?
            if (!is_integer(lambda))
            {
              cout << "\n  Warning in unsigned CYukawaCouplings::AddCoupling(...) : Check the basis of the lattice T^6. Return 0." << endl;
              return 0;
            }
          }
        }
        // end: ERROR check
      }

      this->KnownBasis_of_DualLattice.push_back(Basis_of_DualLattice);
      this->KnownBasis_of_Lattice.push_back(Basis_of_Lattice);
      this->KnownDim_DualLattice.push_back(Dim_DualLattice);
      this->KnownUse_DualBasis.push_back(Use_DualBasis);
    }
    // end: create a basis of the D U A L  L A T T I C E

    // begin: run through the conjugacy classes
    Number_Of_ConjugacyClasses.assign(Order, 1);
    for (l = 0; l < Order; ++l)
    {
      const vector<unsigned> &internalIndex = Fields[PossibleFieldCoupling[l]].GetInternalIndex();
      Number_Of_ConjugacyClasses[l] = this->YC_CCs_of_constructing_Elements[internalIndex[0]][internalIndex[1]].size();
    }

    ConjugacyClass.assign(Order, 0);

    // assume that the space group selection rule forbids the coupling
    SG_coupling_allowed = false;
    do
    {
      Product = Identity;

      for (l = 0; l < Order; ++l)
      {
        const vector<unsigned> &internalIndex = Fields[PossibleFieldCoupling[l]].GetInternalIndex();

        SpaceGroup.SG_Multiply(Product, this->YC_CCs_of_constructing_Elements[internalIndex[0]][internalIndex[1]][ConjugacyClass[l]], result);
        Product = result;
      }

      // begin: check again the P O I N T  G R O U P selection rule
      if (!Product.NoTwist())
      {
        cout << "\n  Warning in void CYukawaCouplings::AddCoupling(...): First check of the Point Group selection rule was ok, second check was not. Return 0." << endl;
        return 0;
      }
      // end: check again the P O I N T  G R O U P selection rule

      if (Product.IsZero())
        SG_coupling_allowed = true;
      else
      {
        Product_e.clear();
        SpaceGroup.SG_ReverseVector(Product, Product_e_complex);
        
        for (l = 0; l < 3; ++l)
        {
          Product_e.Push_back(Product_e_complex[l].real());
          Product_e.Push_back(Product_e_complex[l].imag());
        }
        
        lambda_ti = 0.0;
        
        // begin: S P A C E  G R O U P  selection rule
        if (Product_e == NullVector6)
          SG_coupling_allowed = true;
        else
        {
          if (Use_DualBasis)
          {
            Test_Product_e = NullVector6;
  
            // check whether the shift "Product_e" of product of the constructing elements is in the rotated six-torus
            // by: Product_e * Basis_of_DualLattice[ti] = lambda_ti integer ?
  
            SG_coupling_allowed = true;
            for (ti = 0; SG_coupling_allowed && (ti < Dim_DualLattice); ++ti)
            {
              lambda_ti = Product_e * Basis_of_DualLattice[ti];
  
              // is lambda_ti integer?
              if (!is_integer(lambda_ti))
                SG_coupling_allowed = false;
  
              Test_Product_e = Test_Product_e + (Basis_of_Lattice[ti] * lambda_ti);
            }
            if (SG_coupling_allowed && (Test_Product_e != Product_e))
            {
              cout << "\n  Warning in unsigned CYukawaCouplings::AddCoupling(...) : Bsis of the dual lattice is not complete. Return 0." << endl;
              return 0;
            }
          }
          else
            SG_coupling_allowed = false;
        }
        // end: S P A C E  G R O U P  selection rule
      }
    } while (!SG_coupling_allowed && NextNumber(ConjugacyClass, Number_Of_ConjugacyClasses, Order));
    // end: run through the conjugacy classes

    if (SG_coupling_allowed)
    {
      stable_sort(PossibleFieldCoupling.begin(), PossibleFieldCoupling.end());

      bool Coupling_known = false;

      if (Check_if_known)
      {
        const size_t t1 = Vacuum.FieldCouplings.size();
        for (j = 0; !Coupling_known && (j < t1); ++j)
        {
          if (Vacuum.FieldCouplings[j].FieldIndices == PossibleFieldCoupling)
            Coupling_known = true;
        }
      }
      if (!Coupling_known)
      {
        ++counter;

        tmp_YukawaCoupling.LabelOfModuli.clear();
        tmp_YukawaCoupling.ExponentsOfEtas.clear();

        //bool one_is_negative = false;
        for (j = 0; j < number_of_modularsymmetries; ++j)
        {
          exponent = 1;
          for (k = 0; k < Order; ++k)
            exponent += Fields[PossibleFieldCoupling[k]].GetModularWeight(ModularSymmetries[j], Spectrum);

          exponent *= -2;

          tmp_YukawaCoupling.LabelOfModuli.push_back(ModularSymmetries[j].Label);
          tmp_YukawaCoupling.ExponentsOfEtas.push_back(exponent.numerator());
          
          //if (exponent < 0)
          //  one_is_negative = true;
        }

        tmp_YukawaCoupling.FieldIndices     = PossibleFieldCoupling;
        tmp_YukawaCoupling.CouplingStrength = -1.276453 * ((double)(rand() + 1))/((double)(RAND_MAX + 1));

        //if (one_is_negative)
        //  TestFieldCouplings.push_back(tmp_YukawaCoupling);

        Vacuum.FieldCouplings.push_back(tmp_YukawaCoupling);
        /*
        for (j = 0; j < Order; ++j)
        {
          const CField &tmp_Field = Fields[PossibleFieldCoupling[j]];
          cout << tmp_Field.Label << "_" << tmp_Field.Number << " ";
        }
        cout << endl;
        */
      }
    }
      
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  in_combi.close();

  /*if (TestFieldCouplings.size() != 0)
  {
    cout << "couplings with strange moduli dependence:" << endl;
    CPrint Print(Tstandard, &cout);
    Print.PrintCouplings(Vacuum, TestFieldCouplings, false, true);
    cout << endl;
  }*/

  string filename = "tmp_file";
  filename += procID;
  filename += ".txt";

  std::ofstream tmp_out(filename.data());

  this->Save(Orbifold, Vacuum, tmp_out, 0);
  tmp_out.close();

  return counter;
}



/* ##########################################################################
######   SC_criterion(...)                                             ######
######                                                                 ######
######   Version: 09.03.2011                                           ######
########################################################################## */
bool SC_criterion(const vector<RepVector> &SP1, const vector<RepVector> &SP2)
{
  if (SP1.size() > SP2.size())
    return true;
  
  if (SP1.size() < SP2.size())
    return false;
  
  const RepVector &Rep1 = SP1[0];
  const RepVector &Rep2 = SP2[0];
  
  const size_t s1 = Rep1.size();
  if (s1 != Rep2.size())
  {
    cout << "\n  Warning in bool SC_criterion(...): Return false." << endl;
    return false;
  }
  for (unsigned i = 0; i < s1; ++i)
  {
    if (Rep1[i].Dimension > Rep2[i].Dimension)
      return true;
    if (Rep1[i].Dimension < Rep2[i].Dimension)
      return false;
  }
  
  return false;
}



/* ########################################################################################
######   bool CYukawaCouplings::DeleteCouplingsWithStatesTwice(...)                  ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CYukawaCouplings::DeleteCouplingsWithStatesTwice(SConfig &Vacuum, const CPrint &Print, vector<vector<string> > &Known_GG_string, vector<vector<vector<vector<RepVector> > > > &Known_Zero_SortedCoupling, vector<vector<vector<vector<RepVector> > > > &Known_NonZero_SortedCoupling)
{
  const vector<CField>                    &Fields            = Vacuum.Fields;
  const vector<gaugeGroupFactor<double> > &GaugeGroupFactors = Vacuum.SymmetryGroup.GaugeGroup.factor;

  const unsigned NumberOfGGfs = GaugeGroupFactors.size();

  bool is_coupling_zero = false;

  SDimension tmp_One;
  tmp_One.Dimension = 1;
  tmp_One.AdditionalLabel = "";
  const RepVector Singlet(NumberOfGGfs, tmp_One);

  vector<RepVector>          SubPartCoupling;
  vector<vector<RepVector> > SortedCoupling;
  vector<unsigned>           SC_FieldIndices;

  unsigned order = 0;

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;

  size_t q1 = 0;
  size_t s1 = 0;

  vector<string> GG_string;
  const vector<vector<vector<RepVector> > > EmptySet;
  vector<unsigned>::iterator pos;
  unsigned Index = 0;

  // run through all couplings
      for( vector<YukawaCoupling>::iterator CurrentCoupling = Vacuum.FieldCouplings.begin(); CurrentCoupling != Vacuum.FieldCouplings.end(); ++CurrentCoupling)
  {
    is_coupling_zero = false;

    const vector<unsigned> &FieldIndices = CurrentCoupling->FieldIndices;

    order = FieldIndices.size();
    if (order == 0)
    {
      cout << "\n  Warning in bool CYukawaCouplings::DeleteCouplingsWithStatesTwice(...): Order of coupling is zero. Return false." << endl;
      return false;
    }

    SC_FieldIndices.clear();
    SortedCoupling.clear();

    // begin: sort the current coupling
    // run through the non-singlet fields involded in the current coupling
    for (i = 0; i < order; ++i)
    {
      Index = FieldIndices[i];

      const RepVector &Dimensions = Fields[Index].Dimensions;

      // only non-singlets
      if (!AreRepVectorsEqual(Dimensions, Singlet))
      {
        // begin: count how often each field appears in the current coupling
        pos = find(SC_FieldIndices.begin(), SC_FieldIndices.end(), Index);
        if (pos == SC_FieldIndices.end())
        {
          SC_FieldIndices.push_back(Index);

          SubPartCoupling.clear();
          SubPartCoupling.push_back(Dimensions);
          SortedCoupling.push_back(SubPartCoupling);
        }
        else
          SortedCoupling[distance(SC_FieldIndices.begin(), pos)].push_back(Dimensions);
        // end: count how often each field appears in the current coupling
      }
    }
    // end: sort the current coupling

    // now, SC_FieldIndices is a set of field indices of non-singlet fields involved in the current coupling
    // e.g. V_1 has fields index 5
    //      V_2 has fields index 23
    //      V_3 has fields index 42 (and all V's are (4,1)-plets of SU(4)xSU(2))
    //      and assume the coupling V1^2 V_2 V_3
    //      then SC_FieldIndices = (5,23,42)
    // furthermore, SortedCoupling contains the coupling in terms of the representations
    // e.g. SortedCoupling = {{(4,1),(4,1)}, {(4,1)}, {(4,1)}}

    // if one field appears more than once, the coupling may be zero
    s1 = SortedCoupling.size();
    for (i = 0; !is_coupling_zero && (i < s1); ++i)
    {
      if (SortedCoupling[i].size() > 1)
        is_coupling_zero = true;
    }

    // coupling may be zero because some fields with non-trivial representation appear more than once
    if (is_coupling_zero)
    {
      vector<unsigned> UsedParts(1,i);
      vector<bool> ChargedUnderGG(NumberOfGGfs, false);

      vector<vector<RepVector> > OptimizedSortedCoupling;

      is_coupling_zero = false;
      // run through the "dangerous" parts of the couplings
      for (i = 0; !is_coupling_zero && (i < s1); ++i)
      {
        const vector<RepVector> &SubPart = SortedCoupling[i];

        // only if this SubPart is "dangerous"
        if (SubPart.size() > 1)
        {
          UsedParts.assign(1,i);
          ChargedUnderGG.assign(NumberOfGGfs, false);

          for (j = 0; j < NumberOfGGfs; ++j)
          {
            if (SubPart[0][j].Dimension != 1)
              ChargedUnderGG[j] = true;
          }
          // begin: search for representations which are "linked" to the one of "SubPart"
          bool new_link_found = true;
          while (new_link_found)
          {
            new_link_found = false;

            for (j = 0; j < SortedCoupling.size(); ++j)
            {
              if (find(UsedParts.begin(), UsedParts.end(), j) == UsedParts.end())
              {
                const RepVector &Compare2SubPart = SortedCoupling[j][0];

                bool charged = false;
                for (k = 0; k < NumberOfGGfs; ++k)
                {
                  if (ChargedUnderGG[k] && (Compare2SubPart[k].Dimension != 1))
                    charged = true;
                }
                if (charged)
                {
                  for (k = 0; k < NumberOfGGfs; ++k)
                  {
                    if (Compare2SubPart[k].Dimension != 1)
                    {
                      UsedParts.push_back(j);
                      ChargedUnderGG[k] = true;
                      new_link_found = true;
                    }
                  }
                }
              }
            }
          }
          // end: search for representations which are "linked" to the one of "SubPart"
          OptimizedSortedCoupling.clear();

          RepVector         NonTrivialRep;
          vector<RepVector> NonTrivialPart;
          bool charged = false;

          // run through all parts of "SortedCoupling"
          for (j = 0; j < SortedCoupling.size(); ++j)
          {
            const vector<RepVector> &Part = SortedCoupling[j];
            charged = false;

            // begin: is "Part" charged?
            for (k = 0; k < NumberOfGGfs; ++k)
            {
              if (ChargedUnderGG[k] && (Part[0][k].Dimension != 1))
              {
                charged = true;
                break;
              }
            }
            // end: is "Part" charged?

            if (charged)
            {
              // begin: get the non-singlet part of "Part"
              NonTrivialPart.clear();
              for (k = 0; k < Part.size(); ++k)
              {
                const RepVector &Dimensions = Part[k];
                NonTrivialRep.clear();
                for (l = 0; l < NumberOfGGfs; ++l)
                {
                  if (ChargedUnderGG[l])
                    NonTrivialRep.push_back(Dimensions[l]);
                }
                NonTrivialPart.push_back(NonTrivialRep);
              }
              // end: get the non-singlet part of "Part"
              OptimizedSortedCoupling.push_back(NonTrivialPart);
            }
          }
          stable_sort(OptimizedSortedCoupling.begin(), OptimizedSortedCoupling.end(), SC_criterion);

          GG_string.clear();
          for (j = 0; j < NumberOfGGfs; ++j)
          {
            if (ChargedUnderGG[j])
              GG_string.push_back(GaugeGroupFactors[j].algebra);
          }

          vector<vector<string> >::iterator pos = find(Known_GG_string.begin(), Known_GG_string.end(), GG_string);
          if (pos == Known_GG_string.end())
          {
            Known_GG_string.push_back(GG_string);
            Known_Zero_SortedCoupling.push_back(EmptySet);
            Known_NonZero_SortedCoupling.push_back(EmptySet);
            q1 = Known_GG_string.size()-1;
          }
          else
            q1 = distance(Known_GG_string.begin(), pos);

          vector<vector<vector<RepVector> > > &Zero_SortedCoupling    = Known_Zero_SortedCoupling[q1];
          vector<vector<vector<RepVector> > > &NonZero_SortedCoupling = Known_NonZero_SortedCoupling[q1];

          // is the coupling known to be zero?
          if (find(Zero_SortedCoupling.begin(), Zero_SortedCoupling.end(), OptimizedSortedCoupling) != Zero_SortedCoupling.end())
          {
            is_coupling_zero = true;
            break;
          }
          else
          {
            // is the coupling known to be non-zero?
            if (find(NonZero_SortedCoupling.begin(), NonZero_SortedCoupling.end(), OptimizedSortedCoupling) != NonZero_SortedCoupling.end())
              is_coupling_zero = false;
            else
            {
              // otherwise, user has to decide whether coupling vanishes or not
              // begin: print coupling
              (*Print.out) << "\nnew kind of coupling: ";
              for (j = 0; j < order; ++j)
              {
                Print.PrintLabel(Fields[FieldIndices[j]], Vacuum.use_Labels);
                (*Print.out) << " ";
              }
              (*Print.out) << "\nrelevant part of the coupling: ";
              for (j = 0; j < OptimizedSortedCoupling.size(); ++j)
              {
                for (k = 0; k < OptimizedSortedCoupling[j].size(); ++k)
                  Print.PrintRep(OptimizedSortedCoupling[j][k]);
                (*Print.out) << " | ";
              }
              (*Print.out) << endl;
              // end: print coupling
              (*Print.out) << "Coupling zero? (y/n)" << endl;

              string tmp;
              cin >> tmp;
              if (tmp == "y")
              {
                Zero_SortedCoupling.push_back(OptimizedSortedCoupling);
                is_coupling_zero = true;
                (*Print.out) << "zero\n" << endl;
                break;
              }
              else
              {
                NonZero_SortedCoupling.push_back(OptimizedSortedCoupling);
                is_coupling_zero = false;
                (*Print.out) << "nonzero\n" << endl;
              }
            }
          }
        }
      }
      // begin: delete zero coupling
      if (is_coupling_zero)
      {
        (*Print.out) << "delete ";
        for (j = 0; j < order; ++j)
        {
          Print.PrintLabel(Fields[FieldIndices[j]], Vacuum.use_Labels);
          (*Print.out) << " ";
        }
        (*Print.out) << endl;

        Vacuum.FieldCouplings.erase(CurrentCoupling--);
      }
    }
    // end: delete zero coupling
  }
  return true;
}



/* ########################################################################################
######   Initiate(const CSpaceGroup &SpaceGroup, const SConfig &Vacuum)              ######
######                                                                               ######
######   Version: 24.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
bool CYukawaCouplings::Initiate(const CSpaceGroup &SpaceGroup, const SConfig &Vacuum)
{
  const vector<CTwistVector> &Twists = SpaceGroup.GetTwists();
  
  const unsigned      M          = SpaceGroup.GetM();
  const unsigned      N          = SpaceGroup.GetN();
  const bool          ZMxZN      = SpaceGroup.IsZMxZN();
  const complexMatrix T6_Lattice = SpaceGroup.GetT6Lattice();

  CSpaceGroupElement ZM_Twist = SpaceGroup.GetSG_Generators_Twist()[0];
  CSpaceGroupElement ZN_Twist;
  if (ZMxZN)
    ZN_Twist = SpaceGroup.GetSG_Generators_Twist()[1];

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  unsigned m = 0;
  unsigned n = 0;
  unsigned p = 0;

  const vector<vector<CSpaceGroupElement> > &SectorsOfConstructingElements = SpaceGroup.GetSectors(); 
  const size_t s1 = SectorsOfConstructingElements.size();
  size_t s2 = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: create the information about fixed tori
  {
    CVector localTwist(4);

    // WhereIsFixedTorus[i] = false means: i-th torus is not a fixed torus
    // WhereIsFixedTorus[i] = true  means: i-th torus is a fixed torus
    vector<bool> WhereIsFixedTorus(ComplexLatticeDim, false);

    // run through the untwisted and all twisted sectors
    for (i = 0; i < s1; ++i)
    {
      const vector<CSpaceGroupElement> &Elements = SectorsOfConstructingElements[i];

      if (Elements.size() == 0)
      {
        cout << "\n  Warning in bool CYukawaCouplings::Initiate(...) : The sector has no fixed branes. Return false." << endl;
        return false;
      }
      const CSpaceGroupElement &ConstructingElement = Elements[0];
      
      localTwist = Twists[0] * ConstructingElement.Get_m();
      if (ZMxZN)
        localTwist += Twists[1] * ConstructingElement.Get_n();	//hacking here!!! not optimized for ZMxZNxZK

      for (j = 0; j < ComplexLatticeDim; ++j)
      {
        if (is_integer(localTwist[j + 1]))
          WhereIsFixedTorus[j] = true;
        else
          WhereIsFixedTorus[j] = false;
      }

      this->YC_WhereAreFixedTori.push_back(WhereIsFixedTorus);
    }
  }
  // end: create the information about fixed tori
  ////////////////////////////////////////////////////////////////////////////////////////////////////////


  const CSpaceGroupElement Identity;
  CSpaceGroupElement result;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: create "CSpaceGroupElement" of all twists for conjugation of the constructing element
  vector<vector<CSpaceGroupElement> > ZMxZN_Twists;
  vector<vector<CSpaceGroupElement> > ZMxZN_inverse_Twists;
  {
    CSpaceGroupElement Twist;
    vector<CSpaceGroupElement> vector_of_Twists;
    vector<CSpaceGroupElement> vector_of_inverse_Twists;

    for (k = 0; k < M; ++k)
    {
      vector_of_Twists.clear();
      vector_of_inverse_Twists.clear();
      for (l = 0; l < N; ++l)
      {
        Twist = Identity;

        if (k != 0)
        {
          for (i = 0; i < k; ++i)
          {
            SpaceGroup.SG_Multiply(Twist, ZM_Twist, result);
            Twist = result;
          }
        }
        if (l != 0)
        {
          for (i = 0; i < l; ++i)
          {
            SpaceGroup.SG_Multiply(Twist, ZN_Twist, result);
            Twist = result;
          }
        }
        vector_of_Twists.push_back(Twist);
        vector_of_inverse_Twists.push_back(SpaceGroup.SG_Inverse(Twist));
      }
      ZMxZN_Twists.push_back(vector_of_Twists);
      ZMxZN_inverse_Twists.push_back(vector_of_inverse_Twists);
    }
  }
  this->YC_ZMxZN_Twists_size = M * N;
  // end: create "CSpaceGroupElement" of all twists for conjugation of the constructing element
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: create the lattices (1 - \theta^k\theta^l) e_i
  //        Basis_of_Rotated_Lattices_For_kl_Twisted_Sectors[k][l] contains
  //        the lattice for the (k,l) twisted sector;
  {
    vector<vector<CVector> > Twisted_Sectors;
    vector<CVector> Basis_of_Rotated_Lattices;

    complexVector e_i;
    complexVector e_i_prime(ComplexLatticeDim, complex<double>(0,0));
    CVector NullVector6(LatticeDim);
    CVector tmp(LatticeDim);

    for (k = 0; k < M; ++k)
    {
      const vector<CSpaceGroupElement> &k_Twists = ZMxZN_Twists[k];

      Twisted_Sectors.clear();
      for (l = 0; l < N; ++l)
      {
        const CSpaceGroupElement &kl_Twist = k_Twists[l];

        Basis_of_Rotated_Lattices.clear();
        for (i = 0; i < LatticeDim; ++i)
        {
          e_i = T6_Lattice[i];
          SpaceGroup.SG_Multiply(kl_Twist, e_i, e_i_prime);

          for (j = 0; j < ComplexLatticeDim; ++j)
            e_i[j] -= e_i_prime[j];

          tmp = Convert_complexVector2CVector(e_i);

          if (tmp != NullVector6)
            Basis_of_Rotated_Lattices.push_back(tmp);
        }

        Twisted_Sectors.push_back(Basis_of_Rotated_Lattices);
      }
      this->YC_Basis_of_Rotated_Lattices_For_kl_Twisted_Sectors.push_back(Twisted_Sectors);
    }
  }
  // end: create the lattices (1 - \theta^k\theta^l) e_i
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // checked by hand for Z3
  const vector<CSpaceGroupElement> &SG_AllNonStandardShifts = SpaceGroup.GetSG_AllNonStandardShifts();
  size_t s3 = SG_AllNonStandardShifts.size();
  if (s3 == 0)
  {
    cout << "\n  Warning in bool CYukawaCouplings::Initiate(...) : \"SG_AllNonStandardShifts\" is empty. Return false." << endl;
    return false;
  }
  vector<CSpaceGroupElement> SG_All_inverse_NonStandardShifts(s3);
  for (m = 0; m < s3; ++m)
    SG_All_inverse_NonStandardShifts[m] = SpaceGroup.SG_Inverse(SG_AllNonStandardShifts[m]);
  
  vector<vector<CSpaceGroupElement> > tmp2_CCs_of_constructing_Elements;
  vector<CSpaceGroupElement>          tmp1_CCs_of_constructing_Elements;
  CSpaceGroupElement                  shifted_constr_Element;
  CSpaceGroupElement                  rotated_constr_Element, result2;

  bool                  ConstructingElementHasFP = true;
  bool                  AddFixedPoint = true;
  vector<FPCoordinates> FixedPoints;
  FPCoordinates         FixedPoint;
  FPCoordinates         TransformedFixedPoint;
  
  size_t s4 = 0;
  
  // run through all sectors of constructing elements
  for (i = 0; i < s1; ++i)
  {
    tmp2_CCs_of_constructing_Elements.clear();

    const vector<CSpaceGroupElement> &SectorOfConstructingElements = SectorsOfConstructingElements[i];
    s2 = SectorOfConstructingElements.size();

    // run through all constructing elements of the i-th sector
    for (j = 0; j < s2; ++j)
    {
      tmp1_CCs_of_constructing_Elements.clear();

      const CSpaceGroupElement &ConstructingElement = SectorOfConstructingElements[j];
      
      FixedPoints.clear();
      ConstructingElementHasFP = SpaceGroup.SG_FixedPoint(ConstructingElement, FixedPoint);
      if (ConstructingElementHasFP)
        FixedPoints.push_back(FixedPoint);

      tmp1_CCs_of_constructing_Elements.push_back(ConstructingElement);

      // shift using non-standard lattice vectors
      for (m = 0; m < s3; ++m)
      {
        SpaceGroup.SG_Multiply(SG_AllNonStandardShifts[m], ConstructingElement, result);
        SpaceGroup.SG_Multiply(result, SG_All_inverse_NonStandardShifts[m], shifted_constr_Element);

        // rotate using \theta^k
        for (k = 0; k < M; ++k)
        {
          const vector<CSpaceGroupElement> &k_Twists         = ZMxZN_Twists[k];
          const vector<CSpaceGroupElement> &k_inverse_Twists = ZMxZN_inverse_Twists[k];
          
          // rotate using \omega^l
          for (l = 0; l < N; ++l)
          {
            SpaceGroup.SG_Multiply(k_Twists[l], shifted_constr_Element, result);
            SpaceGroup.SG_Multiply(result, k_inverse_Twists[l], rotated_constr_Element);

            // shift using non-standard lattice vectors
            for (n = 0; n < s3; ++n)
            {
              SpaceGroup.SG_Multiply(SG_AllNonStandardShifts[n], rotated_constr_Element, result);
              SpaceGroup.SG_Multiply(result, SG_All_inverse_NonStandardShifts[n], result2);

              if (ConstructingElementHasFP && SpaceGroup.SG_FixedPoint(result2, TransformedFixedPoint))
              {
                AddFixedPoint = true;
                
                s4 = FixedPoints.size();
                for (p = 0; AddFixedPoint && (p < s4); ++p)
                {
                  if (SpaceGroup.SG_DifferenceInLattice(FixedPoints[p], TransformedFixedPoint))
                    AddFixedPoint = false;
                }
                if (AddFixedPoint)
                {
                  FixedPoints.push_back(TransformedFixedPoint);
                  tmp1_CCs_of_constructing_Elements.push_back(result2);
                }
              }
              else
                tmp1_CCs_of_constructing_Elements.push_back(result2);
            }
          }
        }
      }
      if (tmp1_CCs_of_constructing_Elements.size() == 0)
        tmp1_CCs_of_constructing_Elements.push_back(ConstructingElement);

      tmp2_CCs_of_constructing_Elements.push_back(tmp1_CCs_of_constructing_Elements);
    }
    this->YC_CCs_of_constructing_Elements.push_back(tmp2_CCs_of_constructing_Elements);
  }
  // end: create the conjugacy classes for the space group selection rule
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  this->GI.LoadTensorProducts(Vacuum.SymmetryGroup.GaugeGroup);
    
  return true;
}



/* ########################################################################################
######   void CYukawaCouplings::Save(...) const                                      ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
void CYukawaCouplings::Save(const COrbifold &Orbifold, const SConfig &Vacuum, std::ostream &out, unsigned Order) const
{
  const COrbifoldGroup &OrbifoldGroup = Orbifold.OrbifoldGroup;
  CPrint Print(Tstandard, &out);

  // E8 x E8 or Spin32/Z2
  const SelfDualLattice Lattice = OrbifoldGroup.GetShift(0).Lattice;
  if (Lattice == E8xE8)
    out << "E8xE8\n";
  else
  if (Lattice == Spin32)
    out << "Spin32\n";

  Print.PrintRational(OrbifoldGroup.GetShift(0), Lattice);
  out << "\n";
  Print.PrintRational(OrbifoldGroup.GetShift(1), Lattice);
  out << "\n";
  Print.PrintWilsonLines(OrbifoldGroup.GetWilsonLines(), false);
  out << flush;

  unsigned i = 0;
  unsigned j = 0;
  
  const size_t s1 = Vacuum.FieldCouplings.size();
  size_t s2 = 0;

  if (Order == 0)
  {
    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &FieldCoupling = Vacuum.FieldCouplings[i].FieldIndices;

      s2 = FieldCoupling.size();

      for (j = 0; j < s2; ++j)
        out << FieldCoupling[j] << " ";
      out << "\n";
    }
  }
  else
  {
    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &FieldCoupling = Vacuum.FieldCouplings[i].FieldIndices;

      if (FieldCoupling.size() == Order)
      {
        for (j = 0; j < Order; ++j)
          out << FieldCoupling[j] << " ";
        out << "\n";
      }
    }
  }
  out << flush;
}



/* ########################################################################################
######   unsigned CYukawaCouplings::Load(...)                                        ######
######                                                                               ######
######   Version: 24.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
unsigned CYukawaCouplings::Load(const COrbifold &Orbifold, SConfig &Vacuum, std::ifstream &in)
{
  vector<string> linesOfData = return_next_lines(in,9);

  if (linesOfData.size() != 9)
  {
    cout << "\n  Warning in unsigned CYukawaCouplings::Load(...): File is empty. Return 0." << endl;
    return 0;
  }

  const vector<CField>  &Fields         = Vacuum.Fields;
  const COrbifoldGroup  &OrbifoldGroup  = Orbifold.OrbifoldGroup;
  const CSpaceGroup     &SpaceGroup     = OrbifoldGroup.GetSpaceGroup();
  const vector<CSector> &Spectrum       = Orbifold.GetSectors();
  
  const vector<SModularSymmetry> &ModularSymmetries = SpaceGroup.ModularSymmetries;
  const size_t number_of_modularsymmetries = ModularSymmetries.size();
  
  rational<int> exponent = -1;
  
  // E8 x E8 or Spin32/Z2
  SelfDualLattice Lattice = E8xE8;
  if (linesOfData[0] == "Spin32")
    Lattice = Spin32;
  else
  if (linesOfData[0] == "UNSPECIFIED_LATTICE")
  {
    cout << "\n  Warning in unsigned CYukawaCouplings::Load(...): Define the even and self-dual lattice. Return 0." << endl;
    return 0;
  }

  if (Lattice != OrbifoldGroup.GetShift(0).Lattice)
  {
    cout << "\n  Warning in unsigned CYukawaCouplings::Load(...): T^16 lattice of the file differs from the original definition. Return 0." << endl;
    return 0;
  }

  size_t t1 = 0;
  vector<rational<int> > rat_vector;
  CVector                com_vector(16);
  string::size_type loc1 = 0;
  string string_line = "";
  unsigned i = 0;
  unsigned j = 0;
  for (i = 1; i < 9; ++i)
  {
    string_line = linesOfData[i];
    
    // begin: remove brackets
    loc1 = 0;
    while (loc1 != string::npos)
    {
      loc1 = string_line.find("(", 0);
      if (loc1 != string::npos)
        string_line.erase(loc1,1);

      loc1 = string_line.find(")", 0);
      if (loc1 != string::npos)
        string_line.erase(loc1,1);
    }
    // end: remove brackets
    convert_string_to_vector_of_rational(string_line, rat_vector);
          
    if (rat_vector.size() != 16)
    {
      cout << "\n  Warning in unsigned CYukawaCouplings::Load(...): Vector does not have length 16. Return 0." << endl;
      return 0;
    }
    for (j = 0; j < 16; ++j)
      com_vector[j] = (double)rat_vector[j].numerator()/(double)rat_vector[j].denominator();

    bool error = false;
    if (i == 1)
    {
      if (com_vector != OrbifoldGroup.GetShift(0))
        error = true;
    }
    else
    if (i == 2)
    {
      if (com_vector != OrbifoldGroup.GetShift(1))
        error = true;
    }
    else
    {
      if (com_vector != OrbifoldGroup.GetWilsonLines().GetWilsonLine(i-3))
        error = true;
    }
    if (error)
    {
      cout << "\n  Warning in unsigned CYukawaCouplings::Load(...): Data of the model changed. Return 0." << endl;
      return 0;
    }
  }

  unsigned counter = 0;
  unsigned temp = 0;
  size_t Order = 1;

  vector<unsigned> currentData;

  YukawaCoupling tmp_YukawaCoupling;

  srand ( time(NULL) );

  if (Vacuum.FieldCouplings.size() != 0)
  {
    bool Coupling_known = false;
    t1 = 0;

    while (getline(in, string_line) )
    {
      currentData.clear();
      std::istringstream line(string_line);

      while (line >> temp)
        currentData.push_back(temp);

      Coupling_known = false;
      t1 = Vacuum.FieldCouplings.size();
      for (i = 0; !Coupling_known && (i < t1); ++i)
      {
        if (Vacuum.FieldCouplings[i].FieldIndices == currentData)
          Coupling_known = true;
      }
      if (!Coupling_known)
      {
        ++counter;

        tmp_YukawaCoupling.LabelOfModuli.clear();
        tmp_YukawaCoupling.ExponentsOfEtas.clear();

        Order = currentData.size();
        for (i = 0; i < number_of_modularsymmetries; ++i)
        {
          exponent = 1;
          for (j = 0; j < Order; ++j)
            exponent += Fields[currentData[j]].GetModularWeight(ModularSymmetries[i], Spectrum);

          exponent *= -2;

          tmp_YukawaCoupling.LabelOfModuli.push_back(ModularSymmetries[i].Label);
          tmp_YukawaCoupling.ExponentsOfEtas.push_back(exponent.numerator());
        }

        tmp_YukawaCoupling.CouplingStrength = -1.276453 * ((double)(rand() + 1))/((double)(RAND_MAX + 1));
        tmp_YukawaCoupling.FieldIndices     = currentData;

        Vacuum.FieldCouplings.push_back(tmp_YukawaCoupling);
      }
    }
  }
  else
  {
    while (getline(in, string_line) )
    {
      currentData.clear();
      std::istringstream line(string_line);

      while (line >> temp)
        currentData.push_back(temp);

      ++counter;

      tmp_YukawaCoupling.LabelOfModuli.clear();
      tmp_YukawaCoupling.ExponentsOfEtas.clear();

      Order = currentData.size();
      for (i = 0; i < number_of_modularsymmetries; ++i)
      {
        exponent = 1;
        for (j = 0; j < Order; ++j)
          exponent += Fields[currentData[j]].GetModularWeight(ModularSymmetries[i], Spectrum);

        exponent *= -2;

        tmp_YukawaCoupling.LabelOfModuli.push_back(ModularSymmetries[i].Label);
        tmp_YukawaCoupling.ExponentsOfEtas.push_back(exponent.numerator());
      }

      tmp_YukawaCoupling.CouplingStrength = -1.276453 * ((double)(rand() + 1))/((double)(RAND_MAX + 1));
      tmp_YukawaCoupling.FieldIndices     = currentData;

      Vacuum.FieldCouplings.push_back(tmp_YukawaCoupling);
    }
  }
  return counter;
}
