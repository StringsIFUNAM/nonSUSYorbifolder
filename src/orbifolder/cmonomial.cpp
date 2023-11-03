#include "cmonomial.h"
#include "cgaugeinvariance.h"
#include "cprompt.h"

using std::cout;
using std::endl;
using std::exit;



/* ########################################################################################
######   CMonomial(unsigned NumberOfU1s)                                             ######
######                                                                               ######
######   Version: 03.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) NumberOfU1s : number of U(1) factors of the gauge group                  ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CMonomial object.                                          ######
######################################################################################## */
CMonomial::CMonomial(unsigned NumberOfU1s)
  : U1Charges(NumberOfU1s)
{
  this->VEV_TurnedOn = false;
}



/* ########################################################################################
######   ~CMonomial()                                                                ######
######                                                                               ######
######   Version: 03.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CMonomial object.                                  ######
######################################################################################## */
CMonomial::~CMonomial()
{
}



/* ########################################################################################
######   operator*=(const CMonomial &Monomial2)                                      ######
######                                                                               ######
######   Version: 21.09.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Monomial2 : a CMonomial object to be multiplied                          ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the multiplication of two monomials, e,g,                           ######
######     M_1 = \phi_1^n_1 ... \phi_N^n_N and M_2 = \phi_1^m_1 ... \phi_N^m_N       ######
######     M_1 * M_2 = \phi_1^(n_1+m_1) ... \phi_N^(n_N+m_N)                         ######
######   and stpres the result in this monomial (e.g. in M_1).                       ######
######################################################################################## */
void CMonomial::operator*=(const CMonomial &Monomial2)
{
  if (this->VEV_TurnedOn || Monomial2.VEV_TurnedOn)
    cout << "Warning in void CMonomial::operator*=(const CMonomial &Monomial2): Multiplying monomials whose VEVs are turned on." << endl;

  unsigned i = 0;
  const size_t s1 = Monomial2.GaugeEquivalentFields.size();
  for (i = 0; i < s1; ++i)
  {
    const vector<unsigned> &GaugeEquivalentFields2_i = Monomial2.GaugeEquivalentFields[i];

    vector<vector<unsigned> >::iterator pos = find(this->GaugeEquivalentFields.begin(), this->GaugeEquivalentFields.end(), GaugeEquivalentFields2_i);
    if (pos == this->GaugeEquivalentFields.end())
    {
      this->GaugeEquivalentFields.push_back(GaugeEquivalentFields2_i);
      this->Exponents.push_back(Monomial2.Exponents[i]);
    }
    else
      this->Exponents[distance(this->GaugeEquivalentFields.begin(), pos)] += Monomial2.Exponents[i];
  }

  this->U1Charges += Monomial2.U1Charges;
}



/* ########################################################################################
######   CheckGaugeInvariance(const COrbifold &Orbifold,  ...) const                 ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold                     : the corresponding COrbifold object        ######
######   2) Vacuum                       : the SConfig object "Vacuum" contains the  ######
######                                     gauge group                               ######
######   3) demand_D0                    : check U(1)_anom gauge invariance?         ######
######   4) check_non_Abelian_invariance : check non-Abelien gauge invariance?       ######
######   output:                                                                     ######
######   return value                    : is this monomial gauge invariant (maybe   ######
######                                     up to the anomalous U(1))?                ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this monomial is gauge invariant using the gauge group       ######
######   stored in "Vacuum".                                                         ######
######################################################################################## */
bool CMonomial::CheckGaugeInvariance(const COrbifold &Orbifold, const SConfig &Vacuum, bool demand_D0, bool check_non_Abelian_invariance) const
{
  const double prec = 0.0001;
  unsigned i = 0;

  const double D0_FI_term = Vacuum.SymmetryGroup.D0_FI_term;

  size_t s1 = this->U1Charges.GetSize();
  // run through the U(1) charges of the monomial
  for (i = 0; i < s1; ++i)
  {
    // if the i-th U(1) charge is non-zero
    if (fabs(this->U1Charges[i]) > prec)
    {
      // only the first U(1) charge might be non-zero
      if ((i != 0) || (D0_FI_term == 0) || (demand_D0 && ((this->U1Charges[0] * D0_FI_term) > prec)))
        return false;
    }
  }
  if (!check_non_Abelian_invariance)
    return true;

  vector<unsigned> FieldCoupling;
  
  s1 = this->GaugeEquivalentFields.size();
  for (i = 0; i < s1; ++i)
    FieldCoupling.insert(FieldCoupling.end(), this->Exponents[i], this->GaugeEquivalentFields[i][0]);
  
  // N O N - A B E L I A N  G A U G E Invariance selection rule
  CGaugeInvariance GI;
  GI.LoadTensorProducts(Vacuum.SymmetryGroup.GaugeGroup);

  if (!GI.CheckInvariance(Orbifold, Vacuum, FieldCoupling))
    return false;

  return true;
}



/* ########################################################################################
######   ContainsField(unsigned FieldIndex) const                                    ######
######                                                                               ######
######   Version: 20.09.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldIndex : a field-index giving the position in vector<CField> &Fields ######
######   output:                                                                     ######
######   return value  : is field with index "FieldIndex" contained in the current   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this monimial contains the field specified by the field      ######
######   index "FieldIndex".                                                         ######
######################################################################################## */
bool CMonomial::ContainsField(unsigned FieldIndex) const
{
  const size_t s1 = GaugeEquivalentFields.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    const vector<unsigned> &GaugeEquivalentFields_i = GaugeEquivalentFields[i];

    if (find(GaugeEquivalentFields_i.begin(), GaugeEquivalentFields_i.end(), FieldIndex) != GaugeEquivalentFields_i.end())
      return true;
  }
  return false;
}



/* ########################################################################################
######   ReadMonomial(string input, const vector<CField> &Fields,unsigned use_Labels)######
######                                                                               ######
######   Version: 29.09.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input        : input string, e.g. "D_1 bD_4 S_1 S_2^4" or                ######
######                     with brackets "D_1 bD_4 S_1 (S_2)^4                       ######
######   2) Fields       : vector of all (massless) fields of the current vacuum     ######
######   3) use_Labels   : specifies which label should be used                      ######
######   output:                                                                     ######
######   return bool     : does the input string contain a valid monomial?           ######
###########################################################################################
######   description:                                                                ######
######   Read the monomial from the string "input".                                  ######
######################################################################################## */
bool CMonomial::ReadMonomial(string input, const vector<CField> &Fields, unsigned use_Labels)
{
  // begin: erase the brackets the from string "input"
  string::size_type loc1 = input.find_first_of("()");
  while (loc1 != string::npos)
  {
    input = input.substr(0, loc1) + input.substr(loc1+1, string::npos);
    loc1 = input.find_first_of("()");
  }
  // end: erase the brackets from the string "input"

  this->U1Charges.assign(this->U1Charges.size(), 0);

  this->Exponents.clear();
  this->GaugeEquivalentFields.clear();
  this->VEV_TurnedOn = false;

  string tmp_string = "";
  vector<string> FieldLabels;
  unsigned Exp = 0;

  loc1 = 0;
  string::size_type loc2 = 0;
  string::size_type loc3 = 0;
  string::size_type loc4 = 0;
  vector<unsigned> tmp_Indices;

  while (loc1 != string::npos)
  {
    // find the brackets { and } that create a new set of gauge-equivalent fields
    loc1 = input.find(" ", loc2);
    loc3 = input.find("{", loc2);

    // if the next bracket { is before the next space, a set of gauge-equivalent fields follows
    if (loc3 < loc1)
    {
      loc4 = input.find("}", loc3);
      if (loc4 == string::npos)
        return false;

      loc1 = input.find(" ", loc4);
      tmp_string = input.substr(loc3, loc1 - loc3);

      // begin: erase the brackets from the string "tmp_string1"
      loc4 = tmp_string.find_first_of("{}");
      while (loc4 != string::npos)
      {
        tmp_string = tmp_string.substr(0, loc4) + tmp_string.substr(loc4+1, string::npos);
        loc4 = tmp_string.find_first_of("{}");
      }
      // end: erase the brackets from the string "tmp_string1"
    }
    else
      tmp_string = input.substr(loc2, loc1 - loc2);

    // begin: find the exponent
    loc3 = tmp_string.find("^", 0);
    if (loc3 != string::npos)
    {
      string tmp = tmp_string.substr(loc3+1, string::npos);
      if (tmp.find_first_not_of("0123456789") != string::npos)
        return false;

      Exp = (unsigned)atoi(tmp.c_str());
      if (Exp == 0)
        return false;

      tmp_string = tmp_string.substr(0, loc3);
    }
    else
      Exp = 1;
    // end: find the exponent

    FieldLabels.clear();
    global_DecomposeString(tmp_string, FieldLabels, " ");
    
    // find the field indices
    tmp_Indices = global_GetIndicesOnlyFieldWithNumber(FieldLabels, use_Labels, Fields);
    if (tmp_Indices.size() == 0)
      return false;

    this->GaugeEquivalentFields.push_back(tmp_Indices);
    this->Exponents.push_back(Exp);
    this->U1Charges += (Fields[this->GaugeEquivalentFields.rbegin()->at(0)].U1Charges * Exp);

    loc2 = loc1 + 1;
  }
  return true;
}



/* ########################################################################################
######   SetVEVs(double VEV_Scale, vector<CField> &Fields)                           ######
######                                                                               ######
######   Version: 03.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEV_Scale : e.g. the VEVs for "D_1 bD_4 S_1 S_2^4" will be               ######
######                  <D_1> = <bD_4> = <S_1> = VEV_Scale and                       ######
######                  <S_2> = 2 VEV_Scale                                          ######
######   2) Fields    : list of CField objects; the VEVs of the fields contained     ######
######                  in this monomial will be set (using "VEV_Scale")             ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Turn on the VEVs for all fields contained in this monomial.The VEV depends  ######
######   on the overall scale "VEV_Scale" and the fields' exponent in the monomial.  ######
######################################################################################## */
void CMonomial::SetVEVs(double VEV_Scale, vector<CField> &Fields)
{
  this->VEV_TurnedOn = (fabs(VEV_Scale) > 0.000001);

  const size_t s1 = this->GaugeEquivalentFields.size();
  size_t s2 = 0;

  unsigned i = 0;
  unsigned j = 0;

  unsigned tmp = 0;
  for (i = 0; i < s1; ++i)
  {
    const vector<unsigned> &GaugeEquivalentFields_i = GaugeEquivalentFields[i];
    tmp = this->Exponents[i];

    s2 = GaugeEquivalentFields_i.size();
    for (j = 0; j < s2; ++j)
    {
      CField &tmp_Field = Fields[GaugeEquivalentFields_i[j]];

      if (this->VEV_TurnedOn)
        tmp_Field.VEVs[0] = VEV_Scale * sqrt((double)tmp);
      else
        tmp_Field.VEVs.SetToZero();
    }
  }
}
