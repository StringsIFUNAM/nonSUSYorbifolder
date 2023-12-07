
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <vector>

#include "globalfunctions.h"
#include "corbifold.h"
#include "clinalg.h"
#include "clatticevector.h"








#define CHECKERROR true

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using boost::rational;

const string AllCharacters  = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()-_=+[{]};:',<.>/?|~ ";



/* ########################################################################################
######   global_ReplaceString(string &input, const string &replace_from, ...)        ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input                 : look for the substring "replace_from" in "input" ######
######                              and replace it by "replace_to"                   ######
######   2) replace_from          : string to be replaced                            ######
######   3) replace_to            : string to be inserted instead of "replace_from"  ######
######                                                                               ######
######   output:                                                                     ######
######   -                                                                           ######
######   return:                                                                     ######
######   bool                     : was something replaced?                          ######
######################################################################################## */
bool global_ReplaceString(string &input, const string &replace_from, const string &replace_to)
{
  unsigned counter = 0;

  string::size_type loc1 = input.find(replace_from, 0);
  while (loc1 != string::npos)
  {
    input.erase(loc1, replace_from.size());
    input.insert(loc1, replace_to);

    ++counter;
    if (counter > 500)
    {
      cout << "\n  Warning in bool global_ReplaceString(...) : \"replace_from\" found too often. Return true." << endl;
      return true;
    }

    loc1 = input.find(replace_from, 0);
  }
  return (counter != 0);
}



/* ########################################################################################
######   global_DecomposeString(const string &input, vector<string> &output, ...)    ######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input                 : split the string "input" into its substrings and ######
######                              store them in "output"                           ######
######   3) separator             : a single character that seperates two substrings ######
######                              e.g. " " or ","                                  ######
######   output:                                                                     ######
######   2) output                : vector of strings containing the result          ######
######   return:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void global_DecomposeString(const string &input, vector<string> &output, const string separator)
{
  #ifdef CHECKERROR
  if (separator.size() != 1)
  {
    cout << "\n  Warning in void global_DecomposeString(...) : the \"separator\" is ill-defined. Return." << endl;
    return;
  }
  #endif

  string::size_type loc1 = 0;
  string::size_type loc2 = 0;
  while (loc1 != string::npos)
  {
    loc1 = input.find(separator, loc2);
    output.push_back(input.substr(loc2, loc1 - loc2));
    loc2 = loc1 + 1;
  }
}



/* ########################################################################################
######   global_ConvertVectorOfString2VectorOfUnsigned(...)                          ######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input                 : vector of strings, each string should be an      ######
######                              unsigned integer                                 ######
######   output:                                                                     ######
######   2) output                : vector of unsigned integers containing the       ######
######                              result                                           ######
######   return:                                                                     ######
######   bool                     : function terminated correctly?                   ######
######################################################################################## */
bool global_ConvertVectorOfString2VectorOfUnsigned(const vector<string> &input, vector<unsigned> &output)
{
  string tmp = "";
  for (vector<string>::const_iterator it = input.begin(); it < input.end(); ++it)
  {
    tmp = *it;
    if (tmp.find_first_not_of("0123456789") == string::npos)
      output.push_back((unsigned)atoi(tmp.c_str()));
    #ifdef CHECKERROR
    else
    {
      cout << "\n  Warning in void global_ConvertVectorOfString2VectorOfUnsigned(...) : string contains non-unsigned-integer. Return false." << endl;
      return false;
    }
    #endif
  }
  return true;
}



/* ########################################################################################
######   global_GetIndices(const vector<string> &FieldLabels, unsigned use_Labels...)######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldLabels           : vector of strings, each string is a field label, ######
######                              e.g. either "F_1" or just "F"                    ######
######   2) use_Labels            : which labels are used                            ######
######   3) Fields                : vector of CField objects, where to look for      ######
######                              "FieldLabels"                                    ######
######   output:                                                                     ######
######   -                                                                           ######
######   return:                                                                     ######
######   vector<unsigned>         : the indices of the fields with labels            ######
######                              "FieldLabels"                                    ######
######################################################################################## */
vector<unsigned> global_GetIndices(const vector<string> &FieldLabels, unsigned use_Labels, const vector<CField> &Fields)
{
  vector<unsigned> result;

  string label1 = "";
  string label2 = "";

  const size_t f1 = Fields.size();
  for (unsigned i = 0; i < f1; ++i)
  {
    const CField &Field = Fields[i];

    label1 = Field.Labels[use_Labels];

    label2 = label1;
    label2 += "_";
    std::ostringstream os;
    os << Field.Numbers[use_Labels];
    label2 += os.str();

    if ((find(FieldLabels.begin(), FieldLabels.end(), label1) != FieldLabels.end()) || (find(FieldLabels.begin(), FieldLabels.end(), label2) != FieldLabels.end()))
      result.push_back(i);
  }
  return result;
}



/* ########################################################################################
######   global_GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels, ...)######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldLabels           : vector of strings, each string is a field label, ######
######                              here only (e.g.) "F_1" is possible               ######
######   2) use_Labels            : which labels are used                            ######
######   3) Fields                : vector of CField objects, where to look for      ######
######                              "FieldLabels"                                    ######
######   output:                                                                     ######
######   -                                                                           ######
######   return:                                                                     ######
######   vector<unsigned>         : the indices of the fields with labels            ######
######                              "FieldLabels"                                    ######
######################################################################################## */
vector<unsigned> global_GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels, unsigned use_Labels, const vector<CField> &Fields)
{
  vector<unsigned> result;

  string label = "";

  const size_t f1 = Fields.size();
  for (unsigned i = 0; i < f1; ++i)
  {
    const CField &Field = Fields[i];

    label = Field.Labels[use_Labels];
    label += "_";
    std::ostringstream os;
    os << Field.Numbers[use_Labels];
    label += os.str();

    if (find(FieldLabels.begin(), FieldLabels.end(), label) != FieldLabels.end())
      result.push_back(i);
  }
  return result;
}



/* ########################################################################################
######   global_LoadSGElement(string &input, CSpaceGroupElement &SGElement)          ######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input                 : a string that contains 3 integers and            ######
######                              6 rational numbers                               ######
######   output:                                                                     ######
######   2) SGElement             : a CSpaceGroupElement object to store the result  ######
######   return:                                                                     ######
######   bool                     : function terminated correctly?                   ######
######################################################################################## */
bool global_LoadSGElement(const string &input, CSpaceGroupElement &SGElement)
{
  bool SGElement_OK = true;
  rationalVector RationalVector;
  convert_string_to_vector_of_rational(input, RationalVector);

  if (RationalVector.size() != 9) {						//hacking here!!!
    SGElement_OK = false;
  }
  else
  {
    if (RationalVector[0].denominator() != 1)
      SGElement_OK = false;
    else
    {
      if (RationalVector[1].denominator() != 1)
        SGElement_OK = false;
    }
  }

  if (!SGElement_OK)
  {
    cout << "\n  Warning in bool global_LoadSGElement(...) : space group element ill specified. Return false." << endl;
    return false;
  }

  SGElement.Set_m(RationalVector[0].numerator());		//hacking here!!!
  SGElement.Set_n(RationalVector[1].numerator());
  SGElement.Set_k(RationalVector[2].numerator());

  for (unsigned k = 0; k < 6; ++k)
    SGElement.Set_n_alpha(k, RationalVector[k+3]);		//hacking here!!!

  return true;
}



/* ##########################################################################
######   rational<int> D2Rat(const double x)                           ######
######   Ist die routine falsch? Saul war's!                           ######
######   Version: 18.7.2006                                            ######
########################################################################## */
rational<int> D2Rat(const double x)
{
  double tmp=0.0;

  for(int i=1; i<10000; ++i)
  {
    tmp = x * i;
    if(is_integer(tmp))
      return rational<int>((int)roundf(tmp),i);
  }

  return rational<int>(((int)roundf(10000*x)),10000);
}


/* ##########################################################################
######   NewPrompt                                                     ######
######   Added the ability to navigate in the prompt                   ######
######   Version: 06.11.2023                                           ######
########################################################################## */

//








istream& GetSaveLine(istream &in, string &result)
{
  istream &tmp = getline(in, result);

  // begin: remove " " at the beginning
  string::size_type loc1 = result.find_first_not_of(" ");
  if (loc1 == string::npos)
  {
    result = "";
    return tmp;
  }
  result = result.substr(loc1, string::npos);
  // end: remove " " at the beginning

  // begin: remove unwanted characters
  while (true)
  {
    string::size_type loc1 = result.find_first_not_of(AllCharacters);
    if (loc1 == string::npos)
      break;
    result.erase(loc1,1);
  }
  // end: remove unwanted characters

  // begin: remove " " at the end
  loc1 = result.find_last_not_of(" ");
  if (loc1 != string::npos)
    result.resize(loc1+1);
  // end: remove " " at the end

  return tmp;
}


bool IsSinglet(const RepVector &Dimensions)
{
  const size_t s1 = Dimensions.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if (Dimensions[i].Dimension != 1)
      return false;
  }
  return true;
}



bool CC(const RepVector &Dimensions, RepVector &DimensionsCC)
{
  const size_t s1 = Dimensions.size();
  if (s1 == 0)
    return false;
  
  for (unsigned i = 0; i < s1; ++i)
  {
    if ((Dimensions[i].Dimension != 1) && (Dimensions[i].Dimension != 2))
    {
      DimensionsCC[i].Dimension = -Dimensions[i].Dimension;
      DimensionsCC[i].AdditionalLabel = Dimensions[i].AdditionalLabel;
    }
    else
    {
      DimensionsCC[i].Dimension = Dimensions[i].Dimension;
      DimensionsCC[i].AdditionalLabel = Dimensions[i].AdditionalLabel;
    }
  }
  return true;
}


bool AreHighestWeights_CC(const vector<int> &HighestWeight_DL_1, const vector<int> &HighestWeight_DL_2, const gaugeGroupFactor<double> &GGFactor)
{
  if ((GGFactor.rank != HighestWeight_DL_1.size()) || (GGFactor.rank != HighestWeight_DL_2.size()))
  {
    cout << "\n  Warning in bool AreHighestWeights_CC(...) : the rank of the gauge group and the size of the weights does not match. Return false." << endl;
    return false;
  }
  if (GGFactor.algebra[0] == 'A')
  {
    vector<int> CC_HighestWeight_DL_1 = HighestWeight_DL_1;

    if (GGFactor.rank != 1)
      reverse(CC_HighestWeight_DL_1.begin(), CC_HighestWeight_DL_1.end());

    return (CC_HighestWeight_DL_1 == HighestWeight_DL_2);
  }
  else
  if (GGFactor.algebra[0] == 'D')
  {
    unsigned rank = GGFactor.rank;
    vector<int> tmp_HighestWeight_DL(rank, 0);

    // singlet rep.
    if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
      return true;

    // vector rep.
    tmp_HighestWeight_DL.assign(rank, 0);
    tmp_HighestWeight_DL[0] = 1;

    if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
      return true;

    // spinor rep.
    if ((rank == 4) || (rank == 6) || (rank == 8))
    {
      tmp_HighestWeight_DL.assign(rank, 0);
      tmp_HighestWeight_DL[rank-1] = 1;
      
      if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
        return true;
      
      tmp_HighestWeight_DL.assign(rank, 0);
      tmp_HighestWeight_DL[rank-2] = 1;
      
      if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
        return true;
    }
    else
    {
      tmp_HighestWeight_DL.assign(rank, 0);
      tmp_HighestWeight_DL[rank-1] = HighestWeight_DL_1[rank-2];
      tmp_HighestWeight_DL[rank-2] = HighestWeight_DL_1[rank-1];
      if (tmp_HighestWeight_DL == HighestWeight_DL_2)
        return true;
    }

    // adjoint rep.
    tmp_HighestWeight_DL.assign(rank, 0);
    tmp_HighestWeight_DL[1] = 1;

    if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
      return true;

    return false;
  }
  else
  if (GGFactor.algebra[0] == 'E')
  {
    vector<int> tmp_HighestWeight_DL(GGFactor.rank, 0);

    // singlet rep.
    if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
      return true;

    switch (GGFactor.rank)
    {
      case 6:
      {
        // 27-plet (complex)
        tmp_HighestWeight_DL.assign(6, 0);
        tmp_HighestWeight_DL[0] = HighestWeight_DL_1[4];
        tmp_HighestWeight_DL[4] = HighestWeight_DL_1[0];
        if (tmp_HighestWeight_DL == HighestWeight_DL_2)
          return true;

        // 78 plet (real)
        tmp_HighestWeight_DL.assign(6, 0);
        tmp_HighestWeight_DL[5] = 1;
        if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
          return true;

        return false;
      }
      case 7:
      {
        // 56-plet (real)
        tmp_HighestWeight_DL.assign(7, 0);
        tmp_HighestWeight_DL[5] = 1;
        if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
          return true;

        // 133-plet (real)
        tmp_HighestWeight_DL.assign(7, 0);
        tmp_HighestWeight_DL[0] = 1;
        if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
          return true;

        return false;
      }
      case 8:
      {
        // 248-plet (real)
        tmp_HighestWeight_DL.assign(8, 0);
        tmp_HighestWeight_DL[6] = 1;
        if ((tmp_HighestWeight_DL == HighestWeight_DL_1) && (tmp_HighestWeight_DL == HighestWeight_DL_2))
          return true;

        return false;
      }
      default:
      {
        cout << "\n  Warning in bool AreHighestWeights_CC(...) : case \"" << GGFactor.algebra << "\" not defined. Return false." << endl;
        return false;
      }
    }
  }
  cout << "\n  Warning in bool AreHighestWeights_CC(...) : case not defined. Return false." << endl;
  return false;
}


bool operator==(const RepVector &Dimensions1, const RepVector &Dimensions2)
{
  return AreRepVectorsEqual(Dimensions1, Dimensions2);
}

bool AreRepVectorsEqual(const RepVector &Dimensions1, const RepVector &Dimensions2)
{
  const size_t s1 = Dimensions1.size();
  if (s1 != Dimensions2.size())
  {
    cout << "\n  Warning in bool AreRepVectorsEqual(...) : RepVectors have different length. Return false." << endl;
    return false;
  }
  
  for (unsigned i = 0; i < s1; ++i)
  {
    if ((Dimensions1[i].Dimension != Dimensions2[i].Dimension) || (Dimensions1[i].AdditionalLabel != Dimensions2[i].AdditionalLabel))
      return false;
  }
  return true;
}


bool AreRepVectorsEqual(const SSymmetryGroup &SymmetryGroup, const RepVector &Dimensions1, const RepVector &Dimensions2)
{
  if (Dimensions1.size() != Dimensions2.size())
  {
    cout << "\n  Warning in bool AreRepVectorsEqual(...) : RepVectors have different length. Return false." << endl;
    return false;
  }

  unsigned index = 0;
  const size_t s1 = SymmetryGroup.observable_sector_GGs.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    index = SymmetryGroup.observable_sector_GGs[i];
    if ((Dimensions1[index].Dimension != Dimensions2[index].Dimension) || (Dimensions1[index].AdditionalLabel != Dimensions2[index].AdditionalLabel))
      return false;
  }
  return true;
}



unsigned Multiplicity(const SSymmetryGroup &SymmetryGroup, const RepVector &Dimensions)
{
  unsigned mult = 1;
  const size_t s1 = Dimensions.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if (find(SymmetryGroup.observable_sector_GGs.begin(), SymmetryGroup.observable_sector_GGs.end(), i) == SymmetryGroup.observable_sector_GGs.end())
      mult *= abs(Dimensions[i].Dimension);
  }
  return mult;
}



bool AreU1ChargesEqual(const SSymmetryGroup &SymmetryGroup, const CVector &U1Charges1, const CVector &U1Charges2)
{
  if (U1Charges1.size() != U1Charges2.size())
  {
    cout << "\n  Warning in bool AreRepVectorsEqual(...) : U(1) charge vectors have different length. Return false." << endl;
    return false;
  }
  const double prec(0.0001);

  unsigned index = 0;
  const size_t s1 = SymmetryGroup.observable_sector_U1s.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    index = SymmetryGroup.observable_sector_U1s[i];
    if (fabs(U1Charges1[index] - U1Charges2[index]) > prec)
      return false;
  }
  return true;
}



unsigned ggT(unsigned m, unsigned n)
{
  unsigned z;
  while (n != 0)
  {
    z=n;
    n=m%n;
    m=z;
  }
  return m;
}

unsigned kgV(unsigned m, unsigned n)
{
  return (m*n/ggT(m,n));
}


/* ########################################################################################
######   void RecursiveCounting(...)                                                 ######
######                                                                               ######
######   Version: 14.08.2008                                                         ######
######################################################################################## */
void RecursiveCounting(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<vector<unsigned> > &AllNumbers)
{
  const size_t number_of_Digits = currentNumber.size();

  // runs through the positions of the current number
  for (unsigned i = 0; i < currentMaxDigit; ++i)
  {
    currentNumber[currentPosition] = i;

    unsigned nextPosition = currentPosition + 1;

    if (nextPosition < number_of_Digits)
      RecursiveCounting(currentNumber, nextPosition, MaxDigits[nextPosition], MaxDigits, AllNumbers);
    else
      AllNumbers.push_back(currentNumber);
  }
}




/* ########################################################################################
######   void RecursiveCounting(...)                                                 ######
######                                                                               ######
######   Version: 26.09.2006                                                         ######
######################################################################################## */
void RecursiveCounting(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, const vector<vector<unsigned> > &vec_adjoining_Positions, vector<vector<unsigned> > &AllNumbers)
{
  const size_t number_of_Digits = currentNumber.size();
  
  // runs through the positions of the current number
  for (unsigned i = 0; i < currentMaxDigit; ++i)
  {
    currentNumber[currentPosition] = i;
    
    unsigned nextPosition = currentPosition + 1;
    unsigned nextMaxDigit = 0;

    if (nextPosition < number_of_Digits)
    {
      const unsigned pod = MaxDigits[nextPosition];

      bool from_same_DigitPool = true;
      const size_t s1 = vec_adjoining_Positions.size();

      if (s1 == 0)
        from_same_DigitPool = false;

      for (unsigned j = 0; from_same_DigitPool && (j < s1); ++j)
      {
        const vector<unsigned> &adjoining_Positions = vec_adjoining_Positions[j];
        if (adjoining_Positions[nextPosition] != adjoining_Positions[currentPosition])
          from_same_DigitPool = false;
      }
        
      if (from_same_DigitPool)
      {
        if (i+1 < pod)
          nextMaxDigit = i+1;
        else
          nextMaxDigit = pod;
      }
      else
        nextMaxDigit = pod;

      RecursiveCounting(currentNumber, nextPosition, nextMaxDigit, MaxDigits, vec_adjoining_Positions, AllNumbers);
    }
    else
      AllNumbers.push_back(currentNumber);
  }
}


/* ########################################################################################
######   void RecursiveCounting_to_files(...)                                        ######
######                                                                               ######
######   Version: 02.10.2006                                                         ######
######################################################################################## */
void RecursiveCounting_to_files(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<std::ofstream*> &outs, unsigned long &counter)
{
  const size_t number_of_Digits = currentNumber.size();
  
  // runs through the positions of the current number
  for (unsigned i = 0; i < currentMaxDigit; ++i)
  {
    currentNumber[currentPosition] = i;

    unsigned nextPosition = currentPosition + 1;

    if (nextPosition < number_of_Digits)
      RecursiveCounting_to_files(currentNumber, nextPosition, MaxDigits[nextPosition], MaxDigits, outs, counter);
    else
    {
      ++counter;
      unsigned pos = counter/10000000;
      if (pos >= outs.size())
      {
        cout << "Error in void RecursiveCounting_to_files(...): Not enough files. Return." << endl;
        return;
      }
      std::ostream &out = *(outs[pos]);
      for (unsigned j = 0; j < number_of_Digits; ++j)
        out << currentNumber[j] << " ";
      out << "\n";
    }
  }
}



/* ########################################################################################
######   void NextNumber(...)                                                        ######
######                                                                               ######
######   Version: 13.11.2008                                                         ######
######################################################################################## */
bool NextNumber(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const size_t &number_of_Digits)
{
  for (int i = number_of_Digits-1; i >= 0; --i)
  {
    if (++Number[i] < MaxDigits[i])
      return true;
    else
      Number[i] = 0;
  }
  return false;
}



/* ########################################################################################
######   void NextNumber(...)                                                        ######
######                                                                               ######
######   Version: 13.11.2008                                                         ######
######                                                                               ######
######   starts to count at "0 ... 0 1" when "Number = 0 ... 0"                      ######
######################################################################################## */
bool NextNumber(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const size_t &number_of_Digits, const unsigned &OnlyChangeUpToDigit)
{
  for (unsigned i = number_of_Digits - 1; i > OnlyChangeUpToDigit; --i)
  {
    if (++Number[i] < MaxDigits[i])
      return true;
    else
      Number[i] = 0;
  }
  if (++Number[OnlyChangeUpToDigit] < MaxDigits[OnlyChangeUpToDigit])
    return true;
  else
    Number[OnlyChangeUpToDigit] = 0;
  return false;
}



/* ########################################################################################
######   void NextNumberNoTwins(...)                                                 ######
######                                                                               ######
######   Version: 24.10.2008                                                         ######
######################################################################################## */
bool NextNumberNoTwins(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const vector<vector<unsigned> > &PositionsOfUnequalDigits, const unsigned &OnlyChangeUpToDigit)
{
  unsigned i = 0;
  bool go_on = true;

  const size_t s1 = Number.size();
  const size_t s2 = PositionsOfUnequalDigits.size();

  while (NextNumber(Number, MaxDigits, s1, OnlyChangeUpToDigit))
  {
    go_on = true;
    for (i = 0; go_on && (i < s2); ++i)
    {
      if (Number[PositionsOfUnequalDigits[i][0]] == Number[PositionsOfUnequalDigits[i][1]])
        go_on = false;
    }
    if (go_on)
      return true;
  }
  return false;
}



void round_double_to_int(const double &d, int &i)
{
  const double tmp_floor = floor(d);

  if (fabs(tmp_floor - d) < 0.49999)
  {
    i = (int)tmp_floor;
    return;
  }

  i = (int)ceil(d);
  return;
}


int round_double_to_int(const double &d)
{
  const double tmp_floor = floor(d);

  if (fabs(tmp_floor - d) < 0.49999)
    return (int)tmp_floor;

  return (int)ceil(d);
}


bool is_integer(const double &d)
{
  return ((fabs(floor(d) - d) < 0.00001) || (fabs(ceil(d) - d) < 0.00001));
}


bool is_even(const double &d)
{
  const double d2 = d / 2.0;
  return ((fabs(floor(d2) - d2) < 0.00001) || (fabs(ceil(d2) - d2) < 0.00001));
}


// round double to double
void RoundDouble(double &tmp)
{
  const double prec = 0.000001;

  // if integer
  if ((fabs(floor(tmp) - tmp) < prec) || (fabs(ceil(tmp) - tmp) < prec))
  {
    const double tmp_floor = floor(tmp);

    if (fabs(tmp_floor - tmp) < 0.499999)
    {
      tmp = tmp_floor;
      return;
    }

    tmp = ceil(tmp);
    return;
  }

  double num = 2.0 * tmp;
  for (unsigned i = 2; i < 5000; ++i)
  {
    if(fabs(round(num) - num) < prec)
    {
      tmp = round(num)/(double)i;
      return;
    }
    num += tmp;
  }

  cout << "\n  Warning void RoundDouble(double &tmp) : cannot round " << tmp << "." << endl;
}



/* ########################################################################################
######   RoundCharge(long double &U1_Charge) const                                   ######
######                                                                               ######
######   Version: 21.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) U1_Charge : the U(1) charge as a long double that will be rounded        ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
######################################################################################## */
bool RoundCharge(long double &U1_Charge)
{
  const double prec = 0.0001;
  long double copy_U1_Charge  = 0.0;
  long double floor_U1_Charge = 0.0;
  long double diff            = 0.0;

  unsigned denominator = 2;
  unsigned numerator   = 1;

  // begin: round the charge
  copy_U1_Charge  = U1_Charge;
  floor_U1_Charge = floor(U1_Charge);
  diff            = U1_Charge - floor_U1_Charge;

  if ((diff > -0.0001) && (diff <  0.0001))
    U1_Charge = floor_U1_Charge;
  else
  if ((diff >  0.9999) && (diff <  1.0001))
    U1_Charge = floor_U1_Charge + 1.0;
  else
  if ((diff < -0.9999) && (diff > -1.0001))
    U1_Charge = floor_U1_Charge - 1.0;
  else
  {
    for (denominator = 2; denominator <= 800; ++denominator)
    {
      for (numerator = 1; numerator < denominator; ++numerator)
      {
        if (fabs(diff - ((double)numerator/(double)denominator)) < prec)
        {
          U1_Charge = floor_U1_Charge + (double)numerator/(double)denominator;

          if (fabs(copy_U1_Charge - U1_Charge) > prec)
          {
            cout << "Warning in bool RoundCharge(...) const: Check the rounded U(1) charge. Return false." << endl;
            return false;
          }
          return true;
        }
      }
    }

    cout << "Warning in bool RoundCharge(...) const: U1 charge not known: \n";
    cout << setw(15) << setprecision(12) << "U1_Charge        = " << U1_Charge << "\nfloor(U1_Charge) = " << floor_U1_Charge << "\ndiff             = " << diff << "\n";
    return false;
  }
  // end: round the charge
  if (fabs(copy_U1_Charge - U1_Charge) > prec)
  {
    cout << "Warning in bool RoundCharge(...) const: Check the rounded U(1) charge. Return false." << endl;
    return false;
  }
  return true;
}




/* ##########################################################################
######   rational<CHugeInt> D2HugeInt(const double x)                  ######
######                                                                 ######
######   Version: 30.3.2009                                            ######
########################################################################## */
rational<CHugeInt> D2HugeInt(const double x)
{
  const double prec = 0.000001;
  double tmp = 0.0;

  for(unsigned i = 1; i < 500; ++i)
  {
    tmp = x * i;
    if(fabs(roundf(tmp) - tmp) < prec)
      return rational<CHugeInt>((CHugeInt)((int)roundf(tmp)), CHugeInt(i)); //return rational<CHugeInt>((CHugeInt)((int)roundf(tmp)),i);
  }
  cout << "Error in globalfunctions:D2HugeInt(...): no conversion from double to CHugeInt found. Return 0." << endl;
  return rational<CHugeInt>(CHugeInt(0), CHugeInt(1)); //return rational<CHugeInt>(0,1);
}



bool clean_string(string &input)
{
  // begin: remove " " at the beginning
  string::size_type pos = input.find_first_not_of(" ");
  if (pos == string::npos)
  {
    input = "";
    return true;
  }
  input = input.substr(pos, string::npos);
  // end: remove " " at the beginning

  while (global_ReplaceString(input, ",", " "))
  {}

  while (global_ReplaceString(input, "  ", " "))
  {}

  return true;
}



// a clean string has no brackets and a single space as seperator
bool convert_clean_string_to_vector_of_int(string input, vector<int> &output)
{
  output.clear();

  string tmp_string1 = "";
  string tmp_string2 = "";
  string::size_type pos1 = 0;
  string::size_type pos2 = 0;
  string::size_type last_pos = 0;

  int number = 0;
  int exp    = 0;

  unsigned i = 0;

  while (pos1 != string::npos)
  {
    pos1 = input.find(" ", last_pos);

    tmp_string1 = input.substr(last_pos, pos1 - last_pos);
    pos2 = tmp_string1.find("^", 0);

    // begin: read the integer
    tmp_string2 = tmp_string1.substr(0, pos2);
    if ((tmp_string2 == "") || (tmp_string2.find_first_not_of("+-0123456789", 0) != string::npos))
      return false;

    number = atoi(tmp_string2.c_str());
    // end: read the integer

    if (pos2 == string::npos)
      output.push_back(number);
    else
    {
      // begin: read the exponent
      tmp_string2 = tmp_string1.substr(pos2+1, string::npos);
      if ((tmp_string2 == "") || (tmp_string2.find_first_not_of("0123456789", 0) != string::npos))
        return false;

      exp = (unsigned)atoi(tmp_string2.c_str());
      // end: read the exponent

      for (i = 0; i < exp; ++i)
        output.push_back(number);
    }
    last_pos = pos1 + 1;
  }
  return true;
}


bool convert_string_to_vector_of_int(string input, vector<int> &output)
{
  output.clear();

  const bool b1 = clean_string(input);

  if (input == "")
    return true;

  const bool b2 = convert_clean_string_to_vector_of_int(input, output);

  return (b1 && b2);
}


bool convert_clean_string_to_rational(string input, rational<int> &output)
{
  string::size_type pos1 = input.find("/", 0);
  if (pos1 == string::npos)
  {
    if ((input == "") || (input.find_first_not_of("+-0123456789", 0) != string::npos))
    {
      cout << "\n  Warning in bool convert_clean_string_to_rational(...) : input string ill-defined. Return false." << endl;
      return false;
    }
    output = rational<int>(atoi(input.c_str()),1);
    return true;
  }
  
  string tmp_string = input.substr(0, pos1);
  if ((tmp_string == "") || (tmp_string.find_first_not_of("+-0123456789", 0) != string::npos))
  {
    cout << "\n  Warning in bool convert_clean_string_to_rational(...) : input string ill-defined. Return false." << endl;
    return false;
  }
  int num = atoi(tmp_string.c_str());

  tmp_string = input.substr(pos1+1, string::npos);
  if ((tmp_string == "") || (tmp_string.find_first_not_of("0123456789", 0) != string::npos))
  {
    cout << "\n  Warning in bool convert_clean_string_to_rational(...) : input string ill-defined. Return false." << endl;
    return false;
  }
  int den = atoi(tmp_string.c_str());
  if (den == 0)
  {
    cout << "\n  Warning in bool convert_clean_string_to_rational(...) : input string ill-defined. Return false." << endl;
    return false;
  }

  output = rational<int>(num,den);
  return true;
}



bool convert_string_to_vector_of_rational(string input, rationalVector &output)
{
  output.clear();

  clean_string(input);
  if (input == "")
    return true;

  string::size_type pos1     = 0;
  string::size_type pos2     = 0;
  string::size_type pos3     = 0;
  string::size_type last_pos = 0;

  // input string must be either of the form 
  // "n/d(a, b, c, ...)" where a,b,c,n,d are integers
  // "(n1/d1, a, n2/d2, ...)"
  // "n1/d1, a, n2/d2, ..."
  // or without ","
  pos1 = input.find("/", 0);
  pos2 = input.find("(", 0);

  unsigned case_nr = 0;
  if (pos2 == string::npos)
    case_nr = 3;
  else
  {
    pos3 = input.find(")", pos2);
    if (pos3 != string::npos)
    {
      if (pos1 < pos2)
        case_nr = 1;
      else
        case_nr = 2;
    }
    else
    {
      cout << "\n  Warning in bool convert_string_to_vector_of_rational(...) : input string ill-defined. Return false." << endl;
      return false;
    }
  }
  
  rational<int> rat = 0;
  unsigned exp = 0;
  unsigned i = 0;
  size_t s1 = 0;
  
  string tmp_string1 = "";
  
  if (case_nr == 1)
  {
    // begin: read the common factor
    if (!convert_clean_string_to_rational(input.substr(0, input.find_first_not_of("+-0123456789/", 0)), rat))
    {
      cout << "\n  Warning in bool convert_string_to_vector_of_rational(...) : input string ill-defined. Return false." << endl;
      return false;
    }
    // end: read the common factor
    
    input = input.substr(pos2 + 1, pos3 - pos2 - 1);
     
    vector<int> tmpVector;
    const bool b = convert_clean_string_to_vector_of_int(input, tmpVector);

    s1 = tmpVector.size();
    for (i = 0; i < s1; ++i)
      output.push_back(rat * tmpVector[i]);
    
    return b;
  }

  if ((case_nr == 2) || (case_nr == 3))
  {
    if (case_nr == 2)
      input = input.substr(pos2 + 1, pos3 - pos2 - 1);

    pos1 = 0;
    while (pos1 != string::npos)
    {
      pos1 = input.find(" ", last_pos);

      tmp_string1 = input.substr(last_pos, pos1 - last_pos);

      // begin: read the rational number
      if (!convert_clean_string_to_rational(tmp_string1.substr(0, tmp_string1.find_first_not_of("+-0123456789/", 0)), rat))
      {
        cout << "\n  Warning in bool convert_string_to_vector_of_rational(...) : input string ill-defined. Return false." << endl;
        return false;
      }
      // end: read the rational number

      pos2 = tmp_string1.find("^", 0);
      if (pos2 == string::npos)
        output.push_back(rat);
      else
      {
        // begin: read the exponent
        tmp_string1 = tmp_string1.substr(pos2+1, string::npos);
        if ((tmp_string1 == "") || (tmp_string1.find_first_not_of("0123456789", 0) != string::npos))
        {
          cout << "\n  Warning in bool convert_string_to_vector_of_rational(...) : input string ill-defined. Return false." << endl;
          return false;
        }
        exp = (unsigned)atoi(tmp_string1.c_str());
        // end: read the exponent

        for (i = 0; i < exp; ++i)
          output.push_back(rat);
      }
      last_pos = pos1 + 1;
    }
    return true;
  }

  cout << "\n  Warning in bool convert_string_to_vector_of_rational(...) : case not known. Return false." << endl;
  return false;
}






bool Find_Basis_Of_Orthogonal_Space(const vector<doubleVector> &OriginalSpace, const SelfDualLattice &GaugeLattice, unsigned d, vector<doubleVector> &OrthogonalSpace)
{
  if (OrthogonalSpace.size() != 0)
  {
    cout << "\n  Warning in bool Find_Basis_Of_Orthogonal_Space(...) : basis of orthogonal space not empty. Now cleared." << endl;
    OrthogonalSpace.clear();
  }

  const size_t s1 = OriginalSpace.size();

  unsigned i = 0;

  // gauge group is U(1)^n
  if (s1 == 0)
  {
    doubleVector tmp_vector(d, 0.0);
    for (i = 0; i < d; ++i)
    {
      tmp_vector.assign(d, 0.0);
      tmp_vector[i] = 2.0;
      OrthogonalSpace.push_back(tmp_vector);
    }
    return true;
  }

  if (d != OriginalSpace[0].size())
  {
    cout << "\n  Warning in bool Find_Basis_Of_Orthogonal_Space(...) : dim. of vectors unclear. Return false." << endl;
    return false;
  }

  // begin: convert to CHugeInt
  vector<rational<CHugeInt> > tmp_Vector(d, CHugeInt(0));
  vector<vector<rational<CHugeInt> > > OriginalSpace_CHugeInt(s1, tmp_Vector);

  unsigned j = 0;
  for (i = 0; i < s1; ++i)
  {
    const doubleVector          &ith_OriginalBasisVector = OriginalSpace[i];
    vector<rational<CHugeInt> > &Vector_CHugeInt         = OriginalSpace_CHugeInt[i];

    for (j = 0; j < d; ++j)
      Vector_CHugeInt[j] = D2HugeInt(ith_OriginalBasisVector[j]);
  }
  // end: convert to CHugeInt

  // begin: find orthogonal basis
  CLinAlg<CHugeInt> LA;

  /*cout << "\nOriginalSpace_CHugeInt" << endl;
  for (i = 0; i < OriginalSpace_CHugeInt.size(); ++i)
  {
    for (j = 0; j < OriginalSpace_CHugeInt[i].size(); ++j)
      cout << OriginalSpace_CHugeInt[i][j] << " ";
    cout << endl;
  }
  cout << endl;*/

  vector<vector<rational<CHugeInt> > > OrthogonalSpace_CHugeInt;
  vector<vector<rational<CHugeInt> > > OrthogonalBasis_CHugeInt;

  LA.FindKernel(OriginalSpace_CHugeInt, OrthogonalSpace_CHugeInt, true);
  LA.GramSchmidt(OrthogonalSpace_CHugeInt, OrthogonalBasis_CHugeInt);

  const size_t s2 = OrthogonalBasis_CHugeInt.size();
  for (i = 0; i < s2; ++i)
    LA.ScaleToIntegerVector(OrthogonalBasis_CHugeInt[i]);
  // end: find orthogonal basis

  // begin: convert to double
  doubleVector tmp_vector(d, 0.0);
  OrthogonalSpace.assign(s2, tmp_vector);

  for (i = 0; i < s2; ++i)
  {
    const vector<rational<CHugeInt> > &Vector_CHugeInt = OrthogonalBasis_CHugeInt[i];
    doubleVector                      &BasisVector     = OrthogonalSpace[i];

    for (j = 0; j < d; ++j)
    {
      const rational<CHugeInt> &entry = Vector_CHugeInt[j];
      BasisVector[j] = (double)entry.numerator().ToLongLongInt()/(double)entry.denominator().ToLongLongInt();
    }
  }
  // end: convert to double

  long double sp = 0.0;

  // begin: bring the vectors of "OrthogonalSpace" to lattice vectors
  if ((GaugeLattice != UNSPECIFIED_LATTICE) && (d == 16))
  {
    bool stop = false;
    CLatticeVector LatticeVector(d);

    for (i = 0; i < s2; ++i)
    {
      doubleVector &BasisVector = OrthogonalSpace[i];

      stop = false;

      for (double k = 1.0; !stop && (k < 50); ++k)
      {
        for (j = 0; j < d; ++j)
          LatticeVector[j] = k * BasisVector[j];

        if (GaugeLattice == E8xE8)
        {
          if (LatticeVector.From_E8_Lattice(1) && LatticeVector.From_E8_Lattice(2))
            stop = true;
        }
        else
        {
          if (LatticeVector.From_Spin32_Lattice())
            stop = true;
        }
      }
      if (stop == false)
        cout << "\n  Warning in bool Find_Basis_Of_Orthogonal_Space(...) : basis vector can not be scaled to the lattice." << endl;
      BasisVector = LatticeVector;
    }
  }
  // end: bring the vectors of "OrthogonalSpace" to lattice vectors
    
  if (s2 + s1 != d)
  {
    cout << "\n  Warning in bool Find_Basis_Of_Orthogonal_Space(...) : Return false." << endl;
    cout << "    dim. of orig. space = " << s1 << endl;
    cout << "    dim. of orth. soace = " << s2 << endl;
    cout << "    dim. of space       = " << d << endl;

    for (i = 0; i < s1; ++i)
      cout << OriginalSpace[i] << endl;
    return false;
  }

  unsigned k = 0;

  // begin: check that the vectors of "OrthogonalSpace" are orthogonal to the "OriginalSpace"
  for (i = 0; i < s2; ++i)
  {
    const doubleVector &BasisVector = OrthogonalSpace[i];

    for (j = 0; j < s1; ++j)
    {
      const doubleVector &OriginalVector = OriginalSpace[j];

      // sp = BasisVector x SimpleRoot
      sp = 0.0;

      for (k = 0; k < d; ++k)
        sp += (BasisVector[k] * OriginalVector[k]);

      if (fabs(sp) > 0.0001)
      {
        cout << "\n  Warning in bool Find_Basis_Of_Orthogonal_Space(...) : new vector is not orthogonal to the space. Return false." << endl;
        return false;
      }
    }
  }
  // end: check that the vectors of "OrthogonalSpace" are orthogonal to the "OriginalSpace"
  return true;
}


bool operator==(const FPCoordinates &FP1, const FPCoordinates &FP2)
{
  const double prec = 0.0001;
  for (unsigned i = 0; i < ComplexLatticeDim; ++i)
  {
    if (FP1.FixedTorus[i] != FP2.FixedTorus[i])
      return false;

    if (abs(FP1.Coordinates[i] - FP2.Coordinates[i]) > prec)
      return false;
  }
  return true;
}










