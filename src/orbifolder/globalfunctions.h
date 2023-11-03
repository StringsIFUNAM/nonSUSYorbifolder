
#ifndef CGLOBALFUNCTION_H
#define CGLOBALFUNCTION_H

#include <boost/config.hpp>
#include <boost/rational.hpp>
#include <boost/tokenizer.hpp>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <sstream>
#include <cmath>
#include "chugeint.h"
#include "groupTheory.hpp"

#include <iostream>
#include <fstream>
using namespace std;

using std::vector;
using boost::rational;
using std::string;

struct SDimension;
struct SConfig;
struct SDiscreteSymmetry;
struct SSymmetryGroup;

class CVector;
class CField;
class CSpaceGroupElement;

#ifndef STRUCT_SDIMENSION
#define STRUCT_SDIMENSION
struct SDimension
{
  int    Dimension;
  string AdditionalLabel;
};
#endif

#ifndef STRCUT_FPCOORDINATES
#define STRCUT_FPCOORDINATES
struct FPCoordinates
{
  complexVector Coordinates;
  vector<bool>  FixedTorus;
};
#endif

typedef vector< int > intVector;
typedef vector< intVector > intMatrix;
typedef vector< double > doubleVector;
typedef vector< doubleVector > doubleMatrix;
typedef vector<SDimension> RepVector;

#ifndef ENUM_SELFDUALLATTICE
#define ENUM_SELFDUALLATTICE
enum SelfDualLattice {UNSPECIFIED_LATTICE = 0, E8xE8, Spin32, SO8};
#endif



bool                   clean_string(string &input);
bool                   convert_clean_string_to_vector_of_int(string input, vector<int> &output);
bool                   convert_string_to_vector_of_int(string input, vector<int> &output);
bool                   convert_clean_string_to_rational(string input, rational<int> &output);
bool                   convert_string_to_vector_of_rational(string input, rationalVector &output);
bool                   global_LoadSGElement(const string &input, CSpaceGroupElement &SGElement);

bool                   global_ReplaceString(string &input, const string &replace_from, const string &replace_to);
void                   global_DecomposeString(const string &input, vector<string> &output, const string separator);
bool                   global_ConvertVectorOfString2VectorOfUnsigned(const vector<string> &input, vector<unsigned> &output);

vector<unsigned>       global_GetIndices(const vector<string> &FieldLabels, unsigned use_Labels, const vector<CField> &Fields);
vector<unsigned>       global_GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels, unsigned use_Labels, const vector<CField> &Fields);

istream&               GetSaveLine(istream &in, string &result);

bool                   IsSinglet(const RepVector &Dimensions);
bool                   CC(const RepVector &Dimensions, RepVector &DimensionsCC);
bool                   AreHighestWeights_CC(const vector<int> &HighestWeight_DL_1, const vector<int> &HighestWeight_DL_2, const gaugeGroupFactor<double> &GGFactor);
bool                   AreRepVectorsEqual(const SSymmetryGroup &SymmetryGroup, const RepVector &Dimensions1, const RepVector &Dimensions2);
bool                   AreRepVectorsEqual(const RepVector &Dimensions1, const RepVector &Dimensions2);
bool                   operator==(const RepVector &Dimensions1, const RepVector &Dimensions2);
bool                   AreU1ChargesEqual(const SSymmetryGroup &SymmetryGroup, const CVector &U1Charges1, const CVector &U1Charges2);
unsigned               Multiplicity(const SSymmetryGroup &SymmetryGroup, const RepVector &Dimensions);

bool                   operator==(const FPCoordinates &FP1, const FPCoordinates &FP2);

unsigned               kgV(unsigned m, unsigned n);
unsigned               ggT(unsigned m, unsigned n);

rational<int>          D2Rat(const double x);
rational<CHugeInt>     D2HugeInt(const double x);

void round_double_to_int(const double &d, int &i);
int  round_double_to_int(const double &d);
bool is_integer(const double &d);
bool is_even(const double &d);
void RoundDouble(double &tmp);
bool RoundCharge(long double &U1_Charge);

void RecursiveCounting(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, const vector<vector<unsigned> > &vec_adjoining_Positions, vector<vector<unsigned> > &AllNumbers);
void RecursiveCounting(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<vector<unsigned> > &AllNumbers);
void RecursiveCounting_to_files(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<std::ofstream*> &outs, unsigned long &counter);

bool NextNumberNoTwins(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const vector<vector<unsigned> > &PositionsOfUnequalDigits, const unsigned &OnlyChangeUpToDigit);
bool NextNumber(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const size_t &number_of_Digits, const unsigned &OnlyChangeUpToDigit);
bool NextNumber(vector<unsigned> &Number, const vector<unsigned> &MaxDigits, const size_t &number_of_Digits);

bool Find_Basis_Of_Orthogonal_Space(const vector<doubleVector> &OriginalSpace, const SelfDualLattice &GaugeLattice, unsigned d, vector<doubleVector> &OrthogonalSpace);

#endif
