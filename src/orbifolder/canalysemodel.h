#ifndef CANALYSEMODEL_H
#define CANALYSEMODEL_H

#include "cvector.h"
#include <string>
#include "cspacegroup.h"
#include "cstate.h"
#include "cgaugeinvariance.h"

using std::string;
using std::istream;
using std::ostream;
using std::ifstream;

class CHugeInt;
class CMonomial;
class COrbifold;
class CYukawaCouplings;
class CPrint;

struct SConfig;
struct SDimension;
struct YukawaCoupling;

typedef vector<SDimension> RepVector;

struct PhenoScheme
{
  string            SchemeLabel;
  vector<string>    GaugeGroupFactors;
  vector<string>    GGs_AdditionalLabels;
  vector<RepVector> SetOfDimensions;
  vector<CVector>   SetOfU1Charges;
  vector<string>    Labels;
  vector<SUSYMultiplet> Multiplets; 

  // i specifies which dimension and U(1) charges specify the normalization
  unsigned          IndexOfSetToNormalizeCharges; // equal i
  
  vector<double>    NormalizationOfLengthOfU1Generators;
};


class CAnalyseModel{
public:
  CAnalyseModel();
  ~CAnalyseModel();

  bool AnalyseModel(const COrbifold &Orbifold, const SConfig &OriginalVEVConfig, bool &SM, bool &PS, bool &SU5, vector<SConfig> &AllVEVConfigs, CPrint &Print, unsigned NumberOfGenerations = 3, bool SM_PrintSU5SimpleRoots = false) const;
  
  bool ComputeN2BetaFunctionCoefficient(const COrbifold *Orbifold, const SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, vector<bool> &N2_Sector_Exists, vector<vector<rational<int> > > &BetaFunctionCoefficients) const;
  rational<int> ComputeBetaFunctionCoefficient(const CGaugeIndices &GaugeIndices, const SConfig &VEVConfig, const unsigned factor) const;

  bool GroupAFromB(const gaugeGroupFactor<double> &ggfA, const gaugeGroupFactor<double> &ggfB) const;
    
  bool CreatePhenoScheme(const string &SchemeLabel, PhenoScheme &Scheme) const;
  bool FindPositionOfGaugeGroup(const PhenoScheme &Scheme, const SConfig &OriginalVEVConfig, vector<SConfig> &GoodVEVConfigs) const;
  bool AutoCreateLabels(const PhenoScheme &Scheme, SConfig &VEVConfig) const;

  bool Labels_Create(istream &in, SConfig &VEVConfig, CPrint &Print, const vector<SUSYMultiplet>  &Multiplet, bool info = true) const; 
  bool Labels_Load(string ifilename, SConfig &VEVConfig) const;
  
  bool FindCouplings(const SConfig &VEVConfig, const vector<string> &FieldLabels, vector<YukawaCoupling> &Result, bool EffectiveCoupling = false) const;
  
  bool SM_FindPositionOfGaugeGroupFromGG(const SConfig &OriginalVEVConfig, const vector<vector<double> > &GG_SimpleRoots, vector<SConfig> &SMVEVConfigs) const;
  int  SM_NetNumberOfGenerations(const SConfig &SMVEVConfig) const;
  bool SM_GetHypercharges(const SConfig &SMVEVConfig, vector<CVector> &Hypercharges, vector<vector<CVector> > &PossibleRootsOfSU5, const vector<CVector> &Roots10D) const;
  bool SM_CheckVectorlikeness(const SConfig &SMVEVConfig, unsigned NumberOfGenerations, bool normalizeY, bool info, unsigned &HiggsNo) const;
  bool SM_GetProtonHexalityFromSO12(const COrbifold &Orbifold, const SConfig &SMVEVConfig, const vector<CVector> &Roots10D, vector<CVector> &ProtonHexalities) const;

  int  PS_NetNumberOfGenerations(const SConfig &PSVEVConfig) const;
  bool PS_CheckVectorlikeness(SConfig &PSVEVConfig, CPrint &Print, unsigned NumberOfGenerations = 3, bool info = false) const;
  bool PS_GetProtonHexality(const SConfig &PSVEVConfig, const vector<CVector> &Roots10D, vector<CVector> &ProtonHexalities);

  int  SU5_NetNumberOfGenerations(const SConfig &SU5VEVConfig) const;
  bool SU5_GetFlippedU1(const SConfig &SU5VEVConfig, const vector<CVector> &Roots10D, vector<CVector> &FlippedU1s) const;
  bool SU5_CheckVectorlikeness(const COrbifold &Orbifold, SConfig &SU5VEVConfig, unsigned pos_of_flippedU1, CPrint &Print, bool info = false) const;

  bool AreSomeFixedBranesEmpty(const COrbifold &Orbifold, const SConfig &VEVConfig) const;
  void FindEmptyFixedBranes(const COrbifold &Orbifold, const SConfig &VEVConfig, vector<CSpaceGroupElement> &EmptyFixedBranes) const;
  bool GetAllBlowUpModePerFixedPoint(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, vector<vector<unsigned> > &BlowupModes) const;
  bool GetOneBlowUpModePerFixedPoint(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, vector<unsigned> &BlowupModes) const;

  bool AccidentalU1Charges_Find(SConfig &VEVConfig, const vector<unsigned> &FieldIndicesWithZeroAccCharge, CPrint &Print, bool info = false) const;
  bool AccidentalU1Charges_Load(SConfig &VEVConfig, ifstream &in);
  bool AccidentalU1Charges_Save(const SConfig &VEVConfig, ostream &out);
  
  void FTerms_Compute(SConfig &VEVConfig, unsigned upto_Order_inW) const;
  bool DTerms_FindDMonomials(const COrbifold &Orbifold, SConfig &VEVConfig, const vector<unsigned> &FieldIndices, bool D0_withFI, vector<vector<unsigned> > &GaugeEquivalentFields, vector<unsigned> &AllFieldsInMonomial) const;

  bool SetBmLGenerator(const COrbifold &Orbifold,SConfig &VEVConfig, const CVector &BmLGenerator, bool NeedsToBeNonAnomalous) const;

  bool SetVev(SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const vector<unsigned> &FieldComponents, double VEV, bool random, vector<vector<double> > *U1ChargesOfVevFields = NULL);
  bool FindUnbrokenSymmetries(const COrbifold &Orbifold, const SConfig &OriginalVEVConfig, const vector<CVector> &Weights, SConfig &VEVConfig, CPrint &Print, bool info = false);
  bool FindInvariantConfig(const COrbifold &Orbifold, const CVector &Background, const SConfig &OriginalVEVConfig, SConfig &NewVEVConfig);
};

#endif
