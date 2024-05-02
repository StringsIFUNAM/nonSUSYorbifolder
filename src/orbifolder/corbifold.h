
#ifndef CORBIFOLD_H
#define CORBIFOLD_H
#include <vector>
#include <string>

#include "corbifoldgroup.h"
#include "cfixedpoint.h"
#include "csector.h"
#include "corbifoldcore.h"
#include "canalysemodel.h"
#include "cfield.h"
#include "cgaugeinvariance.h"

using std::vector;
using std::string;

class OrbifoldGroup;
class CFixedBrane;
class CSector;

#ifndef STRUCT_YUKAWACOUPLING
#define STRUCT_YUKAWACOUPLING
struct YukawaCoupling
{
  vector<unsigned>          FieldIndices;
  vector<vector<unsigned> > ComponentCoupling;
  double                    CouplingStrength;
  
  vector<string>            LabelOfModuli;
  vector<int>               ExponentsOfEtas;
};
#endif

#ifndef TYPEDEF_CGAUGEGROUP
#define TYPEDEF_CGAUGEGROUP
typedef gaugeGroup<double> CGaugeGroup;
#endif

#ifndef STRUCT_PID
#define STRUCT_PID
struct PID
{
  vector<pid_t>             PIDs;
  vector<string>            PID_Commands;
  vector<unsigned>          PID_JobIndices;
  vector<time_t>            PID_StartingTimes;
  vector<bool>              PID_Done;
  vector<string>            PID_Filenames;
};
#endif

#ifndef STRUCT_SSYMMETRYGROUP
#define STRUCT_SSYMMETRYGROUP
struct SSymmetryGroup
{
  CGaugeGroup      GaugeGroup;

  vector<unsigned> observable_sector_GGs;
  vector<unsigned> observable_sector_U1s;
  vector<string>   GGs_AdditionalLabels;
  vector<string>   U1s_AdditionalLabels;

  unsigned         Position_of_and_in_GaugeGroup;
  bool             IsFirstU1Anomalous;
  double           D0_FI_term;

  CVector          BmLGenerator;
};
#endif

#ifndef STRUCT_SCONFIG
#define STRUCT_SCONFIG
struct SConfig
{
  string   ConfigLabel;
  unsigned ConfigNumber;

  vector<CVector> InvariantSupercharges;

  SSymmetryGroup                  SymmetryGroup;
  vector<CField>                  Fields;

  vector<YukawaCoupling>          FieldCouplings;

  vector<string>                  NamesOfSetsOfFields;
  vector<vector<unsigned> >       SetsOfFields;

  vector<vector<YukawaCoupling> > FTerms;

  unsigned                        use_Labels;

  unsigned						  HiggsNo;					//hacking here!!!

  PID pid;
};
#endif

class COrbifold{
public: 
  COrbifold();
  COrbifold(const COrbifoldGroup &OrbifoldGroup);
  COrbifold(const COrbifoldGroup &OrbifoldGroup, const vector<CSector> &CoreSpectrum);

  ~COrbifold();

  bool                   Reset(bool CreateModelIndependentPart, bool ResetShifts, bool ResetWilsonLines, bool ResetTorsion);
  bool                   Create();

  bool                   CheckAnomaly(SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, CPrint &Print, bool info = false, double AddFactor = 1.0) const;
  bool                   CheckDiscreteAnomaly(const SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, CSpaceGroupElement &anomalous_element, CPrint &Print, bool info) const;

  bool                   FindGaugeGroup(const CFixedBrane &UntwistedSector, SConfig &VEVConfig, bool CreateAnomalousU1Generator = true);

  bool                   Config_SetU1Direction(SConfig &NewU1VEVConfig, const CVector &U1Direction, unsigned pos_of_U1 = 0) const;
  bool                   Config_ComputeNewU1Charges(SConfig &NewU1VEVConfig) const;

  size_t                 GetNumberOfSectors() const { return this->Sectors.size();};

  CSector               &AccessSector(const unsigned &i);
  CFixedBrane           &AccessFixedBrane(const CSpaceGroupElement &Element, bool &FixedBraneFound);

  const vector<CSector> &GetSectors() const;
  const CSector         &GetSector(const unsigned &i) const;
  const CheckStatus     &GetCheckStatus() const {return this->Orbifold_CheckStatus;};

  /*
  bool CreateLocalFixedBrane(const SpaceGroup_Element &Element_Label, bool WilsonLines_In_FixedTori, bool PrintDetails = false, bool printUSector = false, const CGaugeGroup *GaugeGroup4D = NULL, int pos_of_Y = -1);
  void DecomposeLocalFixedBrane(const SpaceGroup_Element &Element_Label);
*/

  COrbifoldGroup         OrbifoldGroup;

  SConfig                StandardConfig;
  SConfig                TachyonicStandardConfig;

private :
  bool                   Config_Clear(SConfig &VEVConfig, string &Label);
  bool                   CheckCPPartner() const;
  bool                   CreateRepresentations();
  bool                   TachyonicCreateRepresentations();

  vector<CSector>        Sectors;

  CheckStatus            Orbifold_CheckStatus;
};

#endif

