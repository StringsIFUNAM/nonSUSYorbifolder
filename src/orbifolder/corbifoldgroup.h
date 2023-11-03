
#ifndef CORBIFOLDGROUP_H
#define CORBIFOLDGROUP_H

#include <fstream>
#include <iostream>

#include "cspacegroup.h"
#include "cvector.h"
#include "corbifoldgroupelement.h"
#include "cwilsonline.h"
#include "cwilsonlines.h"
#include "cshiftvector.h"


using std::vector;
using std::string;

class CPrint;
class CRandomModel;

#ifndef STRUCT_GDT_PARAMETERS
#define STRUCT_GDT_PARAMETERS
struct GDT_Parameters
{
  bool UseGDT;
  rational<int> a;
  rationalVector b,c;
  rationalMatrix d;
};
#endif

class COrbifoldGroup {
public:
// member functions
  COrbifoldGroup();
  COrbifoldGroup(const CSpaceGroup &SpaceGroup);
  COrbifoldGroup(const CSpaceGroup &SpaceGroup, const vector<CShiftVector> &Shifts, const CWilsonLines &Wilsonlines);
  ~COrbifoldGroup();

  bool                                          LoadOrbifoldGroup(std::ifstream &in, string &ProgramFilename);
  bool                                          CreateRandom(const CRandomModel &RandomModel, bool CreateBrotherModel);
  bool                                          CreateRandom_shifts(const CRandomModel &RandomModel, bool CreateBrotherModel);
  bool                                          CreateRandom_internal(const CRandomModel &RandomModel, bool CreateBrotherModel);
  bool                                          CreateSubGroup(CPrint &Print, vector<unsigned> &TwistFactors1, vector<unsigned> &TwistFactors2, bool use_WL_in_FixedTori, COrbifoldGroup &SubOrbifoldGroup) const;

  bool                                          Reset(bool CreateModelIndependentPart, bool ResetShifts, bool ResetWilsonLines, bool ResetTorsion);
  bool                                          CreateModelDependentPart(CPrint &Print, bool PrintWarningModularInvariance = true);

  bool                                          GetShiftVector(const CSpaceGroupElement &constructing_Element, CShiftVector &result) const;
  bool                                          GetTwistVector(const CSpaceGroupElement &constructing_Element, CTwistVector &result) const;

  void                                          LocalGaugeGroupOnFixedTorus(std::ostream &out) const;
  void                                          PrintToFile(std::ostream &out) const;

  void                                          SetDiscreteTorsionToZero();
  bool                                          SetDiscreteTorsion(const GDT_Parameters &DiscreteTorsion);

  CShiftVector                                 &AccessShift(const unsigned &i);
  CWilsonLines                                 &AccessWilsonLines();
  CSpaceGroup                                  &AccessSpaceGroup();

  SelfDualLattice                               GetLattice() const {return this->Shifts[0].GetLattice();};
  const COrbifoldGroupElement                  &GetElement(const unsigned &i) const;
  const vector<COrbifoldGroupElement>          &GetElements() const;
  const vector<COrbifoldGroupElement>          &GetCentralizer(const unsigned &i) const;
  const CShiftVector                           &GetShift(const unsigned &i) const;
  const vector<CShiftVector>                   &GetShifts() const;
  const CWilsonLines                           &GetWilsonLines() const;
  const GDT_Parameters                         &GetDiscreteTorsion() const;
  const unsigned                               &GetOrderZM() const {return this->SpaceGroup.GetM();};
  const unsigned                               &GetOrderZN() const {return this->SpaceGroup.GetN();};
  
  const unsigned                               &GetOrderZK() const {return this->SpaceGroup.GetK();};  //I added on June23 night for ZMxZNxZK
  
  
  const CSpaceGroup                            &GetSpaceGroup() const;
  const vector<CVector>                        &GetInvariantSupercharges() const;
  const int                                    &GetNumberOfSupersymmetry() const;

  const CheckStatus                            &GetSpaceGroup_CheckStatus() const {return this->SpaceGroup.GetCheckStatus();};
  const CheckStatus                            &GetModelIndependent_CheckStatus() const {return this->ModelIndependent_CheckStatus;};
  const CheckStatus                            &GetModularInvariance_CheckStatus() const {return this->ModularInvariance_CheckStatus;};
  const CheckStatus                            &GetOrbifoldGroup_CheckStatus() const {return this->OrbifoldGroup_CheckStatus;};
  const CheckStatus                            &GetWilsonLines_CheckStatus() const {return this->WilsonLines.GetCheckStatus();};

// member variables
  CWilsonLine                                   FreelyActingWilsonLine;
  bool                                          UseFreelyActingWL;
  string                                        Label;
  vector<CVector>                               LoadedU1Generators;

private:
// member functions
  bool                                          CreateModelIndependentPart();
  bool                                          CheckModularInvariance(CPrint &Print, bool info);

// member variables
  vector<COrbifoldGroupElement>                 Elements;
  vector<vector<COrbifoldGroupElement> >        Centralizer;

  vector<CShiftVector>                          Shifts;
  CWilsonLines                                  WilsonLines;
  GDT_Parameters                                DiscreteTorsion;

  CSpaceGroup                                   SpaceGroup;
  vector<CVector>                               InvariantSupercharges;
  int                                           NumberOfSupersymmetry;

  CheckStatus                                   ModelIndependent_CheckStatus;
  CheckStatus                                   ModularInvariance_CheckStatus;
  CheckStatus                                   OrbifoldGroup_CheckStatus;
};

#endif
