
#ifndef CSPACEGROUP_H
#define CSPACEGROUP_H

#include <vector>
#include <complex>

#include "ctwist.h"
#include "cspacegroupelement.h"

using std::vector;
using std::complex;
using std::string;

class CPrint;


#ifndef TYPEDEF_COMPLEXVECTOR
#define TYPEDEF_COMPLEXVECTOR
typedef vector<complex<double> > complexVector;
#endif

#ifndef TYPEDEF_COMPLEXMATRIX
#define TYPEDEF_COMPLEXMATRIX
typedef vector<complexVector> complexMatrix;
#endif

#ifndef STRCUT_FPCOORDINATES
#define STRCUT_FPCOORDINATES
struct FPCoordinates
{
  complexVector Coordinates;
  vector<bool>  FixedTorus;
};
#endif

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

#ifndef STRUCT_SDISCRETESYMMETRY
#define STRUCT_SDISCRETESYMMETRY
struct SDiscreteSymmetry
{
  unsigned           Order;
  rational<int>      SuperpotentialCharge;
  rationalVector     ChargeOperator;
  CSpaceGroupElement SG_SymmetryGenerator;
};
#endif

#ifndef STRUCT_SMODULARSYMMETRY
#define STRUCT_SMODULARSYMMETRY
struct SModularSymmetry
{
  string         Label;
  rational<int>  Const;
  unsigned       Index;
  rationalVector ChargeOperator; // 3 components: Index-component of q_sh, N and \bar{N}
};
#endif

class CSpaceGroup {
public: 
  CSpaceGroup();
  CSpaceGroup(const unsigned &M, const unsigned &N, const unsigned &K);
  ~CSpaceGroup();

  bool                                       Check();
  void                                       SG_NotChecked() {this->SpaceGroup_CheckStatus = NotChecked;};

  void                                       PrintPointGroup(std::ostream &out) const;
  bool                                       PrintGeometryFile(std::ostream &out) const;

  void                                       Clear();
  bool                                       SetOrder(const unsigned &M, const unsigned &N, const unsigned &K);
  bool                                       CreateSubGroup(CPrint &Print, vector<unsigned> &TwistFactors1, vector<unsigned> &TwistFactors2, bool use_WL_in_FixedTori, CSpaceGroup &SubGroup) const;
  bool                                       LoadSpaceGroup(string ifilename, bool reload = false);
  void                                       SetSpaceGroup(const CSpaceGroup &SpaceGroup);

  bool                                       SG_Commute(const CSpaceGroupElement &SGElement1, const CSpaceGroupElement &SGElement2) const;
  bool                                       SG_FixedPoint(const CSpaceGroupElement &SGElement) const;
  bool                                       SG_FixedPoint(const CSpaceGroupElement &SGElement, FPCoordinates &result) const;
  bool                                       SG_Multiply(const CSpaceGroupElement &SGElement1, const CSpaceGroupElement &SGElement2, CSpaceGroupElement &result) const;
  bool                                       SG_Multiply(const CSpaceGroupElement &SGElement, const complexVector &LatticeVector, complexVector &result) const;
  bool                                       SG_Multiply(const CSpaceGroupElement &SGElement, const FPCoordinates &FixedPoint, FPCoordinates &result) const;
  bool                                       SG_DifferenceInLattice(const FPCoordinates &fp1, const FPCoordinates &fp2) const;
  bool                                       SG_FromSameConjugationClass(const CSpaceGroupElement &Element1, const CSpaceGroupElement &Element2) const;
  bool                                       SG_ReverseVector(const CSpaceGroupElement &SGElement, complexVector &ReverseVector) const;
  CSpaceGroupElement                         SG_Inverse(const CSpaceGroupElement &SGElement) const;

  const string                              &GetAdditional_Label() const {return this->additional_label;};
  size_t                                     GetNumberOfElementsOfSector(unsigned sector) const { return this->Sectors.at(sector).size();};
  size_t                                     GetNumberOfSectors() const { return this->Sectors.size();};
  const CSpaceGroupElement                  &GetElement(unsigned sector, unsigned i) const { return this->Sectors.at(sector).at(i);};
  const vector<unsigned>                    &GetWL_AllowedOrders() const { return this->WL_AllowedOrders;};
  const vector<vector<unsigned> >           &GetWL_Relations() const { return this->WL_Relations;};
  const unsigned                            &GetM() const { return this->M; }
  const unsigned                            &GetN() const { return this->N; }
  const unsigned                            &GetK() const { return this->K; }
  const complexMatrix                       &GetT6Lattice() const { return this->T6_Lattice;}
  const bool                                &IsZMxZN() const { return this->ZMxZN; }
  const bool                                &IsZMxZNxZK() const { return this->ZMxZNxZK; }
  const CTwistVector                        &GetTwist(const unsigned &i) const {return this->Twists[i];};
  const vector<CTwistVector>                &GetTwists() const {return this->Twists;};
  const vector<vector<CSpaceGroupElement> > &GetCentralizer() const {return this->Centralizer;};
  const vector<vector<CSpaceGroupElement> > &GetSectors() const  {return this->Sectors;};

  const vector<CSpaceGroupElement>          &GetSG_Generators_Twist() const {return this->SG_Generators_Twist;};
  const vector<CSpaceGroupElement>          &GetSG_Generators_Shift() const {return this->SG_Generators_Shift;};
  const vector<CSpaceGroupElement>          &GetSG_AllNonStandardShifts() const {return this->SG_AllNonStandardShifts;};

  const CheckStatus                         &GetCheckStatus() const {return this->SpaceGroup_CheckStatus;};

  string                                     additional_label;
  string                                     GeometryFilename;
  string                                     lattice_label;

  vector<vector<double> >                    ShiftsWL_ScalarProductFactors;

  vector<SDiscreteSymmetry>                  DiscreteNonRSymmetries;
  vector<SDiscreteSymmetry>                  DiscreteRSymmetries;
  vector<SModularSymmetry>                   ModularSymmetries;

private:
  bool                                       CreateDualTorus();
  bool                                       CreateTwistMatrices();
  bool                                       CreateSetOfElements();
  bool                                       CreateConstructingElements();
  bool                                       CreateCentralizerElements();
  bool                                       FindRelationsOfWilsonLines();
  bool                                       CreateElementsForConjugacyClass();

  vector<CTwistVector>                       Twists;

  vector<vector<rationalMatrix> >            TwistMatrices;
  complexMatrix                              T6_Lattice;
  vector<CVector>                            DualT6_Lattice;

  vector<CSpaceGroupElement>                 SG_Generators_Twist;
  vector<CSpaceGroupElement>                 SG_Generators_Shift;

  //Sectors contains all constructing elements sorted by sectors
  vector<vector<CSpaceGroupElement> >        Sectors;
  vector<vector<CSpaceGroupElement> >        Centralizer;

  vector<CSpaceGroupElement>                 SG_SetOfElements;

  vector<CSpaceGroupElement>                 SG_AllRotations;
  vector<CSpaceGroupElement>                 SG_AllNonStandardShifts;

  vector<vector<unsigned> >                  WL_Relations;
  vector<unsigned>                           WL_AllowedOrders;

  // Z_M x Z_N orbifold
  bool ZMxZN;
  bool ZMxZNxZK;							//hacking here!!!

  unsigned M;
  unsigned N;
  unsigned K;

  CheckStatus                                SpaceGroup_CheckStatus;
};

#endif
