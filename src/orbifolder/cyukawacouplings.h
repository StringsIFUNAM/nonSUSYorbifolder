#ifndef CYUKAWACOUPLINGS_H
#define CYUKAWACOUPLINGS_H

#include "cspacegroupelement.h"
#include "cgaugeinvariance.h"
#include "canalysemodel.h"

class COrbifold;
struct Field;
struct YukawaCoupling;


class CYukawaCouplings {
public:
  CYukawaCouplings();
  ~CYukawaCouplings();

  bool                      Initiate(const CSpaceGroup &SpaceGroup, const SConfig &Vacuum);

  unsigned                  Load(const COrbifold &Orbifold, SConfig &Vacuum, std::ifstream &in);
  void                      Save(const COrbifold &Orbifold, const SConfig &Vacuum, std::ostream &out, unsigned Order) const;

  unsigned                  AddCoupling(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &FieldLabels, const vector<string> *AllowedFields = NULL);
  void                      CreatePreCouplings(const COrbifold &Orbifold, const SConfig &Vacuum, const vector<unsigned> &MaxDigits, const vector<unsigned> &vec_adjoining_Positions, const vector<vector<unsigned> > &PoolsOfFields, unsigned PoolIndex, const vector<unsigned> &NumbersOfFieldsFromPool);
  bool                      DeleteCouplingsWithStatesTwice(SConfig &Vacuum, const CPrint &Print, vector<vector<string> > &Known_GG_string, vector<vector<vector<vector<RepVector> > > > &Known_Zero_SortedCoupling, vector<vector<vector<vector<RepVector> > > > &Known_NonZero_SortedCoupling);

  bool                      AutoCreateMassMatrix(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &Labels, unsigned maxOrderOfSinglets, const string &SingletLabel, const vector<string> *AllowedFields = NULL);
  bool                      AutoCreateFullMassMatrix(const COrbifold &Orbifold, SConfig &Vacuum, const vector<string> &Labels, unsigned maxOrderOfSinglets, const string &SingletLabel, const vector<string> *AllowedFields = NULL);

  vector<vector<unsigned> > KnownCouplingsCSector;
  vector<vector<CVector> >  KnownBasis_of_DualLattice;
  vector<vector<CVector> >  KnownBasis_of_Lattice;
  vector<unsigned>          KnownDim_DualLattice;
  vector<bool>              KnownUse_DualBasis;

  std::ostream              *tmp_file;

  // YC = Yukawa coupling
  vector<vector<bool> >                        YC_WhereAreFixedTori;
  vector<vector<vector<CVector> > >            YC_Basis_of_Rotated_Lattices_For_kl_Twisted_Sectors;
  size_t                                       YC_ZMxZN_Twists_size;
  vector<vector<vector<CSpaceGroupElement> > > YC_CCs_of_constructing_Elements;

  CGaugeInvariance                             GI;
};

#endif
