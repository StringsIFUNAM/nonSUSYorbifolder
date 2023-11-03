#ifndef CSECTOR_H
#define CSECTOR_H

#include "cfixedpoint.h"
#include "coscillator.h"
#include "chalfstate.h"

//! CSector.
/*!
A CSector object contains the local twist "Twist", the set of fixed branes "FixedBranes", the oscillator excitations and the right-movers.
 */

using std::vector;

class CFixedBrane;

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

class CSector {
public: 
// member functions
  CSector();
  ~CSector();
  
  bool                                  Create(unsigned m, unsigned n, unsigned k, const CTwistVector &constructing_Twist, const vector<CTwistVector> &ZMxZNxZK_Twists);
  bool                                  AddFixedBrane(const CSpaceGroupElement &SGElement, const unsigned &Index_SGElement, const string &FixedBraneLabel);
  
  bool                                  GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const;
  bool                                  SectorHasLeftchiralRightmover() const;

  size_t                                GetNumberOfFixedBranes() const { return this->FixedBranes.size();};
  size_t                                GetNumberOfMasslessRightMovers() const { return this->MasslessRightMovers.size();};
  size_t                                GetNumberOfRightMovers() const { return this->RightMovers.size();};
  
  CFixedBrane                          &AccessFixedBrane(const unsigned &i);

  const CTwistVector                   &GetTwist() const {return this->Twist;};
  const vector<unsigned>               &Get_mnk() const {return this->Sector;};
  const unsigned                       &Get_m() const {return this->Sector[0];};
  const unsigned                       &Get_n() const {return this->Sector[1];};
  const unsigned                       &Get_k() const {return this->Sector[2];};

  const CFixedBrane                    &GetFixedBrane(const unsigned &i) const;

  const vector<S_OscillatorExcitation> &GetLM_Excitations() const;
  const vector<vector<S_OscillatorExcitation> > &TachyonicGetLM_Excitations() const;
  const vector<CModedOscillator>       &GetLM_all_Oscillators() const;
  const CModedOscillator               &GetLM_Oscillator(const unsigned &i) const;
  const vector<S_OscillatorExcitation> &GetRM_Excitations() const;
  const vector<CModedOscillator>       &GetRM_all_Oscillators() const;
  const CModedOscillator               &GetRM_Oscillator(const unsigned &i) const;
  
  const CMasslessHalfState             &GetMasslessRightMover(const unsigned &i) const;
  const vector<CMasslessHalfState>     &GetMasslessRightMovers() const;
  const CHalfState                     &GetRightMover(const unsigned &i) const;
  const vector<CHalfState>             &GetRightMovers() const;
  const vector<CHalfState>             &TachyonicGetRightMovers() const;
  const vector<CTachyonHalfState>	   &GetRTachyons() const;			//hacking here!!!
  vector<CTwistVector>		   		   Twists;							//hacking here!!!
  
private:
// member functions
  bool                                  CreateOscillators(const vector<CTwistVector> &ZMxZN_Twists);
  bool                                  CreateOscillatorExcitations();
  void                                  RecursiveCreate_LM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index);
  void  	  	  	  	  	  	  	  	RecursiveCreate_tachyonicLM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index, const double &M_R);
  void                                  RecursiveCreate_RM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index);
  bool                                  CreateMasslessRightMover();
  bool                                  CreateTachyonicRightMover();
  bool                                  SortByEigenvalue(const vector<CTwistVector> &ZMxZN_Twists);
  bool                                  TachyonicSortByEigenvalue(const vector<CTwistVector> &ZMxZN_Twists);

// member variables
  CTwistVector                          	Twist;
  vector<unsigned>                      	Sector;
  vector<CFixedBrane>                   	FixedBranes;

  vector<S_OscillatorExcitation>        	LM_Excitations;
  vector<S_OscillatorExcitation>	    	TachyonicLM_Excitations_tmp;
  vector<vector<S_OscillatorExcitation> >	TachyonicLM_Excitations;		//first vector index keeps track of tachyonic mass level
  vector<CModedOscillator>              	LM_all_Oscillators;
  vector<S_OscillatorExcitation>        	RM_Excitations;
  vector<CModedOscillator>              	RM_all_Oscillators;

  vector<CMasslessHalfState>            	MasslessRightMovers;
  vector<CHalfState>                    	RightMovers;
  vector<CHalfState>                    	TachyonicRightMovers;
  vector<CTachyonHalfState> 				TachyonMassRightMovers;

  CheckStatus                           	Sector_CheckStatus;
};

#endif
