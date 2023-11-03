#ifndef CPRINT_H
#define CPRINT_H

#include "groupTheory.hpp"
#include "cvector.h"
#include "coscillator.h"
#include "corbifoldgroupelement.h"
#include "corbifold.h"
#include "cstate.h"

class COrbifoldGroupElement;
class COrbifoldGroup;
class COrbifold;
class CSector;
class CFixedBrane;
class CState;
class CHalfState;
class CMasslessHalfState;
class CModedOscillator;
class CVector;
class CWilsonLines;
class CMonomial;
class CMassMatrix;
class CSpectrum;
class CInequivalentModels;

enum OutputType {Tstandard = 0, Tmathematica, Tlatex};

struct SSymmetryGroup;
struct SConfig;
struct GDT_Parameters;
struct SDiscreteSymmetry;
struct YukawaCoupling;
struct SSpectrum;


#ifndef TYPEDEF_CGAUGEGROUP
#define TYPEDEF_CGAUGEGROUP
typedef gaugeGroup<double> CGaugeGroup;
#endif

struct OutputStrings
{
  string cbegin;
  string cend;
  string endofset;
  string prelabel;
  string separator;
  string set_open;
  string set_close;
  string start_index;
  string underscore;
  string vector_open;
  string vector_close;
};


class CPrint{
public:
	// member functions
	CPrint(OutputType output_type, std::ostream *out);
	~CPrint();

	bool              SetOutputType(const OutputType &output_type);
	const OutputType &GetOutputType() const {return this->output_type;};

	bool     PrintSSpectrum(const SSpectrum &Spectrum, bool PrintU1 = false, int Position_of_and_in_GaugeGroup = -1) const;
	bool     PrintInequivalentModels(const vector<vector<CSpectrum> > &InequivalentModels) const;
	bool     PrintAllInequivalentModels(const CInequivalentModels &InequivalentModels) const;

	bool     PrintDiscreteSymmetry(const SDiscreteSymmetry &DiscreteSymmetry) const;
	bool     PrintModularSymmetry(const SModularSymmetry &ModularSymmetry) const;
	void     PrintDiscreteTorsion(const GDT_Parameters &DiscreteTorsion) const;

	bool     PrintCouplings(const SConfig &VEVConfig, const vector<YukawaCoupling> &FieldCouplings, bool EffectivePotential = false, bool PrintModuli = false) const;
	bool     PrintEigenvalues(const vector<COrbifoldGroupElement> &Centralizer, const vector<double> &Eigenvalues) const;
	bool     PrintFTerms(const SConfig &VEVConfig) const;
	bool     PrintGaugeGroupFactor(const SConfig &VEVConfig, const unsigned factor) const;
	void     PrintGaugeGroup(const SConfig &VEVConfig, bool print_selection = false) const;
	void     PrintHalfState(const CHalfState &HalfState, const CMasslessHalfState &MasslessHalfState, const vector<CModedOscillator> &all_Oscillators) const;
	void     PrintTachyonHalfState(const CHalfState &HalfState, const CTachyonHalfState &TachyonHalfState, const vector<CModedOscillator> &all_Oscillators) const;
	void     PrintIdentifiedWilsonLines(const vector<vector<unsigned> > &IdentifiedWilsonLines) const;

	bool     PrintLabel(const CField &Field, unsigned useLabel = 0) const;
	void     PrintListOfCharges(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const string &LabelOfList = "lstFields") const;
	bool     PrintLeftMovers(const CFixedBrane &FixedBrane, const vector<CModedOscillator> &all_Oscillators, const COrbifoldGroupElement &constructing_Element) const;

	void     PrintMasslessHalfState(const CMasslessHalfState &MasslessHalfState, const vector<CModedOscillator> &all_Oscillators, const SelfDualLattice &Lattice) const;
	bool     PrintMasslessLeftMovers(const CFixedBrane &FixedBrane, const vector<CModedOscillator> &all_Oscillators, const COrbifoldGroupElement &constructing_Element) const;
	bool     PrintMasslessRightMovers(const CSector &Sector) const;
	void     PrintMassMatrix(const CMassMatrix &MassMatrix, const SConfig &VEVConfig, unsigned max_order) const;
	void     PrintMassMatrixInfo(const CMassMatrix &MassMatrix, const SConfig &VEVConfig) const;
	unsigned PrintMonomial(const CMonomial &Monomial, const SConfig &VEVConfig, bool print_gauge_equivalent_fields) const;

	void     PrintOrbifoldGroup(const COrbifoldGroup &OrbifoldGroup) const;
	void     PrintOscillator(const CModedOscillator &Oscillator) const;
	void     PrintOscillatorExcitations(const CSector &Sector) const;

	void     PrintRational(const rational<int> &rat, bool indent = true) const;
	void     PrintRational(const doubleVector &Vector, SelfDualLattice Lattice = UNSPECIFIED_LATTICE, bool indent = true) const;
	void     PrintRational(const CVector &Vector, SelfDualLattice Lattice = UNSPECIFIED_LATTICE, bool indent = true) const;
	void     PrintRational(const CWilsonLines &WilsonLines, bool indent = true) const;
	bool     PrintRep(const RepVector &Dimensions) const;
	bool     PrintRep(const RepVector &Dimensions, const SSymmetryGroup &SymmetryGroup) const;
	bool     PrintRep(const CField &Field, const SSymmetryGroup &SymmetryGroup, bool U1_Charges = true) const;
	bool     PrintRightMovers(const CSector &Sector) const;

	bool     PrintSDimension(const SDimension &dim) const;
	bool     PrintSGElement(const CSpaceGroupElement &SGElement, bool formatted = true) const;
	void     PrintShift(const CShiftVector &ShiftVector) const;
	void     PrintSimpleRoots(const CGaugeGroup &GaugeGroup, int ith_simpleroot) const;
	bool     PrintSortedFieldLabels(const vector<unsigned> &FieldIndices, const SConfig &VEVConfig) const;
	void     PrintState(const CState &State, const CSector &Sector, const CFixedBrane &FixedBrane, const SConfig &VEVConfig) const;
	void     PrintTachyonicState(const CState &State, const CSector &Sector, const CFixedBrane &FixedBrane, const SConfig &VEVConfig, const bool &Rexcited) const;
	bool     PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const;
	bool     PrintTachyonicStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const;
	bool     PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, bool PrintInternalInformation = false) const;
	bool     PrintSpectrum(const vector<unsigned> &FieldIndices, const SConfig &VEVConfig, bool print_Labels = false) const;
	bool     PrintSpectrum(const CSpectrum &Spectrum) const;
//	bool     PrintSummary(const CFixedBrane &FixedBrane, const COrbifoldGroup &OrbifoldGroup, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet = LeftChiral, bool print_Labels = false) const;
    bool     PrintSummary(const CFixedBrane &FixedBrane, const COrbifoldGroup &OrbifoldGroup, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet, bool print_Labels = false) const;
//  bool     PrintSummary(const CSector &Sector, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet = LeftChiral, bool print_Labels = false) const;  
    bool     PrintSummary(const CSector &Sector, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet, bool print_Labels = false) const;  

//	void     PrintSummaryOfVEVConfig(const SConfig &VEVConfig, const SUSYMultiplet &Multiplet = LeftChiral, bool print_Labels = false) const;
	void     PrintSummaryOfVEVConfig(const SConfig &VEVConfig, const vector<SUSYMultiplet> &Multiplet, bool print_Labels = false) const;									//hacking here!!!
//	bool     PrintSummaryOfSectors(const COrbifold &Orbifold, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet = LeftChiral, bool print_Labels = false) const;
	bool     PrintSummaryOfSectors(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<SUSYMultiplet> &Multiplet, bool print_Labels = false) const;
//	void     PrintSummaryOfFixedBranes(const COrbifold &Orbifold, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet = LeftChiral, bool print_Labels = false) const;
    void     PrintSummaryOfFixedBranes(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<SUSYMultiplet> &Multiplet, bool print_Labels = false) const; //sept21

	bool     PrintSuperpotential(const SConfig &VEVConfig, bool EffectivePotential = false, bool PrintModuli = false) const;
	bool     PrintSUSYType(const SUSYMultiplet &Multiplet) const;

	void     PrintTwist(const CTwistVector &TwistVector) const;
	bool     PrintU1Charges(const CVector &U1Charges, const SSymmetryGroup &SymmetryGroup) const;
	void     PrintU1Directions(const SSymmetryGroup &SymmetryGroup, int ith_generator = -1) const;
	void     PrintVector(const CVector &Vector, SelfDualLattice Lattice = UNSPECIFIED_LATTICE) const;
	void     PrintWilsonLines(const CWilsonLines &WilsonLines, bool PrintInfo = true) const;

	bool     TexSpectrum(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const vector<unsigned> &print_Labels);

	// member variables
	std::ostream *out;

	string cbegin;
	string cend;
	string endofset;
	string prelabel;
	string separator;
	string set_open;
	string set_close;
	string start_index;
	string underscore;
	string vector_open;
	string vector_close;

private:
	OutputType    output_type;

	OutputStrings TstandardStrings;
	OutputStrings TmathematicaStrings;
	OutputStrings TlatexStrings;
};

#endif
