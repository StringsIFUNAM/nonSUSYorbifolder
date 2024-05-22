#include "cprint.h"
#include "corbifold.h"
#include "cspectrum.h"
#include "globalfunctions.h"
#include "cinequivalentspectra.h"
#include "cspectrum.h"


/* ########################################################################################
######   CPrint()                                                                    ######
######                                                                               ######
######   Version: 23.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) output_type : can be "Tstandard", "Tmathematica" or "Tlatex"             ######
######   2) out         : for example &cout                                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
CPrint::CPrint(OutputType output_type, std::ostream *out)
: out(out), output_type(output_type)
{
	this->TmathematicaStrings.cbegin = "(* ";
	this->TmathematicaStrings.cend   = " *)";
	this->TmathematicaStrings.endofset     = ";";
	this->TmathematicaStrings.prelabel     = "fld";
	this->TmathematicaStrings.separator    = ", ";
	this->TmathematicaStrings.set_open     = "{";
	this->TmathematicaStrings.set_close    = "}";
	this->TmathematicaStrings.start_index  = "";
	this->TmathematicaStrings.underscore   = "";
	this->TmathematicaStrings.vector_open  = "{";
	this->TmathematicaStrings.vector_close = "}";

	this->TstandardStrings.cbegin = "";
	this->TstandardStrings.cend   = "";
	this->TstandardStrings.endofset     = "";
	this->TstandardStrings.prelabel     = "";
	this->TstandardStrings.separator    = ", ";
	this->TstandardStrings.set_open     = "";
	this->TstandardStrings.set_close    = "";
	this->TstandardStrings.start_index  = "_";
	this->TstandardStrings.underscore   = "_";
	this->TstandardStrings.vector_open  = "(";
	this->TstandardStrings.vector_close = ")";

	this->TlatexStrings.cbegin = "";
	this->TlatexStrings.cend   = "";
	this->TlatexStrings.endofset     = "";
	this->TlatexStrings.prelabel     = "";
	this->TlatexStrings.separator    = ", ";
	this->TlatexStrings.set_open     = "";
	this->TlatexStrings.set_close    = "";
	this->TlatexStrings.start_index  = "_";
	this->TlatexStrings.underscore   = "_";
	this->TlatexStrings.vector_open  = "\\left(";
	this->TlatexStrings.vector_close = "\\right)";

	this->SetOutputType(output_type);
	(*this->out) << std::setiosflags(std::ios::fixed);
	(*this->out) << setprecision(2);
}



/* ########################################################################################
######   ~CPrint()                                                                   ######
######                                                                               ######
######   Version: 23.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
CPrint::~CPrint()
{
}



bool CPrint::PrintSSpectrum(const SSpectrum& Spectrum, bool PrintU1, int Position_of_and_in_GaugeGroup) const
{
	unsigned i = 0;
	unsigned j = 0;

	(*this->out) << "  ";
	const size_t g1 = Spectrum.GaugeGroup.size();
	for (i = 0; i < g1; ++i)
	{
		const unsigned G = Spectrum.GaugeGroup[i];
		const unsigned R = Spectrum.Rank[i];

		if (i != 0)
		{
			if (i == Position_of_and_in_GaugeGroup)
				(*this->out) << "and ";
			else
				(*this->out) << "+ ";
		}

		if (G == 1)
			(*this->out) << "SU(" << R+1 << ") ";
		else
			if (G == 2)
				(*this->out) << "SO(" << 2*R << ") ";
			else
				if (G == 3)
					(*this->out) << "E" << R << " ";
	}
	(*this->out) << "\n";

	size_t s1 = Spectrum.Spectrum.size();
	for (i = 0; i < s1; ++i)
	{
		const SRepresentation &CurentRep = Spectrum.Spectrum[i];

		(*this->out) << "  " << setw(3) << CurentRep.Multiplicity << " (";
		const vector<int> &Rep = CurentRep.Representation;
		if (g1 == 0)
			(*this->out) << "  1";
		else
		{
			for (j = 0; j < g1; ++j)
			{
				if (j != 0)
					(*this->out) << ",";
				(*this->out) << setw(3) << Rep[j];
			}
		}
		(*this->out) << ") ";

		if (PrintU1)
			this->PrintRational(D2Rat(CurentRep.U1Charge), false);

		(*this->out) << "\n";
	}
	//(*this->out) << "  Number of lines: " << s1 << "\n";
	(*this->out) << flush;
	return true;
}


bool CPrint::PrintInequivalentModels(const vector<vector<CSpectrum> > &InequivalentModels) const
{
	const size_t s1 = InequivalentModels.size();
	size_t s2 = 0;

	unsigned i = 0;
	unsigned j = 0;

	for (i = 0; i < s1; ++i)
	{
		const vector<CSpectrum> &SetOfInequivalentModels = InequivalentModels[i];

		s2 = SetOfInequivalentModels.size();
		if (s2 != 0)
		{
			(*this->out) << "\n  " << this->cbegin << "==============================================================================" << this->cend << "\n";
			(*this->out) << "  " << this->cbegin << "Non-empty set of inequivalent models:" << this->cend << "\n\n";
			for (j = 0; j < s2; ++j)
			{
				this->PrintSpectrum(SetOfInequivalentModels[j]);
				(*this->out) << endl;
			}
		}
	}
	return true;
}

bool CPrint::PrintAllInequivalentModels(const CInequivalentModels &InequivalentModels) const
{
	(*this->out) << "Case: mult%4 = 0, #ggf =  0" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank00);
	(*this->out) << "Case: mult%4 = 0, #ggf =  1" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank01);
	(*this->out) << "Case: mult%4 = 0, #ggf =  2" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank02);
	(*this->out) << "Case: mult%4 = 0, #ggf =  3" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank03);
	(*this->out) << "Case: mult%4 = 0, #ggf =  4" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank04);
	(*this->out) << "Case: mult%4 = 0, #ggf =  5" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank05);
	(*this->out) << "Case: mult%4 = 0, #ggf =  6" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank06);
	(*this->out) << "Case: mult%4 = 0, #ggf =  7" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank07);
	(*this->out) << "Case: mult%4 = 0, #ggf =  8" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank08);
	(*this->out) << "Case: mult%4 = 0, #ggf =  9" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank09);
	(*this->out) << "Case: mult%4 = 0, #ggf = 10" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank10);
	(*this->out) << "Case: mult%4 = 0, #ggf = 11" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank11);
	(*this->out) << "Case: mult%4 = 0, #ggf = 12" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank12);
	(*this->out) << "Case: mult%4 = 0, #ggf = 13" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank13);
	(*this->out) << "Case: mult%4 = 0, #ggf = 14" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank14);
	(*this->out) << "Case: mult%4 = 0, #ggf = 15" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank15);
	(*this->out) << "Case: mult%4 = 0, #ggf = 16" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c1_Rank16);

	(*this->out) << "Case: mult%4 = 1, #ggf =  0" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank00);
	(*this->out) << "Case: mult%4 = 1, #ggf =  1" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank01);
	(*this->out) << "Case: mult%4 = 1, #ggf =  2" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank02);
	(*this->out) << "Case: mult%4 = 1, #ggf =  3" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank03);
	(*this->out) << "Case: mult%4 = 1, #ggf =  4" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank04);
	(*this->out) << "Case: mult%4 = 1, #ggf =  5" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank05);
	(*this->out) << "Case: mult%4 = 1, #ggf =  6" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank06);
	(*this->out) << "Case: mult%4 = 1, #ggf =  7" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank07);
	(*this->out) << "Case: mult%4 = 1, #ggf =  8" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank08);
	(*this->out) << "Case: mult%4 = 1, #ggf =  9" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank09);
	(*this->out) << "Case: mult%4 = 1, #ggf = 10" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank10);
	(*this->out) << "Case: mult%4 = 1, #ggf = 11" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank11);
	(*this->out) << "Case: mult%4 = 1, #ggf = 12" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank12);
	(*this->out) << "Case: mult%4 = 1, #ggf = 13" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank13);
	(*this->out) << "Case: mult%4 = 1, #ggf = 14" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank14);
	(*this->out) << "Case: mult%4 = 1, #ggf = 15" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank15);
	(*this->out) << "Case: mult%4 = 1, #ggf = 16" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c2_Rank16);

	(*this->out) << "Case: mult%4 = 2, #ggf =  0" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank00);
	(*this->out) << "Case: mult%4 = 2, #ggf =  1" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank01);
	(*this->out) << "Case: mult%4 = 2, #ggf =  2" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank02);
	(*this->out) << "Case: mult%4 = 2, #ggf =  3" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank03);
	(*this->out) << "Case: mult%4 = 2, #ggf =  4" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank04);
	(*this->out) << "Case: mult%4 = 2, #ggf =  5" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank05);
	(*this->out) << "Case: mult%4 = 2, #ggf =  6" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank06);
	(*this->out) << "Case: mult%4 = 2, #ggf =  7" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank07);
	(*this->out) << "Case: mult%4 = 2, #ggf =  8" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank08);
	(*this->out) << "Case: mult%4 = 2, #ggf =  9" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank09);
	(*this->out) << "Case: mult%4 = 2, #ggf = 10" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank10);
	(*this->out) << "Case: mult%4 = 2, #ggf = 11" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank11);
	(*this->out) << "Case: mult%4 = 2, #ggf = 12" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank12);
	(*this->out) << "Case: mult%4 = 2, #ggf = 13" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank13);
	(*this->out) << "Case: mult%4 = 2, #ggf = 14" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank14);
	(*this->out) << "Case: mult%4 = 2, #ggf = 15" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank15);
	(*this->out) << "Case: mult%4 = 2, #ggf = 16" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c3_Rank16);

	(*this->out) << "Case: mult%4 = 3, #ggf =  0" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank00);
	(*this->out) << "Case: mult%4 = 3, #ggf =  1" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank01);
	(*this->out) << "Case: mult%4 = 3, #ggf =  2" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank02);
	(*this->out) << "Case: mult%4 = 3, #ggf =  3" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank03);
	(*this->out) << "Case: mult%4 = 3, #ggf =  4" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank04);
	(*this->out) << "Case: mult%4 = 3, #ggf =  5" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank05);
	(*this->out) << "Case: mult%4 = 3, #ggf =  6" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank06);
	(*this->out) << "Case: mult%4 = 3, #ggf =  7" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank07);
	(*this->out) << "Case: mult%4 = 3, #ggf =  8" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank08);
	(*this->out) << "Case: mult%4 = 3, #ggf =  9" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank09);
	(*this->out) << "Case: mult%4 = 3, #ggf = 10" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank10);
	(*this->out) << "Case: mult%4 = 3, #ggf = 11" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank11);
	(*this->out) << "Case: mult%4 = 3, #ggf = 12" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank12);
	(*this->out) << "Case: mult%4 = 3, #ggf = 13" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank13);
	(*this->out) << "Case: mult%4 = 3, #ggf = 14" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank14);
	(*this->out) << "Case: mult%4 = 3, #ggf = 15" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank15);
	(*this->out) << "Case: mult%4 = 3, #ggf = 16" << endl;
	this->PrintInequivalentModels(InequivalentModels.InequivalentModels_c4_Rank16);

	return true;
}



/* ########################################################################################
######   SetOutputType(const OutputType &output_type)                                ######
######                                                                               ######
######   Version: 17.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) output_type : can be "Tstandard", "Tmathematica" or "Tlatex"             ######
######   output:                                                                     ######
######   return value   : output type changed succesfully?                           ######
######################################################################################## */
bool CPrint::SetOutputType(const OutputType &output_type)
{
	this->output_type = output_type;

	if (this->output_type == Tmathematica)
	{
		this->cbegin       = this->TmathematicaStrings.cbegin;
		this->cend         = this->TmathematicaStrings.cend;
		this->endofset     = this->TmathematicaStrings.endofset;
		this->prelabel     = this->TmathematicaStrings.prelabel;
		this->separator    = this->TmathematicaStrings.separator;
		this->set_open     = this->TmathematicaStrings.set_open;
		this->set_close    = this->TmathematicaStrings.set_close;
		this->start_index  = this->TmathematicaStrings.start_index;
		this->underscore   = this->TmathematicaStrings.underscore;
		this->vector_open  = this->TmathematicaStrings.vector_open;
		this->vector_close = this->TmathematicaStrings.vector_close;
		return true;
	}
	if (this->output_type == Tstandard)
	{
		this->cbegin       = this->TstandardStrings.cbegin;
		this->cend         = this->TstandardStrings.cend;
		this->endofset     = this->TstandardStrings.endofset;
		this->prelabel     = this->TstandardStrings.prelabel;
		this->separator    = this->TstandardStrings.separator;
		this->set_open     = this->TstandardStrings.set_open;
		this->set_close    = this->TstandardStrings.set_close;
		this->start_index  = this->TstandardStrings.start_index;
		this->underscore   = this->TstandardStrings.underscore;
		this->vector_open  = this->TstandardStrings.vector_open;
		this->vector_close = this->TstandardStrings.vector_close;
		return true;
	}
	if (this->output_type == Tlatex)
	{
		this->cbegin       = this->TlatexStrings.cbegin;
		this->cend         = this->TlatexStrings.cend;
		this->endofset     = this->TlatexStrings.endofset;
		this->prelabel     = this->TlatexStrings.prelabel;
		this->separator    = this->TlatexStrings.separator;
		this->set_open     = this->TlatexStrings.set_open;
		this->set_close    = this->TlatexStrings.set_close;
		this->start_index  = this->TlatexStrings.start_index;
		this->underscore   = this->TlatexStrings.underscore;
		this->vector_open  = this->TlatexStrings.vector_open;
		this->vector_close = this->TlatexStrings.vector_close;
		return true;
	}

	return false;
}



/* ########################################################################################
######   PrintEigenvalues(const vector<COrbifoldGroupElement> &Centralizer,...) const######
######                                                                               ######
######   Version: 27.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Centralizer : contains the space group elements that commute with the    ######
######                    constructing element (i.e. the centralizer elements)       ######
######   2) Eigenvalues : the transformation phases of the state with respect to the ######
######                    centralizer elements                                       ######
######   output:                                                                     ######
######   return value   : printed succesfully?                                       ######
######################################################################################## */
bool CPrint::PrintEigenvalues(const vector<COrbifoldGroupElement> &Centralizer, const vector<double> &Eigenvalues) const
{
	(*this->out) << " Eigenvalues:\n\n centralizer \t\t\t eigenvalue\n";
	unsigned i = 0;
	const size_t s1 = Centralizer.size();
	const size_t s2 = Eigenvalues.size();

	if (s2 < s1)
	{
		cout << "\n  Warning in void CPrint::PrintEigenvalues(...) const: sizes of \"Centralizer\" and \"Eigenvalues\" are not equal. Return false." << endl;
		return false;
	}
	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &CElement = Centralizer[i].SGElement;

		(*this->out) << "(" << CElement.Get_m() << ", " << CElement.Get_n() << ") (" << ", " << CElement.Get_k() << ") ("
				<< CElement.Get_n(0) << ", " << CElement.Get_n(1) << ", "
				<< CElement.Get_n(2) << ", " << CElement.Get_n(3) << ", "
				<< CElement.Get_n(4) << ", " << CElement.Get_n(5) << ")";

		(*this->out) << "\t" << Eigenvalues[i] << "\n";
	}
	if (s1 != s2)
	{
		(*this->out) << " Gamma eigenvalues:\n";
		for (i = s1; i < s2; ++i)
			(*this->out) << Eigenvalues[i] << " ";
	}
	(*this->out) << endl;
	return true;
}


/* ########################################################################################
######   PrintGaugeGroupFactor(const SConfig &VEVConfig, const unsigned factor) const######
######                                                                               ######
######   Version: 27.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig : contains the gauge group to be printed                       ######
######   2) factor    : print the "factor"-th gauge group factor                     ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintGaugeGroupFactor(const SConfig &VEVConfig, const unsigned factor) const
{
	if (factor >= VEVConfig.SymmetryGroup.GaugeGroup.factor.size())
	{
		cout << "\n  Warning in bool CPrint::PrintGaugeGroupFactor(...) const : \"factor\" out of range. Return false." << endl;
		return false;
	}

	const gaugeGroupFactor<double> &ggf = VEVConfig.SymmetryGroup.GaugeGroup.factor[factor];

	if (ggf.algebra[0] == 'E')
		(*this->out) << "E_" << ggf.rank;
	else
	{
		if (ggf.algebra[0] == 'A')
			(*this->out) << "SU(" << ggf.rank+1 << ")";
		else
			if (ggf.algebra[0] == 'D')
				(*this->out) << "SO(" << 2 * ggf.rank << ")";
	}

	if (VEVConfig.SymmetryGroup.GGs_AdditionalLabels[factor] != "")
		(*this->out) << "_" << VEVConfig.SymmetryGroup.GGs_AdditionalLabels[factor];

	return true;
}



/* ########################################################################################
######   PrintGaugeGroup(const SConfig &VEVConfig, bool print_selection) const       ######
######                                                                               ######
######   Version: 27.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig       : contains the gauge group to be printed                 ######
######   2) print_selection : if true indicate by brackets [] which gauge group      ######
######                        factor is selected as part of the observable or hidden ######
######                        sector                                                 ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintGaugeGroup(const SConfig &VEVConfig, bool print_selection) const
{
	const CGaugeGroup &GaugeGroup = VEVConfig.SymmetryGroup.GaugeGroup;

	(*this->out) << "  " << this->cbegin << "Gauge group in vev-configuration \"" << VEVConfig.ConfigLabel << VEVConfig.ConfigNumber << "\": ";

	unsigned i = 0;

	if (print_selection)
	{
		const size_t t1 = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();
		const size_t t2 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();

		if ((t1 == 0) && (t2 == 0))
		{
			(*this->out) << "none" << this->cend << "\n" << endl;
			return;
		}

		bool print_brackets = false;
		bool first_factor_printed = false;

		for (i = 0; i < t1; ++i)
		{
			if (!first_factor_printed)
				first_factor_printed = true;
			else
			{
				if (i == VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup)
					(*this->out) << " and ";
				else
					(*this->out) << " x ";
			}

			print_brackets = (find(VEVConfig.SymmetryGroup.observable_sector_GGs.begin(), VEVConfig.SymmetryGroup.observable_sector_GGs.end(), i) == VEVConfig.SymmetryGroup.observable_sector_GGs.end());
			if (print_brackets)
				(*this->out) << "[";

			this->PrintGaugeGroupFactor(VEVConfig, i);

			if (print_brackets)
				(*this->out) << "]";
		}

		if ((t1 != 0) && (t2 != 0))
			(*this->out) << " and ";

		first_factor_printed = false;
		for (i = 0; i < t2; ++i )
		{
			if (!first_factor_printed)
				first_factor_printed = true;
			else
				(*this->out) << " x ";

			print_brackets = (find(VEVConfig.SymmetryGroup.observable_sector_U1s.begin(), VEVConfig.SymmetryGroup.observable_sector_U1s.end(), i) == VEVConfig.SymmetryGroup.observable_sector_U1s.end());

			if (print_brackets)
				(*this->out) << "[";

			(*this->out) << "U(1)_" << i+1;
			if (VEVConfig.SymmetryGroup.U1s_AdditionalLabels[i] != "")
				(*this->out) << "," << VEVConfig.SymmetryGroup.U1s_AdditionalLabels[i];

			if (print_brackets)
				(*this->out) << "]";
		}
		(*this->out) << this->cend;
	}
	else
	{
		const size_t g1 = VEVConfig.SymmetryGroup.observable_sector_GGs.size();
		const size_t u1 = VEVConfig.SymmetryGroup.observable_sector_U1s.size();

		if ((g1 == 0) && (u1 == 0))
		{
			(*this->out) << "none" << this->cend << "\n" << endl;
			return;
		}

		for (i = 0; i < g1; ++i)
		{
			this->PrintGaugeGroupFactor(VEVConfig, VEVConfig.SymmetryGroup.observable_sector_GGs[i]);

			if (i+1 < g1)
			{
				if (i+1 == VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup)
					(*this->out) << " and ";
				else
					(*this->out) << " x ";
			}
		}

		if (u1 != 0)
		{
			if (g1 != 0)
				(*this->out) << " and ";

			(*this->out) << "U(1)";
			if (u1 == 1)
			{
				if (VEVConfig.SymmetryGroup.U1s_AdditionalLabels[VEVConfig.SymmetryGroup.observable_sector_U1s[0]] != "")
					(*this->out) << "_" << VEVConfig.SymmetryGroup.U1s_AdditionalLabels[VEVConfig.SymmetryGroup.observable_sector_U1s[0]];
			}
			else
			{
				if (this->output_type == Tlatex)
					(*this->out) << "^{" << u1 << "}";
				else
					(*this->out) << "^" << u1;
			}
			(*this->out) << this->cend;

			if ((VEVConfig.SymmetryGroup.observable_sector_U1s[0] == 0) && VEVConfig.SymmetryGroup.IsFirstU1Anomalous) {		//hacking here!!!
				(*this->out) << "\n  " << this->cbegin << "First U(1) is anomalous with tr Q_anom = " << setw(5) << setprecision(2) << VEVConfig.SymmetryGroup.D0_FI_term << "." << this->cend;
			}
			/*else {
				(*this->out) << "\n  " << this->cbegin << "No anomalous U(1), tr Q_anom = 0" << this->cend;
			}*/
		}
	}
	(*this->out) << "\n" << endl;
}



/* ########################################################################################
######   PrintHalfState(const CHalfState &HalfState, ...) const                      ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) HalfState         : CHalfState object to be printed                      ######
######   2) MasslessHalfState : the corresponding CMasslessHalfState object that     ######
######                          contains the massless left-moving momenta p_sh       ######
######   3) all_Oscillators   : the set of all oscillators from the sector of this   ######
######                          "HalfState"                                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintHalfState(const CHalfState &HalfState, const CMasslessHalfState &MasslessHalfState, const vector<CModedOscillator> &all_Oscillators) const
{
	unsigned i = 0;
	size_t s1 = 0;

	(*this->out) << "  " << this->cbegin;
	if (HalfState.GetType() == LeftMover)
		(*this->out) << "left";
	else
		if (HalfState.GetType() == RightMover)
			(*this->out) << "right";

	const vector<CVector>        &all_Weights = MasslessHalfState.Weights;
	const S_OscillatorExcitation &Excitation  = MasslessHalfState.Excitation;
	(*this->out) << "-mover: N = " << Excitation.NumberOperator;

	s1 = Excitation.OscillatorIndices.size();
	if (s1 != 0)
	{
		(*this->out) << " with ";
		for (i = 0; i < s1; ++i)
		{
			this->PrintOscillator(all_Oscillators[Excitation.OscillatorIndices[i]]);
			(*this->out) << " ";
		}
	}

	s1 = HalfState.Weights.size();
	(*this->out) << ", #(weights) = " << s1 << this->cend << "\n";

	SelfDualLattice Lattice = E8xE8;

	for (i = 0; i < s1; ++i)
	{
		(*this->out) << "    ";
		this->PrintRational(all_Weights[HalfState.Weights[i]], Lattice);
		(*this->out) << "\n";
	}
}



/* ########################################################################################
######   PrintHalfState(const CHalfState &HalfState, ...) const                      ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) HalfState         : CHalfState object to be printed                      ######
######   2) MasslessHalfState : the corresponding CMasslessHalfState object that     ######
######                          contains the massless left-moving momenta p_sh       ######
######   3) all_Oscillators   : the set of all oscillators from the sector of this   ######
######                          "HalfState"                                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintTachyonHalfState(const CHalfState &HalfState, const CTachyonHalfState &TachyonHalfState, const vector<CModedOscillator> &all_Oscillators) const
{
	unsigned i = 0;
	size_t s1 = 0;

	(*this->out) << "  " << this->cbegin;
	if (HalfState.GetType() == LeftMover)
		(*this->out) << " left";
	else
		if (HalfState.GetType() == RightMover)
			(*this->out) << "  Tachyonic right";

	const vector<CVector>        &all_Weights = TachyonHalfState.Weights;
	const S_OscillatorExcitation &Excitation  = TachyonHalfState.Excitation;
	(*this->out) << "-mover: N = " << Excitation.NumberOperator;

	s1 = Excitation.OscillatorIndices.size();
	if (s1 != 0)
	{
		(*this->out) << " with ";
		for (i = 0; i < s1; ++i)
		{
			this->PrintOscillator(all_Oscillators[Excitation.OscillatorIndices[i]]);
			(*this->out) << " ";
		}
	}

	s1 = HalfState.Weights.size();
	(*this->out) << ", #(weights) = " << s1 << this->cend << "\n";

	SelfDualLattice Lattice = E8xE8;

	for (i = 0; i < s1; ++i)
	{
		(*this->out) << "    ";
		this->PrintRational(all_Weights[HalfState.Weights[i]], Lattice);
		(*this->out) << "\n";
	}
}



/* ########################################################################################
######   PrintIdentifiedWilsonLines(...) const                                       ######
######                                                                               ######
######   Version: 16.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) IdentifiedWilsonLines : sets of indices of Wilson lines that are         ######
######                              identified on the orbifold, e.g. W_1 = W_2 for Z3######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintIdentifiedWilsonLines(const vector<vector<unsigned> > &IdentifiedWilsonLines) const
{
	const size_t i1 = IdentifiedWilsonLines.size();

	(*this->out) << "  " << this->cbegin << "Wilson lines identified on the orbifold: ";
	if (i1 == 0)
	{
		(*this->out) << "none" << this->cend;
		return;
	}
	else
		(*this->out) << this->cend << "\n    ";

	if (i1 != 1)
		(*this->out) << this->set_open;

	for (unsigned i = 0; i < i1; ++i)
	{
		if (i != 0)
			(*this->out) << ", ";

		const vector<unsigned> &tmp_IdWilsonLines = IdentifiedWilsonLines[i];
		(*this->out) << this->set_open;
		(*this->out) << "W" << this->start_index << tmp_IdWilsonLines[0] + 1;

		const size_t i2 = tmp_IdWilsonLines.size();
		for (unsigned j = 1; j < i2; ++j)
			(*this->out) << " = W" << this->start_index << tmp_IdWilsonLines[j] + 1;
		(*this->out) << this->set_close;
	}

	if (i1 != 1)
		(*this->out) << this->set_close;

	(*this->out) << this->endofset;
}



/* ########################################################################################
######   PrintLeftMovers(const CFixedBrane &FixedBrane, ...) const                   ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane           : all left-movers of this fixed point/brane will be ######
######                             printed                                           ######
######   2) all_Oscillators      : the set of all oscillators from the sector of     ######
######                             this "FixedBrane"                                 ######
######   3) constructing_Element : contains the local shift                          ######
######   output:                                                                     ######
######   return value            : printed succesfully?                              ######
######################################################################################## */
bool CPrint::PrintLeftMovers(const CFixedBrane &FixedBrane, const vector<CModedOscillator> &all_Oscillators, const COrbifoldGroupElement &constructing_Element) const
{
	unsigned i = 0;

	(*this->out) << "  fixed brane (";
	(*this->out) << constructing_Element.SGElement.Get_m() << ", " << constructing_Element.SGElement.Get_n() << ", " << constructing_Element.SGElement.Get_k() << ") (";
	for (i = 0; i < 6; ++i)
	{
		(*this->out) << constructing_Element.SGElement.Get_n(i);
		if (i != 5)
			(*this->out) << ", ";
	}
	(*this->out) << ")\n  local ";

	(*this->out) << "Shift = ";
	this->PrintShift(constructing_Element.Shift);
	(*this->out) << "\n";

	const size_t s1 = FixedBrane.GetNumberOfLeftMovers();
	for (i = 0; i < s1; ++i)
	{
		const CHalfState &LeftMover = FixedBrane.GetLeftMover(i);
		this->PrintHalfState(LeftMover, FixedBrane.GetMasslessLeftMover(LeftMover.GetIndex()), all_Oscillators);
		(*this->out) << endl;
	}

	return true;
}



/* ########################################################################################
######   PrintRightMovers(const CSector &Sector) const                               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector    : all right-movers of this (un-)twisted sector will be printed ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintRightMovers(const CSector &Sector) const
{
	unsigned i = 0;

	(*this->out) << "  sector (" << Sector.Get_m() << ", " << Sector.Get_n() << ", " << Sector.Get_k() << ")\n  local Twist = ";
	this->PrintTwist(Sector.GetTwist());
	(*this->out) << "\n";

	const size_t s1 = Sector.GetNumberOfRightMovers();
	for (i = 0; i < s1; ++i)
	{
		const CHalfState &RightMover = Sector.GetRightMover(i);
		this->PrintHalfState(RightMover, Sector.GetMasslessRightMover(RightMover.GetIndex()), Sector.GetRM_all_Oscillators());
		(*this->out) << endl;
	}
	return true;
}



/* ########################################################################################
######   PrintLabel(const CField &Field, unsigned useLabel) const                    ######
######                                                                               ######
######   Version: 18.11.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Field     : field whose label shall be printed                           ######
######   2) useLabel  : print the "useLabel"-th label of "Field"                     ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintLabel(const CField &Field, unsigned useLabel) const
{
	if (useLabel >= Field.Labels.size())
	{
		cout << "\n  Warning in bool CPrint::PrintLabel(...) const: UseLabel-index out of range. Return false." << endl;
		return false;
	}

	if (Field.Labels[useLabel] == "")
		return false;

	(*this->out) << this->prelabel << Field.Labels[useLabel] << this->start_index;

	if (output_type == Tlatex)
		(*this->out) << "{" << Field.Numbers[useLabel] << "}";
	else
		(*this->out) << Field.Numbers[useLabel];

	return true;
}



/* ########################################################################################
######   PrintListOfCharges(const COrbifold &Orbifold, ...) const                    ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold     : the orbifold of the vev-config "VEVConfig"                ######
######   2) VEVConfig    : contains the massless fields                              ######
######   3) FieldIndices : indices of fields to be printed                           ######
######   4) LabelOfList  : assign a label to the list                                ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintListOfCharges(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const string &LabelOfList) const
{
	bool for_mathematica = (this->output_type == Tmathematica);

	const vector<SDiscreteSymmetry> &DiscreteNonRSymmetries = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteNonRSymmetries;
	const vector<SDiscreteSymmetry> &DiscreteRSymmetries    = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteRSymmetries;

	const size_t ds1 = DiscreteNonRSymmetries.size();
	const size_t ds2 = DiscreteRSymmetries.size();

	const bool print_labels = true;
	const bool print_Rcharges     = (ds2 != 0);
	const bool print_localization = (ds1 != 0);

	const size_t f1 = FieldIndices.size();

	const SelfDualLattice Lattice = Orbifold.OrbifoldGroup.GetShift(0).GetLattice();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	bool ListLabelPrinted = false;
	if (LabelOfList != "")
	{
		(*this->out) << "\n  " << LabelOfList << " = ";
		ListLabelPrinted = true;
	}

	if (for_mathematica)
	{
		if (!ListLabelPrinted)
			(*this->out) << "\n  lstFields = ";
		(*this->out) << "{";
	}

	string label1 = "";
	string label2 = "";

	bool firstLine = true;

	for (i = 0; i < f1; ++i)
	{
		const CField &tmp_Field = VEVConfig.Fields[FieldIndices[i]];

		const size_t s1 = tmp_Field.GetNumberOfLMWeights();
		for (j = 0; j < s1; ++j)
		{
			const CVector &tmp_Vector = tmp_Field.GetLMWeight(j, Orbifold.GetSectors());

			if (for_mathematica)
			{
				if (firstLine)
				{
					firstLine = false;
					(*this->out) << "\n  {";
				}
				else
					(*this->out) << ",\n  {";
			}
			else
				(*this->out) << "\n  ";

			this->PrintRational(tmp_Vector, Lattice);
			(*this->out) << this->separator << " ";

			this->PrintRational(tmp_Field.GetRMWeight(0,Orbifold.GetSectors()), SO8); 
			(*this->out) << this->separator << " ";

			if (print_localization)
			{
				(*this->out) << this->vector_open;
				for (k = 0; k < ds1; ++k)
				{
					this->PrintRational(tmp_Field.GetDiscreteCharge(DiscreteNonRSymmetries[k]));
					if (k < ds1-1)
						(*this->out) << this->separator;
				}
				(*this->out) << this->vector_close << this->separator << " ";
			}

			if (print_labels)
				(*this->out) << "\"" << this->prelabel << tmp_Field.Labels[VEVConfig.use_Labels] << this->start_index << tmp_Field.Numbers[VEVConfig.use_Labels] << "\"";
			if (for_mathematica)
				(*this->out) << "}";
		}
	}
	if (for_mathematica)
		(*this->out) << "};";
	(*this->out) << endl;
}



/* ########################################################################################
######   PrintMasslessHalfState(...) const                                           ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) MasslessHalfState : CMasslessHalfState object to be printed              ######
######   2) all_Oscillators   : the set of all oscillators from the corresponding    ######
######                          (un-)twisted sector                                  ######
######   3) Lattice           : E8xE8 or Spin32                                      ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintMasslessHalfState(const CMasslessHalfState &MasslessHalfState, const vector<CModedOscillator> &all_Oscillators, const SelfDualLattice &Lattice) const
{
	unsigned i = 0;
	size_t s1 = 0;

	(*this->out) << "  " << this->cbegin << "massless ";
	if (MasslessHalfState.Type == LeftMover)
		(*this->out) << "left";
	else
		if (MasslessHalfState.Type == RightMover)
			(*this->out) << "right";

	(*this->out) << "-mover:" << this->cend << " N = ";
	this->PrintRational(D2Rat(MasslessHalfState.Excitation.NumberOperator), false);

	s1 = MasslessHalfState.Excitation.OscillatorIndices.size();
	if (s1 != 0)
	{
		(*this->out) << "\n  " << this->cbegin << "with oscillator(s): ";
		for (i = 0; i < s1; ++i)
		{
			(*this->out) << " ";
			this->PrintOscillator(all_Oscillators[MasslessHalfState.Excitation.OscillatorIndices[i]]);
		}
	}

	s1 = MasslessHalfState.Weights.size();
	(*this->out) << ", #(weights) = " << s1 << this->cend << "\n";

	string endofline = "";
	if (this->output_type == Tmathematica)
		endofline = ",";
	else
		if (this->output_type == Tlatex)
			endofline = "\\\\";

	string eqnarray = "";
	if (this->output_type == Tlatex)
	{
		if (s1 == 1)
			(*this->out) << "\\begin{equation}\n";
		else
		{
			eqnarray = "& &";
			(*this->out) << "\\begin{eqnarray}\n";
		}
	}

	for (i = 0; i < s1; ++i)
	{
		(*this->out) << "    " << eqnarray;
		if ((i == 0) && (s1 != 1))
			(*this->out) << this->set_open;

		this->PrintRational(MasslessHalfState.Weights[i], Lattice);

		if (i+1 < s1)
			(*this->out) << endofline << "\n";
	}
	if (s1 != 1)
		(*this->out) << this->set_close;

	(*this->out) << this->endofset << "\n";

	if (this->output_type == Tlatex)
	{
		if (s1 == 1)
			(*this->out) << "\\end{equation}\n";
		else
			(*this->out) << "\\end{eqnarray}\n";
	}
}



/* ########################################################################################
######   PrintMasslessLeftMovers(...) const                                          ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane           : all CMasslessHalfState objects of this fixed      ######
######                             point /brane will be printed                      ######  
######   2) all_Oscillators      : the set of all oscillators from the corresponding ######
######                             (un-)twisted sector                               ######
######   3) constructing_Element : contains the local shift                          ######
######   output:                                                                     ######
######   return value            : printed succesfully?                              ######
######################################################################################## */
bool CPrint::PrintMasslessLeftMovers(const CFixedBrane &FixedBrane, const vector<CModedOscillator> &all_Oscillators, const COrbifoldGroupElement &constructing_Element) const
{
	unsigned i = 0;
	const SelfDualLattice &Lattice = constructing_Element.Shift.Lattice;

	(*this->out) << "  " << this->cbegin << "fixed brane (";
	(*this->out) << constructing_Element.SGElement.Get_m() << ", " << constructing_Element.SGElement.Get_n() << ") (" << ", " << constructing_Element.SGElement.Get_k() << ") (";
	for (i = 0; i < 6; ++i)
	{
		this->PrintRational(constructing_Element.SGElement.Get_n(i), false);
		if (i != 5)
			(*this->out) << ", ";
	}
	(*this->out) << ")" << this->cend << "\n  " << this->cbegin << "local Shift:" << this->cend << "\n";
	if (this->output_type == Tlatex)
		(*this->out) << "\\begin{equation}\n";
	(*this->out) << "    V = ";
	this->PrintRational(constructing_Element.Shift, Lattice);
	(*this->out) << this->endofset << "\n";
	if (this->output_type == Tlatex)
		(*this->out) << "\\end{equation}\n";

	const size_t s1 = FixedBrane.GetNumberOfMasslessLeftMovers();
	for (i = 0; i < s1; ++i)
	{
		this->PrintMasslessHalfState(FixedBrane.GetMasslessLeftMover(i), all_Oscillators, Lattice);
		(*this->out) << "\n";
	}
	(*this->out) << flush;
	return true;
}



/* ########################################################################################
######   PrintMasslessRightMovers(const CSector &Sector) const                       ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector    : all CMasslessHalfState objects of this (un-)twisted sector   ######
######                  will be printed                                              ######  
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintMasslessRightMovers(const CSector &Sector) const
{
	unsigned i = 0;

	(*this->out) << "  " << this->cbegin << "sector (" << Sector.Get_m() << ", " << Sector.Get_n() << ", " << Sector.Get_k() << ")" << this->cend << "\n";
	(*this->out) << "  " << this->cbegin << "local Twist:" << this->cend << "\n";
	if (this->output_type == Tlatex)
		(*this->out) << "\\begin{equation}\n";
	(*this->out) << "    v = ";
	this->PrintRational(Sector.GetTwist(), SO8);
	(*this->out) << this->endofset << "\n";
	if (this->output_type == Tlatex)
		(*this->out) << "\\end{equation}\n";

	const size_t s1 = Sector.GetNumberOfMasslessRightMovers();
	for (i = 0; i < s1; ++i)
	{
		this->PrintMasslessHalfState(Sector.GetMasslessRightMover(i), Sector.GetRM_all_Oscillators(), SO8);
		(*this->out) << "\n";
	}
	(*this->out) << flush;
	return true;
}


/* ########################################################################################
######   PrintOrbifoldGroup(const COrbifoldGroup &OrbifoldGroup) const               ######
######                                                                               ######
######   Version: 25.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup : COrbifoldGroup object to be printed                      ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintOrbifoldGroup(const COrbifoldGroup &OrbifoldGroup) const
{
	OrbifoldGroup.GetSpaceGroup().PrintPointGroup((*this->out));
	(*this->out) << "\n";

	const SelfDualLattice Lattice = OrbifoldGroup.GetShift(0).GetLattice();

	if (Lattice == E8xE8)
		(*this->out) << "of E8 x E8";

	if (Lattice == Spin32)
		(*this->out) << "of spin32/Z2";

	(*this->out) << " heterotic string\n\n";

	if (OrbifoldGroup.GetNumberOfSupersymmetry() != -1)
		(*this->out) << "unbroken susy: N = " << OrbifoldGroup.GetNumberOfSupersymmetry() << "\n";

	if (OrbifoldGroup.GetSpaceGroup().IsZMxZN())
	{
		(*this->out) << "Shift0:\n";
		this->PrintShift(OrbifoldGroup.GetShift(0));
		(*this->out) << "\nShift1:\n";
		this->PrintShift(OrbifoldGroup.GetShift(1));

		(*this->out) << "\nTwist0:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(0));
		(*this->out) << "\nTwist1:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(1));
	}
	else if (OrbifoldGroup.GetSpaceGroup().IsZMxZNxZK())
	{
		(*this->out) << "Shift0:\n";
		this->PrintShift(OrbifoldGroup.GetShift(0));
		(*this->out) << "\nShift1:\n";
		this->PrintShift(OrbifoldGroup.GetShift(1));
		(*this->out) << "\nShift2:\n";
		this->PrintShift(OrbifoldGroup.GetShift(2));

		(*this->out) << "\nTwist0:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(0));
		(*this->out) << "\nTwist1:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(1));
		(*this->out) << "\nTwist2:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(2));
	}
	else
	{
		(*this->out) << "Shift:\n";
		this->PrintShift(OrbifoldGroup.GetShift(0));

		(*this->out) << "\nTwist:\n";
		this->PrintTwist(OrbifoldGroup.GetSpaceGroup().GetTwist(0));
	}
	(*this->out) << "\n\n";
	this->PrintWilsonLines(OrbifoldGroup.GetWilsonLines());
	(*this->out) << endl;

	if (OrbifoldGroup.UseFreelyActingWL)
	{
		(*this->out) << "using freely acting WL:\n";
		this->PrintRational(OrbifoldGroup.FreelyActingWilsonLine, Lattice);
		(*this->out) << endl;
	}
}



/* ########################################################################################
######   PrintOscillator(const CModedOscillator &Oscillator) const                   ######
######                                                                               ######
######   Version: 17.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Oscillator : CModedOscillator object to be printed                       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintOscillator(const CModedOscillator &Oscillator) const
{
	if (this->output_type == Tlatex)
	{
		if (Oscillator.GetType() == LeftMover)
			(*this->out) << "\\tilde";

		(*this->out) << "\\alpha_{";
		this->PrintRational(D2Rat(Oscillator.GetFrequency()), false);
		(*this->out) << "}";

		if (Oscillator.GetComplex())
			(*this->out) << "^{\\overline{" << Oscillator.GetIndex() << "}}";
		else
			(*this->out) << "^{" << Oscillator.GetIndex() << "}";
		return;
	}

	//if (this->output_type == Tstandard)
	{
		if (Oscillator.GetType() == LeftMover)
			(*this->out) << "-";

		(*this->out) << "a_(";
		this->PrintRational(D2Rat(Oscillator.GetFrequency()));
		(*this->out) << ")";

		if (Oscillator.GetComplex())
			(*this->out) << "^(-" << Oscillator.GetIndex() << ")";
		else
			(*this->out) << "^" << Oscillator.GetIndex();
		return;
	}
}



/* ########################################################################################
######   PrintOscillatorExcitations(const CSector &Sector) const                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector : all oscillator excitations of this sector will be printed       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintOscillatorExcitations(const CSector &Sector) const
{
	(*this->out) << "(" << Sector.Get_m() << ", " << Sector.Get_n() << ", " << Sector.Get_k() << ") with local twist: ";
	this->PrintTwist(Sector.GetTwist());
	(*this->out) << "\nexcitations of the left-mover:\n";

	unsigned i = 0;
	unsigned j = 0;

	const vector<S_OscillatorExcitation> &LM_Excitations = Sector.GetLM_Excitations();
	size_t s1 = LM_Excitations.size();
	size_t s2 = 0;

	if (s1 == 0)
		cout << "\n  Warning in void CPrint::PrintOscillatorExcitations(...) const: even the no oscillator case is not known." << endl;

	for (i = 0; i < s1; ++i)
	{
		const S_OscillatorExcitation &LM_Excitation = LM_Excitations[i];
		(*this->out) << "N = " << LM_Excitation.NumberOperator;

		s2 = LM_Excitation.OscillatorIndices.size();
		if (s2 != 0)
		{
			(*this->out) << " with ";
			for (j = 0; j < s2; ++j)
			{
				this->PrintOscillator(Sector.GetLM_Oscillator(LM_Excitation.OscillatorIndices[j]));
				(*this->out) << " ";
			}
		}
		(*this->out) << "\n";
	}

	(*this->out) << "\nexcitations of the right-mover:\n";

	const vector<S_OscillatorExcitation> &RM_Excitations = Sector.GetRM_Excitations();
	s1 = RM_Excitations.size();
	if (s1 == 0)
		cout << "\n  Warning in void CPrint::PrintOscillatorExcitations(...) const: even the no oscillator case is not known." << endl;

	for (i = 0; i < s1; ++i)
	{
		const S_OscillatorExcitation &RM_Excitation = RM_Excitations[i];
		(*this->out) << "N = " << RM_Excitation.NumberOperator;

		s2 = RM_Excitation.OscillatorIndices.size();
		if (s2 != 0)
		{
			(*this->out) << " with ";
			for (j = 0; j < s2; ++j)
			{
				this->PrintOscillator(Sector.GetRM_Oscillator(RM_Excitation.OscillatorIndices[j]));
				(*this->out) << " ";
			}
		}
		(*this->out) << "\n";
	}

	(*this->out) << endl;
}



/* ########################################################################################
######   PrintRational(const rational<int> &rat, bool indent) const                  ######
######                                                                               ######
######   Version: 17.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) rat    : a rational number to be printed                                 ######
######   2) indent : if true indent the output                                       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintRational(const rational<int> &rat, bool indent) const
{
	if (indent)
	{
		if (rat.denominator() == 1)
			(*this->out) << setw(5) << rat.numerator();
		else
		{
			if (this->output_type == Tlatex)
			{
				if (rat < 0)
					(*this->out) << "-\\tfrac{" << -rat.numerator() << "}{" << rat.denominator() << "}";
				else
					(*this->out) << "\\tfrac{" << rat.numerator() << "}{" << rat.denominator() << "}";
			}
			else
				(*this->out) << setw(3) << rat.numerator() << "/" << rat.denominator();
		}
	}
	else
	{
		if (rat.denominator() == 1)
			(*this->out) << rat.numerator();
		else
		{
			if (this->output_type == Tlatex)
			{
				if (rat < 0)
					(*this->out) << "-\\tfrac{" << -rat.numerator() << "}{" << rat.denominator() << "}";
				else
					(*this->out) << "\\tfrac{" << rat.numerator() << "}{" << rat.denominator() << "}";
			}
			else
				(*this->out) << rat.numerator() << "/" << rat.denominator();
		}
	}
}



/* ########################################################################################
######   PrintRational(const CVector &Vector, SelfDualLattice Lattice, ...) const    ######
######                                                                               ######
######   Version: 24.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector  : CVector object to be printed as vector of rational numbers     ######
######   2) Lattice : E8xE8 or Spin32                                                ######
######   3) indent  : if true indent the output                                      ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintRational(const CVector &Vector, SelfDualLattice Lattice, bool indent) const
{
	if (Vector.GetSize() == 0)
		return;

	if (Lattice != UNSPECIFIED_LATTICE)
		(*this->out) << this->vector_open;

	rational<int> tmp;
	unsigned i = 0;
	if ((Vector.GetSize() == 16) && (this->output_type != Tmathematica))
	{
		for (i = 0; i < 8; ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < 7)
				(*this->out) << this->separator;
		}

		if (Lattice == E8xE8)
			(*this->out) << this->vector_close << this->separator << " " << this->vector_open;
		else
			(*this->out) << this->separator;

		for (i = 8; i < 16; ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < 15)
				(*this->out) << this->separator;
		}
	}
	else
	{
		for (i = 0; i < Vector.GetSize(); ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < Vector.GetSize()-1)
				(*this->out) << this->separator;
		}
	}

	if (Lattice != UNSPECIFIED_LATTICE)
		(*this->out) << this->vector_close;
}



/* ########################################################################################
######   PrintRational(const doubleVector &Vector, ...) const                        ######
######                                                                               ######
######   Version: 12.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector  : CVector object to be printed as vector of rational numbers     ######
######   2) Lattice : E8xE8 or Spin32                                                ######
######   3) indent  : if true indent the output                                      ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintRational(const doubleVector &Vector, SelfDualLattice Lattice, bool indent) const
{
	if (Vector.size() == 0)
		return;

	(*this->out) << this->vector_open;

	rational<int> tmp;
	unsigned i = 0;
	if ((Vector.size() == 16) && (this->output_type != Tmathematica))
	{
		for (i = 0; i < 8; ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < 7)
				(*this->out) << this->separator;
		}

		if (Lattice == E8xE8)
			(*this->out) << this->vector_close << " " << this->vector_open;

		for (i = 8; i < 16; ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < 15)
				(*this->out) << this->separator;
		}
	}
	else
	{
		for (i = 0; i < Vector.size(); ++i)
		{
			this->PrintRational(D2Rat(Vector[i]), indent);
			if (i < Vector.size()-1)
				(*this->out) << this->separator;
		}
	}

	(*this->out) << this->vector_close;
}



/* ########################################################################################
######   PrintRational(const CWilsonLines &WilsonLines, bool indent) const           ######
######                                                                               ######
######   Version: 12.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) WilsonLines : to be printed as vectors of rational numbers               ######
######   2) indent      : if true indent the output                                  ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintRational(const CWilsonLines &WilsonLines, bool indent) const
{
	const bool for_mathematica = (this->output_type == Tmathematica);
	const bool for_latex       = (this->output_type == Tlatex);

	if (for_latex)
		(*this->out) << "\\begin{eqnarray}\n";

	for(unsigned i = 0; i < 6; ++i)
	{
		if (for_latex)
			(*this->out) << "  W_{" << i+1 << "} & = & ";
		else
			if (for_mathematica)
				(*this->out) << "  W" << i+1 << " = ";
			else
				(*this->out) << "  W_" << i+1 << " = ";

		this->PrintRational(WilsonLines.GetWilsonLine(i), WilsonLines.GetWilsonLine(i).GetLattice(), indent);

		if (for_latex)
			(*this->out) << "\\\\";
		else
			if (for_mathematica)
				(*this->out) << ";";

		(*this->out) << "\n";
	}
	if (for_latex)
		(*this->out) << "\\end{eqnarray}\n";
}



/* ########################################################################################
######   PrintRep(const CField &Field, ...) const                                    ######
######                                                                               ######
######   Version: 12.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Field         : representation of this CField object to be printed       ######
######   2) SymmetryGroup : contains the choice of observable sector                 ######
######   3) U1_Charges    : print the U(1) charges of the observable sector if true  ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
bool CPrint::PrintRep(const CField &Field, const SSymmetryGroup &SymmetryGroup, bool U1_Charges) const
{
	const bool for_mathematica = (this->output_type == Tmathematica);

	if (for_mathematica)
		(*this->out) << "{";

	this->PrintRep(Field.Dimensions, SymmetryGroup);
	this->PrintSUSYType(Field.Multiplet);

	if (U1_Charges && (SymmetryGroup.observable_sector_U1s.size() != 0))
	{
		if (for_mathematica)
			(*this->out) << ", ";
		else
			(*this->out) << "  U(1): ";
		this->PrintU1Charges(Field.U1Charges, SymmetryGroup);
	}

	if (for_mathematica)
		(*this->out) << "}";

	return true;
}



/* ########################################################################################
######   PrintRep(const RepVector &Dimensions, ...) const                            ######
######                                                                               ######
######   Version: 06.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Dimensions    : dimensions of non-Abelian representations                ######
######   2) SymmetryGroup : contains the choice of observable sector                 ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
bool CPrint::PrintRep(const RepVector &Dimensions, const SSymmetryGroup &SymmetryGroup) const
{
	const size_t s1 = SymmetryGroup.observable_sector_GGs.size();
	if (s1 == 0)
	{
		SDimension dim;
		dim.Dimension = 1;
		dim.AdditionalLabel = "";
		this->PrintSDimension(dim);
		return true;
	}

	// representation is charged with respect to non-abelian gauge groups
	(*this->out) << this->vector_open;
	for (unsigned i = 0; i < s1; ++i)
	{
		if (i != 0) (*this->out) << ",";
		this->PrintSDimension(Dimensions[SymmetryGroup.observable_sector_GGs[i]]);
	}
	(*this->out) << this->vector_close;

	return true;
}



/* ########################################################################################
######   PrintRep(const RepVector &Dimensions) const                                 ######
######                                                                               ######
######   Version: 01.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RepVector : dimensions of non-Abelian representations                    ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintRep(const RepVector &Dimensions) const
{
	const size_t s1 = Dimensions.size();
	if (s1 == 0)
	{
		SDimension dim;
		dim.Dimension = 1;
		dim.AdditionalLabel = "";
		this->PrintSDimension(dim);
		return true;
	}

	// representation is charged with respect to non-abelian gauge groups
	(*this->out) << this->vector_open;
	for (unsigned i = 0; i < s1; ++i)
	{
		if (i != 0) (*this->out) << ",";
		this->PrintSDimension(Dimensions[i]);
	}
	(*this->out) << this->vector_close;

	return true;
}



/* ########################################################################################
######   PrintSDimension(const SDimension &dim) const                                ######
######                                                                               ######
######   Version: 16.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) dim       : dimension of non-Abelian representation                      ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintSDimension(const SDimension &dim) const
{
	if (this->output_type == Tstandard)
	{
		(*this->out) << setw(3) << dim.Dimension;

		if (dim.AdditionalLabel != "")
			(*this->out) << "_" << dim.AdditionalLabel;

		return true;
	}

	if (this->output_type == Tlatex)
	{
		if (dim.Dimension < 0)
			(*this->out) << setw(3) << "\\crep{" << -dim.Dimension << "}";
		else
			(*this->out) << setw(3) << "\\rep{" << dim.Dimension << "}";

		if (dim.AdditionalLabel != "")
			(*this->out) << "_\\text{" << dim.AdditionalLabel << "}";

		return true;
	}

	if (this->output_type == Tmathematica)
	{
		(*this->out) << setw(3) << dim.Dimension;
		(*this->out) << dim.AdditionalLabel;

		return true;
	}

	return false;
}



/* ########################################################################################
######   PrintSGElement(const CSpaceGroupElement &SGElement, bool formatted) const   ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement : CSpaceGroupElement object to be printed                      ######
######   2) formatted : used in Tstandard output style                               ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintSGElement(const CSpaceGroupElement &SGElement, bool formatted) const
{
	if (this->output_type == Tlatex)
	{
		const bool mnotnull = (SGElement.Get_m() != 0);
		const bool nnotnull = (SGElement.Get_n() != 0);
		const bool knotnull = (SGElement.Get_k() != 0);

		(*this->out) << "\\left(";
		if (mnotnull || nnotnull || knotnull)
		{
			if (mnotnull)
			{
				(*this->out) << "\\theta";
				if (SGElement.Get_m() != 1)
					(*this->out) << "^{" << SGElement.Get_m() << "} ";
			}
			if (nnotnull)
			{
				(*this->out) << "\\omega";
				if (SGElement.Get_n() != 1)
					(*this->out) << "^{" << SGElement.Get_n() << "} ";
			}
			if (knotnull)
			{
				(*this->out) << "\\phi";
				if (SGElement.Get_k() != 1)
					(*this->out) << "^{" << SGElement.Get_k() << "} ";
			}
		}
		else
			(*this->out) << "1";

		(*this->out) << ",";
		bool first_printed = false;
		for (unsigned i = 0; i < LatticeDim; ++i)
		{
			if (SGElement.Get_n(i) != 0)
			{
				if (first_printed)
				{
					if (SGElement.Get_n(i) > 0)
						(*this->out) << " +";
				}
				first_printed = true;

				if (SGElement.Get_n(i) < 0)
					(*this->out) << " -";

				if (SGElement.Get_n(i) != 1)
				{
					(*this->out) << " ";
					this->PrintRational(abs(SGElement.Get_n(i)));
				}
				(*this->out) << " e_" << i+1;
			}
		}
		if (!first_printed)
			(*this->out) << " 0";
		(*this->out) << "\\right)";
		return true;
	}

	if (this->output_type == Tmathematica)
	{
		(*this->out) << "{" << SGElement.Get_m() << ", " << SGElement.Get_n() << ", " << SGElement.Get_k() << ", ";
		for (unsigned i = 0; i < LatticeDim; ++i)
		{
			this->PrintRational(SGElement.Get_n(i));
			if (i != LatticeDim-1)
				(*this->out) << ", ";
		}
		(*this->out) << "}";
		return true;
	}

	if (formatted)
	{
		(*this->out) << "(" << SGElement.Get_m() << ", " << SGElement.Get_n() << ", " << SGElement.Get_k() << ") (";
		for (unsigned i = 0; i < LatticeDim; ++i)
		{
			this->PrintRational(SGElement.Get_n(i));
			if (i != LatticeDim-1)
				(*this->out) << ", ";
		}
		(*this->out) << ")";
		return true;
	}

	(*this->out) << SGElement.Get_m() << " " << SGElement.Get_n() << " " << SGElement.Get_k() << " ";
	for (unsigned i = 0; i < LatticeDim; ++i)
	{
		this->PrintRational(SGElement.Get_n(i));
		if (i != LatticeDim-1)
			(*this->out) << " ";
	}
	return true;
}



/* ########################################################################################
######   PrintShift(const CShiftVector &ShiftVector) const                           ######
######                                                                               ######
######   Version: 22.09.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ShiftVector : CShiftVector object to be printed                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintShift(const CShiftVector &ShiftVector) const
{
	this->PrintRational(ShiftVector, ShiftVector.Lattice);
}



/* ########################################################################################
######   PrintSimpleRoots(const CGaugeGroup &GaugeGroup, int ith_simpleroot) const   ######
######                                                                               ######
######   Version: 31.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) GaugeGroup     : contains the simple roots                               ######
######   2) ith_simpleroot : if negative print all simple roots of "GaugeGroup",     ######
######                       else only the "ith_simpleroot"                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintSimpleRoots(const CGaugeGroup &GaugeGroup, int ith_simpleroot) const
{
	const size_t s0 = GaugeGroup.factor.size();
	size_t s1 = 0;
	unsigned i = 0;
	unsigned j = 0;

	size_t number_of_srs = 0;
	for (unsigned i = 0; i < s0; ++i)
		number_of_srs += GaugeGroup.factor[i].simpleroots.size();

	if (number_of_srs == 0)
	{
		(*this->out) << "  " << this->cbegin << "Number of simple roots is zero." << this->cend << "\n";
		return;
	}
	if ((ith_simpleroot >= 0) && (ith_simpleroot >= number_of_srs))
	{
		(*this->out) << "  " << this->cbegin << "Simple root #" << ith_simpleroot+1 << " does not exist. Number of simple roots: " << number_of_srs << this->cend << "\n";
		return;
	}

	CVector tmp;
	unsigned counter = 0;

	if (ith_simpleroot < 0)
	{
		(*this->out) << "  " << this->cbegin << "Simple roots:" << this->cend << "\n";
		if (this->output_type == Tmathematica)
			(*this->out) << "  lsr={\n";

		for (i = 0; i < s0; ++i)
		{
			const vector<vector<double> > &ggf_SimpleRoots = GaugeGroup.factor[i].simpleroots;
			s1 = ggf_SimpleRoots.size();

			for (j = 0; j < s1; ++j)
			{
				++counter;
				tmp = ggf_SimpleRoots[j];

				if (tmp.GetSize() != 16)
				{
					cout << "\n  Warning in void CPrint::PrintSimpleRoots(...) const: Simple root does not have 16 entires. Return." << endl;
					return;
				}

				(*this->out) << "    ";
				this->PrintRational(tmp, E8xE8);

				if (this->output_type == Tmathematica)
				{
					if (counter != number_of_srs)
						(*this->out) << ",\n";
					else
						(*this->out) << "\n  };\n";
				}
				else
					(*this->out) << "\n";
			}
		}
	}
	else
	{
		(*this->out) << "  " << this->cbegin << "simple root (" << ith_simpleroot+1 << "):" << this->cend << "\n    ";
		if (this->output_type == Tmathematica)
			(*this->out) << "lsr" << ith_simpleroot << "=";

		bool stop = false;
		for (i = 0; !stop && (i < s0); ++i)
		{
			const vector<vector<double> > &ggf_SimpleRoots = GaugeGroup.factor[i].simpleroots;
			s1 = ggf_SimpleRoots.size();

			for (j = 0; !stop && (j < s1); ++j)
			{
				if (counter == ith_simpleroot)
				{
					stop = true;
					tmp = ggf_SimpleRoots[j];
				}
				else
					++counter;
			}
		}

		if (tmp.GetSize() != 16)
		{
			cout << "\n  Warning in void CPrint::PrintSimpleRoots(...) const: Simple root does not have 16 entires. Return." << endl;
			return;
		}

		this->PrintRational(tmp, E8xE8);

		if (this->output_type == Tmathematica)
			(*this->out) << ";";
		(*this->out) << "\n";
	}
}



/* ########################################################################################
######   PrintSortedFieldLabels(const vector<unsigned> &FieldIndices, ...) const     ######
######                                                                               ######
######   Version: 04.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldIndices : indices of fields that shall be printed                   ######
######   2) VEVConfig    : contains all fields                                       ######
######   output:                                                                     ######
######   return value    : printed succesfully?                                      ######
######################################################################################## */
bool CPrint::PrintSortedFieldLabels(const vector<unsigned> &FieldIndices, const SConfig &VEVConfig) const
{
	const size_t s1 = FieldIndices.size();
	if (s1 == 0)
		return true;

	vector<string>            KnownLabels;
	vector<vector<unsigned> > KnownIndices;
	vector<unsigned>          tmp_KnownIndices;

	bool     LabelNotKnown = true;
	string   Label = "";
	unsigned Index = 0;

	unsigned i = 0;
	unsigned j = 0;
	size_t s2 = 0;
	size_t s3 = 0;

	for (i = 0; i < s1; ++i)
	{
		const CField &Field = VEVConfig.Fields[FieldIndices[i]];
		Label = Field.Labels[VEVConfig.use_Labels];
		Index = Field.Numbers[VEVConfig.use_Labels];

		LabelNotKnown = true;
		s2 = KnownLabels.size();
		for (j = 0; LabelNotKnown && (j < s2); ++j)
		{
			if (KnownLabels[j] == Label)
			{
				KnownIndices[j].push_back(Index);
				LabelNotKnown = false;
			}
		}
		if (LabelNotKnown)
		{
			KnownLabels.push_back(Label);
			tmp_KnownIndices.clear();
			tmp_KnownIndices.push_back(Index);
			KnownIndices.push_back(tmp_KnownIndices);
		}
	}

	string comma = "";
	if (this->output_type != Tstandard)
		comma = ",";

	s2 = KnownIndices.size();
	for (i = 0; i < s2; ++i)
	{
		Label = KnownLabels[i];
		vector<unsigned>  &tmp = KnownIndices[i];
		sort(tmp.begin(), tmp.end());

		s3 = tmp.size();
		for (j = 0; j < s3; ++j)
		{
			Index = tmp[j];

			(*this->out) << " " << this->prelabel << Label << this->start_index;
			if (output_type == Tlatex)
				(*this->out) << "{" << Index << "}";
			else
				(*this->out) << Index;

			if ((j+2 < s3) && (tmp[j+1] - tmp[j] == 1) && (tmp[j+2] - tmp[j+1] == 1))
			{
				(*this->out) << " - ";
				++j;
				while ((j < s3) && (tmp[j] - tmp[j-1] == 1))
				{
					++j;
				}
				--j;
				Index = tmp[j];
				(*this->out) << this->prelabel << Label << this->start_index;
				if (output_type == Tlatex)
					(*this->out) << "{" << Index << "}";
				else
					(*this->out) << Index;
			}
			if (j+1 < s3)
				(*this->out) << comma;
		}
	}
	return true;
}



/* ########################################################################################
######   PrintSpectrum(const CSpectrum &Spectrum) const                              ######
######                                                                               ######
######   Version: 22.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Spectrum  : CSpectrum object to be printed                               ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintSpectrum(const CSpectrum &Spectrum) const
{
	if (Spectrum.GetNumberOfSortedSpectra() == 0)
		this->PrintSSpectrum(Spectrum.GetSpectrum(), Spectrum.GetUseSpecialU1(), Spectrum.GetPosition_of_and_in_GaugeGroup());
	else
		this->PrintSSpectrum(Spectrum.GetSortedSpectrum(0), Spectrum.GetUseSpecialU1(), Spectrum.GetPosition_of_and_in_GaugeGroup());

	const vector<double> &AdditionalIdentifier = Spectrum.GetAdditionalIdentifier();

	size_t s1 = AdditionalIdentifier.size();
	if (s1 != 0)
	{
		(*this->out) << "\n  additional identifiers: ";
		for (unsigned i = 0; i < s1; ++i)
			(*this->out) << AdditionalIdentifier[i] << " ";
		(*this->out) << "\n";
	}
	return true;
}



/* ########################################################################################
######   PrintSpectrum(const vector<unsigned> &FieldIndices, ...) const              ######
######                                                                               ######
######   Version: 20.06.2013                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldIndices : indices of fields that shall be printed                   ######
######   2) VEVConfig    : contains all fields                                       ######
######   3) print_Labels : print the field labels?                                   ######
######   output:                                                                     ######
######   return value    : printed succesfully?                                      ######
######################################################################################## */
bool CPrint::PrintSpectrum(const vector<unsigned> &FieldIndices, const SConfig &VEVConfig, bool print_Labels) const
{
	bool SortByMSSMFields = false;
	if (print_Labels && (VEVConfig.ConfigLabel == "SMConfig"))
		SortByMSSMFields = true;

	const vector<CField> &Fields = VEVConfig.Fields;

	const size_t f1 = Fields.size();
	const size_t f2 = FieldIndices.size();

	if (f2 == 0)
	{
		(*this->out) << "  " << this->cbegin << "Spectrum is empty." << this->cend << endl;
		return true;
	}

	const SSymmetryGroup &SymmetryGroup = VEVConfig.SymmetryGroup;

	const bool print_U1Charges = (SymmetryGroup.observable_sector_U1s.size() != 0);
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	vector<unsigned> tmp_FieldIndices;

	vector<SUSYMultiplet>     Spec_Multiplets;
	vector<unsigned>          Spec_Multiplicities;
	vector<RepVector>         Spec_Dimensions;
	vector<CVector>           Spec_U1Charges;
	vector<vector<unsigned> > Spec_FieldIndices;

	size_t s1 = 0;
	bool field_not_known = true;

	unsigned Index = 0;
	for (i = 0; i < f2; ++i)
	{
		Index = FieldIndices[i];
		if (Index >= f1)
		{
			cout << "\n  Warning in bool CPrint::PrintSpectrum(..) const: FieldIndex out of range. Return false." << endl;
			return false;
		}
		const CField &Field = Fields[Index];

		field_not_known = true;
		s1 = Spec_Multiplicities.size();
		for (j = 0; field_not_known && (j < s1); ++j)
		{
			if ((Field.Multiplet == Spec_Multiplets[j]) && AreRepVectorsEqual(SymmetryGroup, Spec_Dimensions[j], Field.Dimensions) && AreU1ChargesEqual(SymmetryGroup, Spec_U1Charges[j],Field.U1Charges))
			{
				Spec_Multiplicities[j] += Multiplicity(SymmetryGroup, Field.Dimensions);
				if (print_Labels)
					Spec_FieldIndices[j].push_back(Index);

				field_not_known = false;
			}
		}
		if (field_not_known)
		{
			Spec_Multiplets.push_back(Field.Multiplet);
			Spec_Multiplicities.push_back(Multiplicity(SymmetryGroup, Field.Dimensions));
			Spec_Dimensions.push_back(Field.Dimensions);
			Spec_U1Charges.push_back(Field.U1Charges);

			tmp_FieldIndices.clear();
			if (print_Labels)
				tmp_FieldIndices.push_back(Index);
			Spec_FieldIndices.push_back(tmp_FieldIndices);
		}
	}

	SDimension One;
	One.Dimension = 1;
	One.AdditionalLabel = "";

	vector<unsigned> ExtraLines;
	// begin: for SM configurations sort the table
	if (SortByMSSMFields)
	{
		vector<SUSYMultiplet>     tmp_Spec_Multiplets;
		vector<unsigned>          tmp_Spec_Multiplicities;
		vector<RepVector>         tmp_Spec_Dimensions;
		vector<CVector>           tmp_Spec_U1Charges;
		vector<vector<unsigned> > tmp_Spec_FieldIndices;

		const unsigned use_Labels = VEVConfig.use_Labels;
		s1 = Spec_Multiplicities.size();

		vector<bool> Spec_Found(s1, false);

		vector<string> Order;
		Order.push_back("n");
		Order.push_back("q");    Order.push_back("bq");
		Order.push_back("bu");   Order.push_back("u");
		Order.push_back("bd");   Order.push_back("d");
		Order.push_back("l");    Order.push_back("bl");
		Order.push_back("be");   Order.push_back("e");

		Order.push_back("v");    Order.push_back("bv");
		Order.push_back("w");    Order.push_back("bw");
		Order.push_back("x");    Order.push_back("bx");
		Order.push_back("y");    Order.push_back("by");
		Order.push_back("z");    Order.push_back("bz");

		bool FirstNonMSSMFound = false;
		for (i = 0; i < Order.size(); ++i)
		{
			for (j = 0; j < s1; ++j)
			{
				if ((Spec_Found[j] == false) && (Fields[Spec_FieldIndices[j][0]].Labels[use_Labels] == Order[i]))
				{
					if (i == 0)
						ExtraLines.push_back(0);

					if ((i > 10) && !FirstNonMSSMFound)
					{
						ExtraLines.push_back(tmp_Spec_Multiplets.size()-1);
						FirstNonMSSMFound = true;
					}

					Spec_Found[j] = true;
					tmp_Spec_Multiplets.push_back(Spec_Multiplets[j]);
					tmp_Spec_Multiplicities.push_back(Spec_Multiplicities[j]);
					tmp_Spec_Dimensions.push_back(Spec_Dimensions[j]);
					tmp_Spec_U1Charges.push_back(Spec_U1Charges[j]);
					tmp_Spec_FieldIndices.push_back(Spec_FieldIndices[j]);
					break;
				}
			}
		}
		FirstNonMSSMFound = false;
		for (j = 0; j < s1; ++j)
		{
			if (Spec_Found[j] == false)
			{
				if (!FirstNonMSSMFound)
				{
					ExtraLines.push_back(tmp_Spec_Multiplets.size()-1);
					FirstNonMSSMFound = true;
				}

				Spec_Found[j] = true;
				tmp_Spec_Multiplets.push_back(Spec_Multiplets[j]);
				tmp_Spec_Multiplicities.push_back(Spec_Multiplicities[j]);
				tmp_Spec_Dimensions.push_back(Spec_Dimensions[j]);
				tmp_Spec_U1Charges.push_back(Spec_U1Charges[j]);
				tmp_Spec_FieldIndices.push_back(Spec_FieldIndices[j]);
			}
		}

		Spec_Multiplets     = tmp_Spec_Multiplets;
		Spec_Multiplicities = tmp_Spec_Multiplicities;
		Spec_Dimensions     = tmp_Spec_Dimensions;
		Spec_U1Charges      = tmp_Spec_U1Charges;
		Spec_FieldIndices   = tmp_Spec_FieldIndices;
	}
	// end: for SM configurations sort the table

	s1 = Spec_Multiplicities.size();
	const size_t s2 = Spec_Dimensions.at(0).size();
	const bool print_irrep = (s2 != 0);

	size_t s3 = 0;

	string table_seperator = " ";
	string table_eol = " ";
	string mathmode = "";
	string mathematicacomma = "";

	if (this->output_type == Tmathematica)
	{
		mathematicacomma = ",";
		table_seperator = ", ";
	}
	else
		if (this->output_type == Tlatex)
		{
			table_seperator = " & ";
			table_eol = "\\\\";
			mathmode = "$";

			string tabular_parameter = "|r|";
			string headline          = " \\#";

			if (print_irrep)
			{
				tabular_parameter += "c|";
				headline          += " & irrep";
			}
			if (print_U1Charges)
			{
				tabular_parameter += "l|";
				headline          += " & $\\text{U}(1)$ charges";
			}
			if (print_Labels)
			{
				tabular_parameter += "l|";
				headline          += " & labels ";
			}
			headline += "\\\\\n";

			(*this->out) << "\\begin{longtable}{" << tabular_parameter << "}\n";
			(*this->out) << "\\hline\n" << headline << "\\hline\\hline\n";
			(*this->out) << "\\endhead\n";
			(*this->out) << "\\hline\n";
			(*this->out) << "\\endfoot\n";
		}

	vector<SUSYMultiplet> OrderOfSUSYMultiplets;

	// Non-SuperSymmetric types
	OrderOfSUSYMultiplets.push_back(Gauge);			//hacking here!!!
	OrderOfSUSYMultiplets.push_back(bGauge);
	OrderOfSUSYMultiplets.push_back(Scalar);
	OrderOfSUSYMultiplets.push_back(bScalar);
	OrderOfSUSYMultiplets.push_back(LeftFermi);
	OrderOfSUSYMultiplets.push_back(RightFermi);

	// SuperSymmetric types
	OrderOfSUSYMultiplets.push_back(Gravity);
	OrderOfSUSYMultiplets.push_back(GravityCC);
	OrderOfSUSYMultiplets.push_back(moduli);
	OrderOfSUSYMultiplets.push_back(bmoduli);
	OrderOfSUSYMultiplets.push_back(LCModulus);
	OrderOfSUSYMultiplets.push_back(RCModulus);
	OrderOfSUSYMultiplets.push_back(Vector);
	OrderOfSUSYMultiplets.push_back(VectorCC);
	OrderOfSUSYMultiplets.push_back(LeftChiral);
	OrderOfSUSYMultiplets.push_back(RightChiral);
	OrderOfSUSYMultiplets.push_back(Hyper);
	OrderOfSUSYMultiplets.push_back(Halfhyper);

	bool firsttime = false;
	bool SUSYMultipletPrinted = false;

	if (this->output_type == Tmathematica)
		(*this->out) << "  lstspectrum = " << this->set_open << "\n";

	for (i = 0; i < OrderOfSUSYMultiplets.size(); ++i)
	{
		const SUSYMultiplet &Multiplet = OrderOfSUSYMultiplets[i];

		firsttime = true;
		for (j = 0; j < s1; ++j)
		{
			if (Multiplet == Spec_Multiplets[j])
			{
				if (firsttime && SUSYMultipletPrinted && (this->output_type == Tstandard))
				{
					SUSYMultipletPrinted = false;
					(*this->out) << "\n";
				}
				if (firsttime)
				{
					SUSYMultipletPrinted = true;
					firsttime = false;
				}

				(*this->out) << "  " << this->set_open << setw(3) << Spec_Multiplicities[j] << table_seperator;

				(*this->out) << mathmode;
				if (print_irrep)
					this->PrintRep(Spec_Dimensions[j], SymmetryGroup);
				else
					this->PrintSDimension(One);
				this->PrintSUSYType(Multiplet);
				(*this->out) << mathmode;

				if (print_U1Charges)
				{
					(*this->out) << table_seperator;

					if (this->output_type == Tstandard)
						(*this->out) <<" U(1) : ";

					(*this->out) << mathmode;
					this->PrintU1Charges(Spec_U1Charges[j], SymmetryGroup);
					(*this->out) << mathmode;
				}

				if (print_Labels)
				{
					(*this->out) << table_seperator;

					vector<unsigned> &tmp_FieldIndices = Spec_FieldIndices[j];

					if (this->output_type == Tmathematica)
					{
						(*this->out) << this->set_open;
						s3 = tmp_FieldIndices.size();
						for (k = 0; k < s3; ++k)
						{
							(*this->out) << " ";
							this->PrintLabel(Fields[tmp_FieldIndices[k]], VEVConfig.use_Labels);
							if (k+1 < s3)
								(*this->out) << mathematicacomma;
						}

						(*this->out) << this->set_close;
					}
					else
					{
						(*this->out) << mathmode;
						this->PrintSortedFieldLabels(tmp_FieldIndices, VEVConfig);
						(*this->out) << mathmode;
					}
				}

				(*this->out) << table_eol << this->set_close;

				if (j+1 < s1)
					(*this->out) << mathematicacomma;

				(*this->out) << "\n";

				if ((this->output_type == Tstandard) && (ExtraLines.size() != 0))
				{
					if (find(ExtraLines.begin(), ExtraLines.end(), j) != ExtraLines.end())
						(*this->out) << "\n";
				}
			}
		}
	}

	if (this->output_type == Tmathematica)
		(*this->out) << "  " << this->set_close << this->endofset << "\n";
	else
		if (this->output_type == Tlatex)
			(*this->out) << "\\end{longtable}" << "\n";

	(*this->out) << flush;

	return true;
}



/* ########################################################################################
######   PrintState(const CState &State, const CSector &Sector, ...) const           ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) State        : CState object to be printed                               ######
######   2) Sector       : the corresponding sector                                  ######
######   3) FixedBrane   : the corresponding fixed point/ brane                      ######
######   4) VEVConfig    : contains the fields                                       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintState(const CState &State, const CSector &Sector, const CFixedBrane &FixedBrane, const SConfig &VEVConfig) const
{
	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;

	bool print_Labels = true;
	bool print_Representations = true;

	vector<unsigned> FieldIndices;
	State.GetFieldIndices(VEVConfig.Fields, AnyKind, FieldIndices);
	if (FieldIndices.size() == 0)
		return;
	//print_Representations = false;

	// begin: print localization of the state
	const CSpaceGroupElement &SGElement = FixedBrane.GetSGElement();

	(*this->out) << "=========================================================================\n";
	if ((SGElement.Get_m() == 0) && (SGElement.Get_n() == 0) && (SGElement.Get_k() == 0))
		(*this->out) << "  untwisted Sector\n";
	else
	{
		(*this->out) << "  (" << SGElement.Get_m() << ", " << SGElement.Get_n() << ") twisted Sector ";
		this->PrintSGElement(SGElement, true);
		(*this->out) << "\n";
	}
	if (FixedBrane.GetFixedBraneLabel() != "")
		(*this->out) << "  label: " << FixedBrane.GetFixedBraneLabel() << "\n";
	// end: print localization of the state

	// begin: print list of representations
	if (print_Representations)
	{
		(*this->out) << "-------------------------------------------------------------------------\n";
		this->PrintSpectrum(FieldIndices, VEVConfig, print_Labels);
	}
	// end: print list of representations

	const CMasslessHalfState &MLeftMover  = FixedBrane.GetMasslessLeftMover(State.GetLeftMover().GetIndex());
	//const CMasslessHalfState &MRightMover = FixedBrane.MasslessRightMovers[State.RightMover.GetIndex()];

	(*this->out) << "-------------------------------------------------------------------------\n";

	if (!print_Representations)
		this->PrintHalfState(State.GetLeftMover(), MLeftMover, Sector.GetLM_all_Oscillators());
	else
	{
		(*this->out) << "  " << this->cbegin << "  left-movers: ";

		// begin: print oscillators of left mover
		if (MLeftMover.Excitation.OscillatorIndices.size() != 0)
		{
			(*this->out) << "N = " << MLeftMover.Excitation.NumberOperator << " :";

			const vector<unsigned> &OscillatorIndices = MLeftMover.Excitation.OscillatorIndices;
			s1 = OscillatorIndices.size();
			for (j = 0; j < s1; ++j)
			{
				(*this->out) << " ";
				this->PrintOscillator(Sector.GetLM_Oscillator(OscillatorIndices[j]));
			}
			(*this->out) << " ";
		}
		// end: print oscillators

		// begin: print weights of left mover
		(*this->out) << "#(weights) = " << State.GetLeftMover().Weights.size() << this->cend << "\n";

		SelfDualLattice Lattice = E8xE8;
		s1 = FieldIndices.size();
		for (i = 0; i < s1; ++i)
		{
			const CField &Field = VEVConfig.Fields[FieldIndices[i]];

			const size_t s3 = Field.WeightIndices.size();
			for (j = 0; j < s3; ++j)
			{
				(*this->out) << "    ";
				this->PrintRational(MLeftMover.Weights[Field.WeightIndices[j]], Lattice);
				(*this->out) << ", ";
				this->PrintRep(Field, VEVConfig.SymmetryGroup, false);
				if (print_Labels)
				{
					(*this->out) << "  ";
					this->PrintLabel(Field, VEVConfig.use_Labels);
				}
				(*this->out) << "\n";
			}
		}
		// end: print weights of left mover
	}
	(*this->out) << "-------------------------------------------------------------------------\n";
	this->PrintHalfState(State.GetRightMover(), Sector.GetMasslessRightMover(State.GetRightMover().GetIndex()), Sector.GetRM_all_Oscillators());

	// begin: print elements of centralizer and corresponding eigenvalues
	//(*this->out) << "------------------------------------------------\n";
	//(*this->out) << " Eigenvalues:\n\n centralizer \t\t\t LM \t RM\n";
	/*const size_t s4 = FixedBrane.Centralizer.size();
  for (i = 0; i < s4; ++i)
  {
  const vector<int> &Cent_Label = FixedBrane.Centralizer.at(i).Element_Label;
  (*this->out) << "(" << Cent_Label[0] << ", " << Cent_Label[1] << ") ("
  << Cent_Label[2] << ", " << Cent_Label[3] << ", "
  << Cent_Label[4] << ", " << Cent_Label[5] << ", "
  << Cent_Label[6] << ", " << Cent_Label[7] << ")";

  (*this->out) << "\t" << this->LeftMover.Eigenvalues[i] << "\t" << this->RightMover.Eigenvalues[i] << "\n";
}*/
	// end: print elements of centralizer and corresponding eigenvalues

	// begin: print the gamma phases
	const size_t s5 = State.GetNumberOfGammaPhases();
	if (s5 != 0)
	{
		(*this->out) << "-------------------------------------------------------------------------\n";
		for (i = 0; i < s5; ++i)
			(*this->out) << "  gamma_" << i+1 << " = " << State.GetGammaPhase(i) << "\n";
	}
	// end: print the gamma phases

	(*this->out) << "=========================================================================\n" << endl;
}



/* ########################################################################################
######   PrintState(const CState &State, const CSector &Sector, ...) const           ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) State        : CState object to be printed                               ######
######   2) Sector       : the corresponding sector                                  ######
######   3) FixedBrane   : the corresponding fixed point/ brane                      ######
######   4) VEVConfig    : contains the fields                                       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintTachyonicState(const CState &State, const CSector &Sector, const CFixedBrane &FixedBrane, const SConfig &VEVConfig, const bool &Rexcited) const
{
	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;

	bool print_Labels = true;
	bool print_Representations = true;

	vector<unsigned> FieldIndices;
	State.GetFieldIndices(VEVConfig.Fields, AnyKind, FieldIndices);
	if (FieldIndices.size() == 0)
		return;
	//print_Representations = false;

	// begin: print localization of the state
	const CSpaceGroupElement &SGElement = FixedBrane.GetSGElement();

	(*this->out) << "=========================================================================\n";
	if ((SGElement.Get_m() == 0) && (SGElement.Get_n() == 0) && (SGElement.Get_k() == 0))
		(*this->out) << "  untwisted Sector\n";
	else
	{
		(*this->out) << "  (" << SGElement.Get_m() << ", " << SGElement.Get_n() << ") twisted Sector ";
		this->PrintSGElement(SGElement, true);
		(*this->out) << "\n";
	}
	if (FixedBrane.GetFixedBraneLabel() != "")
		(*this->out) << "  label: " << FixedBrane.GetFixedBraneLabel() << "\n";
	// end: print localization of the state

	// begin: print list of representations
	if (print_Representations)
	{
		(*this->out) << "-------------------------------------------------------------------------\n";
		this->PrintSpectrum(FieldIndices, VEVConfig, print_Labels);
	}
	// end: print list of representations


	(*this->out) << "-------------------------------------------------------------------------\n";

	if (!Rexcited)
	{
		const CMasslessHalfState &MLeftMover  = FixedBrane.TachyonicGetMassLeftMover(State.GetLeftMover().GetIndex());

		if (!print_Representations)
			this->PrintHalfState(State.GetLeftMover(), MLeftMover, Sector.GetLM_all_Oscillators());
		else
		{
			(*this->out) << "  " << this->cbegin << "  left-movers: ";

			// begin: print oscillators of left mover
			if (MLeftMover.Excitation.OscillatorIndices.size() != 0)
			{
				(*this->out) << "N = " << MLeftMover.Excitation.NumberOperator << " :";

				const vector<unsigned> &OscillatorIndices = MLeftMover.Excitation.OscillatorIndices;
				s1 = OscillatorIndices.size();
				for (j = 0; j < s1; ++j)
				{
					(*this->out) << " ";
					this->PrintOscillator(Sector.GetLM_Oscillator(OscillatorIndices[j]));
				}
				(*this->out) << " ";
			}
			// end: print oscillators

			// begin: print weights of left mover
			(*this->out) << "#(weights) = " << State.GetLeftMover().Weights.size() << this->cend << "\n";

			SelfDualLattice Lattice = E8xE8;
			s1 = FieldIndices.size();
			for (i = 0; i < s1; ++i)
			{
				const CField &Field = VEVConfig.Fields[FieldIndices[i]];

				const size_t s3 = Field.WeightIndices.size();
				for (j = 0; j < s3; ++j)
				{
					(*this->out) << "    ";
					this->PrintRational(MLeftMover.Weights[Field.WeightIndices[j]], Lattice);
					(*this->out) << ", ";
					this->PrintRep(Field, VEVConfig.SymmetryGroup, false);
					if (print_Labels)
					{
						(*this->out) << "  ";
						this->PrintLabel(Field, VEVConfig.use_Labels);
					}
					(*this->out) << "\n";
				}
			}
			// end: print weights of left mover
		}
	}
	else
	{
		const CMasslessHalfState &MLeftMover  = FixedBrane.ExcitedTachyonicGetMassLeftMover(State.GetLeftMover().GetIndex());

		if (!print_Representations)
			this->PrintHalfState(State.GetLeftMover(), MLeftMover, Sector.GetLM_all_Oscillators());
		else
		{
			(*this->out) << "  " << this->cbegin << "  left-movers: ";

			// begin: print oscillators of left mover
			if (MLeftMover.Excitation.OscillatorIndices.size() != 0)
			{
				(*this->out) << "N = " << MLeftMover.Excitation.NumberOperator << " :";

				const vector<unsigned> &OscillatorIndices = MLeftMover.Excitation.OscillatorIndices;
				s1 = OscillatorIndices.size();
				for (j = 0; j < s1; ++j)
				{
					(*this->out) << " ";
					this->PrintOscillator(Sector.GetLM_Oscillator(OscillatorIndices[j]));
				}
				(*this->out) << " ";
			}
			// end: print oscillators

			// begin: print weights of left mover
			(*this->out) << "#(weights) = " << State.GetLeftMover().Weights.size() << this->cend << "\n";

			SelfDualLattice Lattice = E8xE8;
			s1 = FieldIndices.size();
			for (i = 0; i < s1; ++i)
			{
				const CField &Field = VEVConfig.Fields[FieldIndices[i]];

				const size_t s3 = Field.WeightIndices.size();
				for (j = 0; j < s3; ++j)
				{
					(*this->out) << "    ";
					this->PrintRational(MLeftMover.Weights[Field.WeightIndices[j]], Lattice);
					(*this->out) << ", ";
					this->PrintRep(Field, VEVConfig.SymmetryGroup, false);
					if (print_Labels)
					{
						(*this->out) << "  ";
						this->PrintLabel(Field, VEVConfig.use_Labels);
					}
					(*this->out) << "\n";
				}
			}
			// end: print weights of left mover
		}
	}

	(*this->out) << "-------------------------------------------------------------------------\n";

	//Assuming at most one excited right-moving tachyonic level per sector
	if (!Rexcited)
		this->PrintTachyonHalfState(State.GetRightMover(), Sector.GetRTachyons()[0], Sector.GetRM_all_Oscillators());
	else
		this->PrintTachyonHalfState(State.GetRightMover(), Sector.GetRTachyons()[1], Sector.GetRM_all_Oscillators());


	// begin: print elements of centralizer and corresponding eigenvalues
	//(*this->out) << "------------------------------------------------\n";
	//(*this->out) << " Eigenvalues:\n\n centralizer \t\t\t LM \t RM\n";
	/*const size_t s4 = FixedBrane.Centralizer.size();
  for (i = 0; i < s4; ++i)
  {
  const vector<int> &Cent_Label = FixedBrane.Centralizer.at(i).Element_Label;
  (*this->out) << "(" << Cent_Label[0] << ", " << Cent_Label[1] << ") ("
  << Cent_Label[2] << ", " << Cent_Label[3] << ", "
  << Cent_Label[4] << ", " << Cent_Label[5] << ", "
  << Cent_Label[6] << ", " << Cent_Label[7] << ")";

  (*this->out) << "\t" << this->LeftMover.Eigenvalues[i] << "\t" << this->RightMover.Eigenvalues[i] << "\n";
}*/
	// end: print elements of centralizer and corresponding eigenvalues

	// begin: print the gamma phases
	const size_t s5 = State.GetNumberOfGammaPhases();
	if (s5 != 0)
	{
		(*this->out) << "-------------------------------------------------------------------------\n";
		for (i = 0; i < s5; ++i)
			(*this->out) << "  gamma_" << i+1 << " = " << State.GetGammaPhase(i) << "\n";
	}
	// end: print the gamma phases

	(*this->out) << "=========================================================================\n" << endl;
}




/* ########################################################################################
######   PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const      ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold  : all states of this orbifold shall be printed                 ######
######   2) VEVConfig : contains the fields                                          ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const
{
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const size_t s1 = Orbifold.GetNumberOfSectors();
	size_t s2 = 0;
	size_t s3 = 0;

	for (i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
		{
			const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

			s3 = FixedBrane.GetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
				this->PrintState(FixedBrane.GetInvariantState(k), Sector, FixedBrane, VEVConfig);
		}
	}
	(*this->out) << flush;

	return true;
}



/* ########################################################################################
######   PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const      ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold  : all states of this orbifold shall be printed                 ######
######   2) VEVConfig : contains the fields                                          ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintTachyonicStates(const COrbifold &Orbifold, const SConfig &VEVConfig) const
{
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const size_t s1 = Orbifold.GetNumberOfSectors();
	size_t s2 = 0;
	size_t s3 = 0;

	for (i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
		{
			const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

			s3 = FixedBrane.TachyonicGetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
				this->PrintTachyonicState(FixedBrane.TachyonicGetInvariantState(k), Sector, FixedBrane, VEVConfig, false);

			s3 = FixedBrane.ExcitedTachyonicGetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
				this->PrintTachyonicState(FixedBrane.ExcitedTachyonicGetInvariantState(k), Sector, FixedBrane, VEVConfig, true);
		}
	}
	(*this->out) << flush;

	return true;
}



/* ########################################################################################
######   PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig, ...) const ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold                 : the corresponding orbifold                    ######
######   2) VEVConfig                : contains all fields                           ######
######   3) FieldIndices             : indices of the fields to be printed           ######
######   4) PrintInternalInformation : print gamma phases and internal indices?      ######
######   output:                                                                     ######
######   return value                : printed succesfully?                          ######
######################################################################################## */
bool CPrint::PrintStates(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, bool PrintInternalInformation) const
{
	const vector<CField> &SV_Fields = Orbifold.StandardConfig.Fields;
	const vector<CField> &Fields    = VEVConfig.Fields;
	
	const size_t f1 = Fields.size();

	string tmp_string1 = "";
	string tmp_string2 = "";

	unsigned i = 0;
	unsigned j = 0;

	const size_t t1 = FieldIndices.size();
	size_t t2 = 0;

	const SSymmetryGroup &St_SymmetryGroup = Orbifold.StandardConfig.SymmetryGroup;

	const bool PrintU1_VEVConfig      = (VEVConfig.SymmetryGroup.observable_sector_U1s.size() != 0);
	const bool PrintU1_StandardConfig = (St_SymmetryGroup.observable_sector_U1s.size() != 0);

	bool Print_VEVConfig = ((VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size() != St_SymmetryGroup.GaugeGroup.u1directions.size())
			|| (VEVConfig.SymmetryGroup.GaugeGroup.factor.size()       != St_SymmetryGroup.GaugeGroup.factor.size())
			|| (VEVConfig.SymmetryGroup.observable_sector_GGs          != St_SymmetryGroup.observable_sector_GGs)
			|| (VEVConfig.SymmetryGroup.observable_sector_U1s          != St_SymmetryGroup.observable_sector_U1s));
	for (i = 0; !Print_VEVConfig && (i < VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size()); ++i)
	{
		if (VEVConfig.SymmetryGroup.GaugeGroup.u1directions[i] != St_SymmetryGroup.GaugeGroup.u1directions[i])
			Print_VEVConfig = true;
	}

	const vector<SDiscreteSymmetry> &DiscreteNonRSymmetries = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteNonRSymmetries;
	const vector<SDiscreteSymmetry> &DiscreteRSymmetries    = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteRSymmetries;
	const vector<SModularSymmetry>  &ModularSymmetries      = Orbifold.OrbifoldGroup.GetSpaceGroup().ModularSymmetries;

	const bool Print_DiscreteNonRSymmetries = (DiscreteNonRSymmetries.size() != 0);
	const bool Print_DiscreteRSymmetries    = (DiscreteRSymmetries.size() != 0);
	const bool Print_ModularSymmetries      = (ModularSymmetries.size() != 0);

	const vector<CSector> &Sectors = Orbifold.GetSectors();
	const SelfDualLattice &Lattice = Orbifold.OrbifoldGroup.GetLattice();

	for (i = 0; i < t1; ++i)
	{
		if (FieldIndices[i] >= f1)
		{
			cout << "\n  Warning in bool CPrint::PrintStates(...) const: \"FieldIndex\" out of range. Return false." << endl;
			return false;
		}
		const CField &tmp_Field = Fields[FieldIndices[i]];

		(*this->out) << "\n    ";
		if (tmp_Field.Multiplet == RightChiral)
			(*this->out) << "CPT partner: ";
		else
			if ((tmp_Field.Multiplet == LCModulus) || (tmp_Field.Multiplet == RCModulus))
				(*this->out) << "modulus: ";
			else
				if ((tmp_Field.Multiplet == Vector) || (tmp_Field.Multiplet == VectorCC))
					(*this->out) << "gauge boson: ";
				else
					if ((tmp_Field.Multiplet == Gravity) || (tmp_Field.Multiplet == GravityCC))
						(*this->out) << "graviton: ";

		this->PrintLabel(tmp_Field, VEVConfig.use_Labels);

		if ((VEVConfig.use_Labels != 0) && (tmp_Field.Labels[0] != tmp_Field.Labels[VEVConfig.use_Labels]) && (tmp_Field.Numbers[0] != tmp_Field.Numbers[VEVConfig.use_Labels]))
		{
			(*this->out) << "  " << this->vector_open;
			this->PrintLabel(tmp_Field, 0);
			(*this->out) << this->vector_close;
		}
		(*this->out) << "\n";

		const CSpaceGroupElement &Element = tmp_Field.SGElement;
		(*this->out) << "  sector (m,n,k)        : " << this->vector_open << Element.Get_m() << this->separator << Element.Get_n() << this->separator << Element.Get_k() << this->vector_close << "\n";
		(*this->out) << "  fixed point n_a     : " << this->vector_open;
		this->PrintRational(Element.Get_n(0), false);
		(*this->out) << this->separator;
		this->PrintRational(Element.Get_n(1), false);
		(*this->out) << this->separator;
		this->PrintRational(Element.Get_n(2), false);
		(*this->out) << this->separator;
		this->PrintRational(Element.Get_n(3), false);
		(*this->out) << this->separator;
		this->PrintRational(Element.Get_n(4), false);
		(*this->out) << this->separator;
		this->PrintRational(Element.Get_n(5), false);
		(*this->out) << this->vector_close << "\n";

		// begin: print discrete space group charges
		if (Print_DiscreteNonRSymmetries)
		{
			(*this->out) << "  space group charges : " << this->vector_open;
			t2 = DiscreteNonRSymmetries.size();
			for (j = 0; j < t2; ++j)
			{
				this->PrintRational(tmp_Field.GetDiscreteCharge(DiscreteNonRSymmetries[j]), false);
				if (j+1 < t2)
					(*this->out) << this->separator;
			}
			(*this->out) << this->vector_close << "\n";
		}
		// end: print discrete space group charges

		(*this->out) << "\n";

		// begin: print (gauge and acc.) charges
		if (Print_VEVConfig)
		{
			(*this->out) << "  rep. in config      : ";
			this->PrintRep(tmp_Field, VEVConfig.SymmetryGroup, PrintU1_VEVConfig);
			(*this->out) << "\n";
		}
		if (tmp_Field.OriginStandardConfig >= 0)
		{
			const unsigned OrigIndex = (unsigned)abs(tmp_Field.OriginStandardConfig);

			(*this->out) << "  representation      : ";
			this->PrintRep(SV_Fields[OrigIndex], St_SymmetryGroup, PrintU1_StandardConfig);

			if (OrigIndex != FieldIndices[i])
			{

				(*this->out) << "  with label ";
				this->PrintLabel(SV_Fields[OrigIndex], Orbifold.StandardConfig.use_Labels);
				(*this->out) << " in orginal vacuum";
			}
			(*this->out) << "\n";
		}

		if (VEVConfig.SymmetryGroup.BmLGenerator.GetSize() == 16)
		{
			(*this->out) << "  U(1)_B-L charge     : ";
			this->PrintRational(D2Rat(tmp_Field.BmLCharge), false);
			(*this->out)<< "\n";
		}

		t2 = tmp_Field.AccU1Charges.size();
		if (t2 != 0)
		{
			(*this->out) << "  acc. U(1) charges   : " << this->vector_open;
			for (j = 0; j < t2; ++j)
				(*this->out) << tmp_Field.AccU1Charges[j] << this->separator;
			(*this->out) << this->vector_close << "\n";
		}
		// end: print (gauge and acc.) charges

		(*this->out) << "\n";
        (*this->out) << "  left-moving p_sh    : ";		
		(*this->out) << "\n";

        size_t t3 = 0;
        t3 = tmp_Field.GetNumberOfLMWeights();
        for (int k = 0; k < t3; ++k)
        {
          this->PrintRational(tmp_Field.GetLMWeight(k, Sectors), Lattice, false);
         (*this->out) << "\n";
        }
        (*this->out) << "\n";
		(*this->out) << "  right-moving q_sh   : ";
		this->PrintRational(tmp_Field.GetRMWeight(0, Sectors), SO8, false); 
		(*this->out) << "\n";

		// begin: print oscillators
		(*this->out) << "  oscillators         : ";
		t2 = tmp_Field.GetNumberOfOscillators(Sectors);
		if (t2 == 0)
			(*this->out) << "no\n";
		else
		{
			for (j = 0; j < t2; ++j)
			{
				this->PrintOscillator(tmp_Field.GetOscillator(j, Sectors));
				(*this->out) << " ";
			}
			(*this->out) << "\n";
		}		
		// end: print oscillators

		(*this->out) << "\n";

		(*this->out) << "  vev                 : " << tmp_Field.VEVs.GetLength() << "\n";

		if (PrintInternalInformation)
		{
			const vector<unsigned> &internalIndex = tmp_Field.GetInternalIndex();
			(*this->out) << "\n";
			// begin: print gamma
			const CState &State = tmp_Field.GetState(Sectors);
			t2 = State.GetNumberOfGammaPhases();
			if (t2 != 0)
			{
				(*this->out) << "  gamma phases        : " << this->vector_open;
				for (j = 0; j < t2; ++j)
				{
					this->PrintRational(D2Rat(State.GetGammaPhase(j)), false);
					if (j+1 < t2)
						(*this->out) << this->separator;
				}
				(*this->out) << this->vector_close << "\n";
			}
			// end: print gamma
			(*this->out) << "  internalIndex       : " << this->vector_open << internalIndex[0] << this->separator << internalIndex[1] << this->separator << internalIndex[2] << this->vector_close << "\n";
			(*this->out) << "  field no.           : " << FieldIndices[i] << "\n";
		}

		if (i + 1 < t1)
			(*this->out) << endl;
	}
	(*this->out) << endl;
	return true;
}


/* ########################################################################################		//hacking here!!!
######   PrintSummaryOfVEVConfig(const SConfig &VEVConfig, ...) const                ######
######                                                                               ######
######   Version: 20.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig    : contains the fields to be printed                         ######
######   2) Multiplet    : print only those fields of "VEVConfig" of this SUSY type  ######
######   3) print_Labels : print the field labels?                                   ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintSummaryOfVEVConfig(const SConfig &VEVConfig, const vector<SUSYMultiplet>  &Multiplet, bool print_Labels) const
{
	this->PrintGaugeGroup(VEVConfig);

	const vector<CField> &Fields = VEVConfig.Fields;
	const size_t f1 = Fields.size();

	//const bool UseAnySUSYKind = (Multiplet == AnyKind);

	vector<unsigned> FieldIndices;
	for (unsigned i = 0; i < f1; ++i)
	{
		for (int j=0; j<Multiplet.size(); j++)
		{
			if ( Fields[i].Multiplet == Multiplet[j] )
				FieldIndices.push_back(i);
		}
	}
	this->PrintSpectrum(FieldIndices, VEVConfig, print_Labels);
}



/* ########################################################################################
######   PrintSummary(const CSector &Sector, const SConfig &VEVConfig, ...) const    ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector       : fields from this CSector object shall be printed          ######
######   2) VEVConfig    : contains the fields                                       ######
######   3) Multiplet    : print only those fields of SUSY type "Multiplet"          ######
######   4) print_Labels : print the field labels?                                   ######
######   output:                                                                     ######
######   return value    : printed succesfully?                                      ######
######################################################################################## */
bool CPrint::PrintSummary(const CSector &Sector, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet, bool print_Labels) const 
{
	if ((Sector.Get_m() == 0) && (Sector.Get_n() == 0) && (Sector.Get_k() == 0))
		(*this->out) << "  " << this->cbegin << "U";
	else
		(*this->out) << "  " << this->cbegin << "T(" << Sector.Get_m() << "," << Sector.Get_n() << "," << Sector.Get_k() << ")";
	(*this->out) << " Sector:" << this->cend;

	vector<unsigned> FieldIndices;
	Sector.GetFieldIndices(VEVConfig.Fields, Multiplet, FieldIndices);

	if (FieldIndices.size() == 0)
	{
		(*this->out) << "  " << this->cbegin << " empty" << this->cend << "\n";
		return true;
	}
	(*this->out) << "\n";

	this->PrintSpectrum(FieldIndices, VEVConfig, print_Labels);

	return true;
}


/* ########################################################################################
######   PrintSummary(const CFixedBrane &FixedBrane, ...) const                      ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane    : fields from this CFixedBrane object shall be printed     ######
######   2) OrbifoldGroup : contains the local shift of this "FixedBrane"            ######
######   3) VEVConfig     : contains the fields                                      ######
######   4) Multiplet     : print only those fields of SUSY type "Multiplet"         ######
######   5) print_Labels  : print the field labels?                                  ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
bool CPrint::PrintSummary(const CFixedBrane &FixedBrane, const COrbifoldGroup &OrbifoldGroup, const SConfig &VEVConfig, const SUSYMultiplet &Multiplet, bool print_Labels) const
{

	const COrbifoldGroupElement &constr_Element = OrbifoldGroup.GetElement(FixedBrane.Getconstructing_Element());

	const CSpaceGroupElement &Element = constr_Element.SGElement;
	(*this->out) << "  -------------------------------------------------------------------------------------------------------------\n";
	(*this->out) << "  sector:      (m,n,k)  = (" << Element.Get_m() << ", " << Element.Get_n() << ", " << Element.Get_k() << ")\n";
	(*this->out) << "  fixed point:  ";

	if (FixedBrane.GetFixedBraneLabel() != "")
		(*this->out) << FixedBrane.GetFixedBraneLabel() << "\n                ";

	(*this->out) << "n_a   = (";
	this->PrintRational(Element.Get_n(0));
	(*this->out) << ", ";
	this->PrintRational(Element.Get_n(1));
	(*this->out) << ", ";
	this->PrintRational(Element.Get_n(2));
	(*this->out) << ", ";
	this->PrintRational(Element.Get_n(3));
	(*this->out) << ", ";
	this->PrintRational(Element.Get_n(4));
	(*this->out) << ", ";
	this->PrintRational(Element.Get_n(5));
	(*this->out) << ")\n";
	(*this->out) << "  -------------------------------------------------------------------------------------------------------------\n";
	(*this->out) << "  V" << this->underscore << "loc = ";
	this->PrintRational(constr_Element.Shift, constr_Element.Shift.GetLattice());
	(*this->out) << this->endofset << "\n";
	(*this->out) << "  -------------------------------------------------------------------------------------------------------------\n";

	vector<unsigned> FieldIndices;
	FixedBrane.GetFieldIndices(VEVConfig.Fields, Multiplet, FieldIndices);

	if (FieldIndices.size() == 0)
		(*this->out) << "  " << this->cbegin << " empty" << this->cend << "\n";
	else
		this->PrintSpectrum(FieldIndices, VEVConfig, print_Labels);

	(*this->out) << "  -------------------------------------------------------------------------------------------------------------\n";

	return true;
}



/* ########################################################################################
######   PrintSummaryOfSectors(const COrbifold &Orbifold, ...) const                 ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold      : the orbifold of the vev-config "VEVConfig"               ######
######   2) VEVConfig     : contains the fields                                      ######
######   3) Multiplet     : print only those fields of SUSY type "Multiplet"         ######
######   4) print_Labels  : print the field labels?                                  ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
bool CPrint::PrintSummaryOfSectors(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<SUSYMultiplet>  &Multiplet, bool print_Labels) const
{
	// print the gauge group
	this->PrintGaugeGroup(VEVConfig);

	// print the sectors
	const size_t s1 = Orbifold.GetNumberOfSectors();
	for (unsigned i = 0; i < s1; ++i)
	{
       for (int j=0; j<Multiplet.size(); j++)   
	   {
		if (!this->PrintSummary(Orbifold.GetSector(i), VEVConfig, Multiplet[j], print_Labels))
			return false;
		(*this->out) << endl;
	   }
	}

	return true;
}

   

/* ########################################################################################
######   PrintSUSYType(const SUSYMultiplet &Multiplet) const                         ######
######                                                                               ######
######   Version: 04.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Multiplet : print the abbreviation of the SUSY type, for example "l" for ######
######                  "LeftChiral"                                                 ######
######   output:                                                                     ######
######   return value : printed succesfully?                                         ######
######################################################################################## */
bool CPrint::PrintSUSYType(const SUSYMultiplet &Multiplet) const
{
	(*this->out) << this->underscore;
	if (output_type == Tmathematica)
		(*this->out) << ",";

	if (output_type == Tlatex)
		(*this->out) << "{";

	switch (Multiplet)
	{
	case Gauge:							//hacking here!!! non-susy labels
		(*this->out) << "g";
		break;
	case bGauge:
		(*this->out) << "bg";
		break;
	case Scalar:
		(*this->out) << "s";
		break;
	case bScalar:
		(*this->out) << "bs";
		break;
	case LeftFermi:
		(*this->out) << "f";
		break;
	case RightFermi:
		(*this->out) << "bf";
		break;
	case moduli:
		(*this->out) << "m";
		break;
	case bmoduli:
		(*this->out) << "bm";
		break;
	case LeftChiral:
		(*this->out) << "l";
		break;
	case RightChiral:
		(*this->out) << "r";
		break;
	case Vector:
		(*this->out) << "v";
		break;
	case VectorCC:
		(*this->out) << "bv";
		break;
	case Hyper:
		(*this->out) << "h";
		break;
	case Halfhyper:
		(*this->out) << "hh";
		break;
	case Gravity:
		(*this->out) << "grav";
		break;
	case GravityCC:
		(*this->out) << "bgrav";
		break;
	case LCModulus:
		(*this->out) << "m";
		break;
	case RCModulus:
		(*this->out) << "bm";
		break;
	case Tachyon:
		(*this->out) << "tach";
		break;
	case bTachyon:
		(*this->out) << "btach";
		break;
	case Excited_Tachyon:
		(*this->out) << "exc_tach";
		break;
	case Excited_bTachyon:
		(*this->out) << "exc_btach";
		break;
	default:
	{
		cout << "\n  Warning in bool CPrint::PrintSUSYType(...) const: SUSYMultiplet not defined. Return false." << endl;
		return false;
	}
	}
	if (output_type == Tlatex)
		(*this->out) << "}";

	return true;
}



/* ########################################################################################
######   PrintTwist(const CTwistVector &TwistVector) const                           ######
######                                                                               ######
######   Version: 26.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) TwistVector : the CTwistVector object to be printed                      ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintTwist(const CTwistVector &TwistVector) const
{
	this->PrintRational(TwistVector);
}


/* ########################################################################################
######   PrintSummaryOfFixedBranes(const COrbifold &Orbifold, ...) const             ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold      : the orbifold of the vev-config "VEVConfig"               ######
######   2) VEVConfig     : contains the fields                                      ######
######   3) Multiplet     : print only those fields of SUSY type "Multiplet"         ######
######   4) print_Labels  : print the field labels?                                  ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
void CPrint::PrintSummaryOfFixedBranes(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<SUSYMultiplet>  &Multiplet, bool print_Labels) const
{
	// print the gauge group
	this->PrintGaugeGroup(VEVConfig);

	vector<unsigned> FieldIndices;

	unsigned j = 0;
	size_t s2 = 0;

	// run through the sectors
	const size_t s1 = Orbifold.GetNumberOfSectors();
	for (unsigned i = 0; i < s1; ++i)
	{
      for (int k=0; k<Multiplet.size(); k++)   
	  {   
		const CSector &Sector = Orbifold.GetSector(i);
		s2 = Sector.GetNumberOfFixedBranes();
		// run through the fixed points of the current sector
		for (j = 0; j < s2; ++j)
		{
			this->PrintSummary(Sector.GetFixedBrane(j), Orbifold.OrbifoldGroup, VEVConfig, Multiplet[k], print_Labels);
			(*this->out) << endl;
		}

      }
	}
}


/* ########################################################################################
######   PrintU1Charges(const CVector &U1Charges, ...) const                         ######
######                                                                               ######
######   Version: 06.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) U1Charges     : a CVector object containing U(1) charges                 ######
######   2) SymmetryGroup : contains the choice of observable sector                 ######
######   output:                                                                     ######
######   return value     : printed succesfully?                                     ######
######################################################################################## */
bool CPrint::PrintU1Charges(const CVector &U1Charges, const SSymmetryGroup &SymmetryGroup) const
{
	const size_t s1 = SymmetryGroup.observable_sector_U1s.size();
	if (s1 == 0)
		return true;

	(*this->out) << this->vector_open;
	for (unsigned j = 0; j < s1; ++j)
	{
		if (j != 0) (*this->out) << ",";
		this->PrintRational(D2Rat(U1Charges[SymmetryGroup.observable_sector_U1s[j]]), true);
	}
	(*this->out) << this->vector_close;

	return true;
}



/* ########################################################################################
######   PrintU1Directions(const SSymmetryGroup &SymmetryGroup, ...) const           ######
######                                                                               ######
######   Version: 06.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SymmetryGroup : contains the U(1) generators                             ######
######   2) ith_generator : if negative print all U(1)s of "SymmetryGroup", else     ######
######                      only the "ith_generator"                                 ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintU1Directions(const SSymmetryGroup &SymmetryGroup, int ith_generator) const
{
	const vector<vector<double> > &U1Directions = SymmetryGroup.GaugeGroup.u1directions;
	const size_t                  number_of_U1s = U1Directions.size();

	if (number_of_U1s == 0)
	{
		(*this->out) << "  " << this->cbegin << "Number of U(1) generators is zero." << this->cend << "\n";
		return;
	}
	if ((ith_generator >= 0) && (ith_generator >= number_of_U1s))
	{
		(*this->out) << "  " << this->cbegin << "U(1) generator #" << ith_generator+1 << " does not exist. Number of U(1)s: " << number_of_U1s << this->cend << "\n";
		return;
	}

	CVector tmp;

	if (ith_generator < 0)
	{
		unsigned i = 0;

		(*this->out) << "  " << this->cbegin << "U(1) generators:" << this->cend << "\n";
		if (this->output_type == Tmathematica)
			(*this->out) << "  lu1={\n";

		for (i = 0; i < number_of_U1s; ++i)
		{
			tmp = U1Directions[i];
			if (tmp.GetSize() != 16)
			{
				cout << "\n  Warning in void CPrint::PrintU1Directions(...) const: U(1) generator does not have 16 entires. Return." << endl;
				return;
			}

			(*this->out) << "    ";
			this->PrintRational(tmp, E8xE8);

			if (this->output_type == Tstandard)
				(*this->out) << "  " << SymmetryGroup.U1s_AdditionalLabels[i];

			if (this->output_type == Tmathematica)
			{
				if (i != number_of_U1s-1)
					(*this->out) << ",\n";
				else
					(*this->out) << "\n  };\n";
			}
			else
				(*this->out) << "\n";
		}
	}
	else
	{
		(*this->out) << "  " << this->cbegin << "U(1) generator (" << ith_generator+1 << "):" << this->cend << "\n    ";
		if (this->output_type == Tmathematica)
			(*this->out) << "lu1" << ith_generator << "=";

		tmp = U1Directions[ith_generator];
		if (tmp.GetSize() != 16)
		{
			cout << "\n  Warning in void CPrint::PrintU1Directions(...) const: U(1) generator does not have 16 entires. Return." << endl;
			return;
		}

		this->PrintRational(tmp, E8xE8);

		if (this->output_type == Tmathematica)
			(*this->out) << ";";
		(*this->out) << "\n";
	}
}



/* ########################################################################################
######   PrintVector(const CVector &Vector, SelfDualLattice Lattice) const           ######
######                                                                               ######
######   Version: 05.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector  : CVector object to be printed                                   ######
######   2) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintVector(const CVector &Vector, SelfDualLattice Lattice) const
{
	if (Vector.GetSize() == 0)
		return;

	if (Lattice != UNSPECIFIED_LATTICE)
		(*this->out) << this->vector_open;

	unsigned i = 0;
	if ((Vector.GetSize() == 16) && (this->output_type != Tmathematica))
	{
		(*this->out) << setw(5) << setprecision(2) << Vector[0];
		for (i = 1; i < 8; ++i)
			(*this->out) << this->separator << setw(5) << setprecision(2) << Vector[i];

		if (Lattice == E8xE8)
			(*this->out) << this->vector_close << " " << this->vector_open;
		else
			(*this->out) << this->separator;

		(*this->out) << setw(5) << setprecision(2) << Vector[8];
		for (i = 9; i < 16; ++i)
			(*this->out) << this->separator << setw(5) << setprecision(2) << Vector[i];
	}
	else
	{
		(*this->out) << setw(5) << setprecision(2) << Vector[0];
		for (i = 1; i < Vector.GetSize(); ++i)
			(*this->out) << this->separator << setw(5) << setprecision(2) << Vector[i];
	}

	if (Lattice != UNSPECIFIED_LATTICE)
		(*this->out) << this->vector_close;
}



/* ########################################################################################
######   PrintWilsonLines(const CWilsonLines &WilsonLines, bool PrintInfo) const     ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) WilsonLines : CWilsonLines object to be printed                          ######
######   2) PrintInfo   : print some additional information                          ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrint::PrintWilsonLines(const CWilsonLines &WilsonLines, bool PrintInfo) const
{
	if (PrintInfo)
	{
		for(unsigned i = 0; i < 6; ++i)
		{
			const CWilsonLine &Wilson = WilsonLines.GetWilsonLine(i);
			if (!Wilson.GetIs_Zero())
			{
				(*this->out) << "  Wilson line W_" << i+1 << ": Order " << Wilson.GetOrder() << "\n";
				this->PrintRational(Wilson, Wilson.GetLattice());
				(*this->out) << "\n";
			}
		}
		return;
	}

	for(unsigned i = 0; i < 6; ++i)
	{
		(*this->out) << "  ";
		const CWilsonLine &Wilson = WilsonLines.GetWilsonLine(i);
		this->PrintRational(Wilson, Wilson.GetLattice());
		(*this->out) << "\n";
	}
	return;
}



/* ########################################################################################
######   TexSpectrum(const COrbifold &Orbifold, const SConfig &VEVConfig, ...)       ######
######                                                                               ######
######   Version: 15.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold     : the orbifold of the vev-config "VEVConfig"                ######
######   2) VEVConfig    : contains all massless fields                              ######
######   3) FieldIndices : indices of fields to be printed                           ######
######   4) print_Labels : for printing multiple field labels: set of indices, for   ######
######                     each field from "FieldIndices" print these field labels   ######
######   output:                                                                     ######
######   return value    : printed succesfully?                                      ######
######################################################################################## */
bool CPrint::TexSpectrum(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const vector<unsigned> &print_Labels)
{
	const OutputType original_output_type = this->output_type;
	const bool output_type_changed = (original_output_type != Tlatex);
	if (output_type_changed)
		this->SetOutputType(Tlatex);

	const CGaugeGroup &GaugeGroup                        = Orbifold.StandardConfig.SymmetryGroup.GaugeGroup;
	const vector<SDiscreteSymmetry> &DiscreteRSymmetries = Orbifold.OrbifoldGroup.GetSpaceGroup().DiscreteRSymmetries;

	const size_t ds = DiscreteRSymmetries.size();
	const bool PrintRCharges = (ds != 0);

	const size_t number_of_factors = GaugeGroup.factor.size();
	const size_t number_of_U1s     = GaugeGroup.u1directions.size();

	unsigned FieldIndex = 0;
	const vector<CField> &Fields = VEVConfig.Fields;
	const size_t f1 = Fields.size();

	if (f1 == 0)
	{
		cout << "\n  Warning in bool CPrint::TexSpectrum(...): No field specified. Return false." << endl;
		return false;
	}
	const size_t number_of_labels = Fields[0].Labels.size();

	unsigned i = 0;

	const size_t p0 = VEVConfig.SymmetryGroup.observable_sector_GGs.size();
	for (i = 0; i < p0; ++i)
	{
		if (VEVConfig.SymmetryGroup.observable_sector_GGs[i] >= number_of_factors)
		{
			cout << "\n  Warning in bool CPrint::TexSpectrum(...): \"VEVConfig.SymmetryGroup.observable_sector_GGs\" out of range. Return false." << endl;
			return false;
		}
	}

	const size_t p1 = VEVConfig.SymmetryGroup.observable_sector_U1s.size();
	for (i = 0; i < p1; ++i)
	{
		if (VEVConfig.SymmetryGroup.observable_sector_U1s[i] >= number_of_U1s)
		{
			cout << "\n  Warning in bool CPrint::TexSpectrum(...): \"VEVConfig.SymmetryGroup.observable_sector_U1s\" out of range. Return false." << endl;
			return false;
		}
	}

	const size_t p2 = print_Labels.size();
	for (i = 0; i < p2; ++i)
	{
		if (print_Labels[i] >= number_of_labels)
		{
			cout << "\n  Warning in bool CPrint::TexSpectrum(...): \"print_Labels\" out of range. Return false." << endl;
			return false;
		}
	}

	const size_t f2 = FieldIndices.size();
	for (i = 0; i < f2; ++i)
	{
		if (FieldIndices[i] >= f1)
		{
			cout << "\n  Warning in bool CPrint::TexSpectrum(...): Field-index out of range. Return false." << endl;
			return false;
		}
	}

	string tabular_parameter = "|l|l|";
	string empty_line_string = " & ";
	string headline          = "sector & irrep ";

	// begin: define some constants and variables
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	unsigned n = 0;
	unsigned o = 0;

	vector<string>::const_iterator pos;

	const size_t s1 = Orbifold.GetNumberOfSectors();
	size_t s2 = 0;
	size_t s3 = 0;
	size_t s4 = 0;

	bool untwisted_sector_empty  = true;
	bool fixedbrane_empty        = true;
	bool first_line_of_sector    = true;
	bool first_printed_sector    = true;

	const CVector NullVector4(4);
	CVector q_momentum_p(4);
	CVector q_momentum_m(4);	
	CVector q_momentum_pf(4); 
	CVector q_momentum_mf(4); 
	// end: define some constants and variables

    const vector<CSector> &Sectors = Orbifold.GetSectors();  
	std::ostringstream oheadline;

	// begin: create some standard strings
	// U(1) charges
	for (i = 0; i < p1; ++i)
	{
		tabular_parameter += "c";
		empty_line_string += " &";
		if ((VEVConfig.SymmetryGroup.observable_sector_U1s[i] == 0) && (Orbifold.StandardConfig.SymmetryGroup.IsFirstU1Anomalous))
			oheadline << "& $q_\\text{anom}$ ";
		else
			oheadline << "& $q_{" << VEVConfig.SymmetryGroup.observable_sector_U1s[i] << "}$ ";
	}
	if (p1 != 0)
		tabular_parameter += "|";

	// Labels
	for (i = 0; i < p2; ++i)
	{
		tabular_parameter += "c";
		empty_line_string += " &";
		oheadline << "& Label " << (print_Labels[i]+1) << " ";
	}
	if (p2 != 0)
		tabular_parameter += "|";

	headline += oheadline.str();
	headline += "\\\\\n";

	empty_line_string += " \\\\\n";
	// end: create some standard strings

	(*this->out) << "{\\footnotesize\n";
	(*this->out) << "\\begin{longtable}{" << tabular_parameter << "}\n";
	(*this->out) << "\\hline\n" << headline << "\\hline\\hline\n";
	(*this->out) << "\\endhead\n";
	(*this->out) << "\\hline\n";
	(*this->out) << "\\endfoot\n";

	// run through the sectors
	for (i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		first_line_of_sector = true;

		// untwisted sector
		if ((Sector.Get_m() == 0) && (Sector.Get_n() == 0) && (Sector.Get_k() == 0))
		{
			const CFixedBrane &FixedBrane = Sector.GetFixedBrane(0);
			for (j = 1; j < 4; ++j)
			{
				q_momentum_p = NullVector4;
				q_momentum_m = NullVector4;
				q_momentum_p[j] = +1;
				q_momentum_m[j] = -1;
				q_momentum_pf = NullVector4;  
				q_momentum_mf = NullVector4;  		
				if(j==1)                      
				{                             
				q_momentum_pf[j-1] = 2.0;      
				q_momentum_pf[j] = 2.0;        
				q_momentum_pf[j+1] = 2.0;       
				q_momentum_pf[j+2] = 2.0;		
			    }                                

				// first check whether this fixed point contains fields
				// that are elements of FieldIndices
				fixedbrane_empty = true;

				s2 = FixedBrane.GetNumberOfInvariantStates();
				for (k = 0; fixedbrane_empty && (k < s2); ++k)
				{
					const CState &InvariantState = FixedBrane.GetInvariantState(k);

					s3 = InvariantState.GetNumberOfFieldIndices();
					for (l = 0; fixedbrane_empty && (l < s3); ++l)
					{
						FieldIndex = InvariantState.GetFieldIndex(l);
						const CField &Field = Fields[FieldIndex];

						if (find(FieldIndices.begin(), FieldIndices.end(), FieldIndex) != FieldIndices.end())
						{
							const CVector &q_Vector = Field.GetRMWeight(0, Sectors);  
                            q_momentum_mf = q_Vector + q_momentum_pf;   
							if ((q_Vector == q_momentum_p) || (q_Vector == q_momentum_m) || (q_momentum_mf[1]==2.5) || (q_momentum_mf[1]==1.5))
							{
								untwisted_sector_empty = false;
								fixedbrane_empty = false;
							}
						}
					}
				}
				if (!fixedbrane_empty)
				{
					// begin: print (k, l; n_alpha)
					(*this->out) << "$U_" << j << "$ ";
					// end: print (k, l; n_alpha)

					for (k = 0; k < f2; ++k)
					{
						const CField &Field = Fields[FieldIndices[k]];

						if (Field.SGElement.NoTwist())
						{
							
							const CVector &q_Vector = Field.GetRMWeight(0, Sectors);  
                            q_momentum_mf = q_Vector + q_momentum_pf;   

							
							if ((q_Vector == q_momentum_p) || (q_Vector == q_momentum_m) || (q_momentum_mf[1]==2.5) || (q_momentum_mf[1]==1.5))
							{
								const RepVector &Dimensions = Field.Dimensions;

								(*this->out) << " & $";
								this->PrintRep(Dimensions, VEVConfig.SymmetryGroup);
								(*this->out) << "$ ";

								// begin: print U(1) charges
								for (o = 0; o < p1; ++o)
								{
									(*this->out) << "& $";
									this->PrintRational(D2Rat(Field.U1Charges[VEVConfig.SymmetryGroup.observable_sector_U1s[o]]));
									(*this->out) << "$ ";
								}
								// end: print U(1) charges

								// begin: print Labels
								for (o = 0; o < p2; ++o)
								{
									(*this->out) << "& $";
									this->PrintLabel(Field, print_Labels[o]);
									(*this->out) << "$ ";
								}
								// end: print Labels

								(*this->out) << "\\\\" << endl;
								//break;
							}
						}
					}
				}
			}
			if (!untwisted_sector_empty)
			{
				first_printed_sector = false;
				(*this->out) << "\\hline" << endl;
			}
		}

		// twisted sector
		else
		{
			const CTwistVector &localTwist = Sector.GetTwist();

			s2 = Sector.GetNumberOfFixedBranes();

			// run through the fixed branes
			for (j = 0; j < s2; ++j)
			{
				const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

				// first check whether this fixed point contains fields
				// that are elements of RepLabels
				fixedbrane_empty = true;

				s3 = FixedBrane.GetNumberOfInvariantStates();
				if (s3 != 0)
				{
					for (k = 0; fixedbrane_empty && (k < s3); ++k)
					{
						const CState &InvariantState = FixedBrane.GetInvariantState(k);

						s4 = InvariantState.GetNumberOfFieldIndices();
						for (l = 0; fixedbrane_empty && (l < s4); ++l)
						{
							if (find(FieldIndices.begin(), FieldIndices.end(), InvariantState.GetFieldIndex(l)) != FieldIndices.end())
								fixedbrane_empty = false;
						}
					}
				}
				if (!fixedbrane_empty)
				{
					if (first_line_of_sector && !first_printed_sector)
					{
						first_line_of_sector = false;
						(*this->out) << "\\hline\n";
					}

					if (first_printed_sector)
						first_printed_sector = false;

					// begin: print localization (m, n, k; n_alpha)   
					const CSpaceGroupElement &Label = FixedBrane.GetSGElement();
					(*this->out) << "$T_{(" << Label.Get_m() << ", " << Label.Get_n() << ", " << Label.Get_k() << ")}^{(";
		
						 this->PrintRational(Label.Get_n(0),false);
                         (*this->out) << ", ";
                         this->PrintRational(Label.Get_n(1),false);
                         (*this->out) << ", ";
              
						 this->PrintRational(Label.Get_n(2),false);
                         (*this->out) << ", ";
                         this->PrintRational(Label.Get_n(3),false);
                         (*this->out) << ", ";
             
						 this->PrintRational(Label.Get_n(4),false);
                         (*this->out) << ", ";
                         this->PrintRational(Label.Get_n(5),false);
                   
					(*this->out) << ")}$ ";                            
					// end: print localization (k, l; n_alpha)    
                    
					for (k = 0; k < f2; ++k)
					{
						const CField &Field = Fields[FieldIndices[k]];

						if (Field.SGElement == Label)
						{
							const RepVector &Dimensions = Field.Dimensions;

							(*this->out) << " & $";
							this->PrintRep(Dimensions, VEVConfig.SymmetryGroup);
							(*this->out) << "$ ";

							// begin: print U(1) charges
							for (n = 0; n < p1; ++n)
							{
								(*this->out) << "& $";
								this->PrintRational(D2Rat(Field.U1Charges[VEVConfig.SymmetryGroup.observable_sector_U1s[n]]));
								(*this->out) << "$ ";
							}
							// end: print U(1) charges

							// begin: print Labels
							for (n = 0; n < p2; ++n)
							{
								(*this->out) << "& $";
								this->PrintLabel(Field, print_Labels[n]);
								(*this->out) << "$ ";
							}
							// end: print Labels

							(*this->out) << "\\\\" << endl;
						}
					}
				}
			}
		}
	}

	(*this->out) << "\\end{longtable}\n\n";
	(*this->out) << "}" << endl;

	if (output_type_changed)
		this->SetOutputType(original_output_type);

	return true;
}

