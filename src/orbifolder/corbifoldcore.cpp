#include "corbifoldcore.h"
#include "cfixedpoint.h"


/* ########################################################################################
######   COrbifoldCore(const COrbifoldGroup &OrbifoldGroup, ...)                     ######
######                                                                               ######
######   Version: 16.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup : the orbifold-group                                       ######
######   2) Spectrum      : the shift and Wilson line independent part of the model  ######
######                      are saved here (i.e. sectors, fixed branes, oscillators  ######
######                      and the right-movers)                                    ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifoldCore object. Creates the model-independent part,  ######
######   being the (un-)twisted sectors, their fixed branes and for each sector      ######
######   the oscillators and the right-movers.                                       ######
######################################################################################## */
COrbifoldCore::COrbifoldCore(const COrbifoldGroup &OrbifoldGroup, vector<CSector> &Spectrum)
{
  if (OrbifoldGroup.GetOrbifoldGroup_CheckStatus() != CheckedAndGood)
  {
    cout << "\n  Warning in COrbifoldCore::COrbifoldCore(...): Orbifold group is ill-defined. Return." << endl;
    this->OrbifoldCore_CheckStatus = CheckedAndFailed;
    return;
  }

  const vector<COrbifoldGroupElement> &constructing_Elements = OrbifoldGroup.GetElements();
  const size_t s1 = constructing_Elements.size();

  if (s1 == 0)
  {
    cout << "\n  Warning in COrbifoldCore::COrbifoldCore(...): Orbifold-group is empty. Return false." << endl;
    this->OrbifoldCore_CheckStatus = CheckedAndFailed;
    return;
  }
  if (Spectrum.size() != 0)
  {
    cout << "\n  Warning in COrbifoldCore::COrbifoldCore(...): \"Spectrum\" is not empty. Now cleared." << endl;
    Spectrum.clear();
  }

  const vector<CTwistVector> &ZMxZNxZK_Twists = OrbifoldGroup.GetSpaceGroup().GetTwists();

  bool SectorNotFound = true;
  string FixedBraneLabel = "";

  size_t s2 = 0;
  unsigned j = 0;
  unsigned counter = 1;
  unsigned pos_index = 0;
  // run through all constructing elements
  for (unsigned i = 0; i < s1; ++i)
  {
    // create new fixed brane, using the constructing elements
    const COrbifoldGroupElement &OG_Element = constructing_Elements[i]; 
    const CSpaceGroupElement    &SG_Element = OG_Element.SGElement;

    // run through the existing (un-)twisted sectors
    SectorNotFound = true;
    s2 = Spectrum.size();
    for (j = 0; SectorNotFound && (j < s2); ++j)
    {
      if ((Spectrum[j].Get_m() == SG_Element.Get_m()) && (Spectrum[j].Get_n() == SG_Element.Get_n()) && (Spectrum[j].Get_k() == SG_Element.Get_k()))
      {
        pos_index = j;
        SectorNotFound = false;
      }
    }
    // if the sector of the new fixed brane is not known create a new one
    if (SectorNotFound)
    {
      pos_index = Spectrum.size();

      CSector new_Sector;
      Spectrum.push_back(new_Sector);

      CSector &Sector = Spectrum[pos_index];
      Sector.Create(SG_Element.Get_m(), SG_Element.Get_n(), SG_Element.Get_k(), OG_Element.Twist, ZMxZNxZK_Twists);
    }

    // begin: add the new fixed brane to its sector
    if ((SG_Element.Get_m() == 0) && (SG_Element.Get_n() == 0) && (SG_Element.Get_k() == 0))
      FixedBraneLabel = "U";
    else
    {
      std::ostringstream os;
      os << counter;
      ++counter;
      FixedBraneLabel = "T" + os.str();
    }

    Spectrum[pos_index].AddFixedBrane(SG_Element, i, FixedBraneLabel);
    // end: add the new fixed brane to its sector
  }
  this->OrbifoldCore_CheckStatus = CheckedAndGood;      
}



/* ########################################################################################
######   COrbifoldCore()                                                             ######
######                                                                               ######
######   Version: 03.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a COrbifoldCore object.                              ######
######################################################################################## */
COrbifoldCore::~COrbifoldCore()
{
}
