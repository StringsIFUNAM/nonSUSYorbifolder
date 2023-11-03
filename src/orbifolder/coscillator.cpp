#include "coscillator.h"
#include "globalfunctions.h"

using std::cout;
using std::endl;
using std::exit;



/* ########################################################################################
######   CModedOscillator(MoversType type, unsigned index, bool indexbar, ...)       ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Type               : "LeftMover" or "RightMover"                         ######
######   2) index              : direction of the string excitation                  ######
######                           (in complex coordinates for index = 1,2,3)          ######
######   3) indexbar           : complex conjugate for index = 0,1,2,3               ######
######   4) constructing_Twist : a CTwistVector object                               ######
######   5) m                  : integer that adds up to the oscillator frequency    ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CModedOscillator object. Creates a left- or right-moving   ######
######   oscillator (depending on "Type") pointing into the direction specified by   ######
######   "index" and "indexbar" and with freuqency specified by "constructing_Twist" ######
######   and "m".                                                                    ######
######################################################################################## */
CModedOscillator::CModedOscillator(MoversType type, unsigned index, bool indexbar, const CTwistVector &constructing_Twist, int m)
{
  // Set the precision
  const double prec = 0.00001;

  this->LadderOperator = NOT_DEF_LAD;
  this->index          = index;
  this->indexbar       = indexbar;
  this->type           = type;

  if ((this->index > 25) || (indexbar && (this->index > 3)) || (type == NOT_DEF_MOVTYPE))
  {
    cout << "Warning in bool CModedOscillator::CModedOscillator(...): Index out of range. Return." << endl;
    this->Oscillator_CheckStatus = CheckedAndFailed;
    return;
  }

  double entry = 0.0;
  if (this->index < 4)
  {
    entry = constructing_Twist[this->index];
    if (fabs( ((int)entry) - entry) > prec)
    {
      entry = fmod(entry, 1.0);
      if (entry < 0)
        ++entry;
    }
  }

  if (this->type == LeftMover)
  {
    if (this->indexbar)
      this->frequency      = m + entry;
    else
      this->frequency      = m - entry;
  }
  else
  if (this->type == RightMover)
  {
    if (this->indexbar)
      this->frequency      = m - entry;
    else
      this->frequency      = m + entry;
  }

  this->NumberOperator = -this->frequency;

  if (this->frequency < -prec)
    this->LadderOperator = CreationOperator;

  if (this->frequency > prec)
    this->LadderOperator = AnnihilationOperator;

  this->Oscillator_CheckStatus = CheckedAndGood;
}



/* ########################################################################################
######   ~CModedOscillator()                                                         ######
######                                                                               ######
######   Version: 01.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CModedOscillator object.                           ######
######################################################################################## */
CModedOscillator::~CModedOscillator()
{
}



/* ########################################################################################
######   operator==(const CVector &Oscillator2) const                                ######
######                                                                               ######
######   Version: 14.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Oscillator2 : a CModedOscillator object                                  ######
######   output:                                                                     ######
######   return value   : are the oscillators equal?                                 ######
###########################################################################################
######   description:                                                                ######
######   Compares if the oscillator "Oscillator2" is equal to this oscillator.       ######
######################################################################################## */
bool CModedOscillator::operator==(const CModedOscillator &Oscillator2) const
{
  if (this->LadderOperator != Oscillator2.LadderOperator)
    return false;
  if (this->NumberOperator != Oscillator2.NumberOperator)
    return false;
  if (this->indexbar != Oscillator2.indexbar)
    return false;
  if (this->frequency != Oscillator2.frequency)
    return false;
  if (this->index != Oscillator2.index)
    return false;
  if (this->type != Oscillator2.type)
    return false;

  return true;
}



/* ########################################################################################
######   &GetComplex() const                                                         ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : does the oscillator have a complex conjugated index?         ######
###########################################################################################
######   description:                                                                ######
######   Returns true if the oscillator has a complex index and false otherwise.     ######
######################################################################################## */
const bool &CModedOscillator::GetComplex() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in double CModedOscillator::GetTransformationProperty() const : CModedOscillator ill-defined. Return false." << endl;

  return this->indexbar;
}



/* ########################################################################################
######   &GetFrequency() const                                                       ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the frequency of the oscillator                              ######
###########################################################################################
######   description:                                                                ######
######   Returns the freuqency of the oscillator.                                    ######
######################################################################################## */
const double &CModedOscillator::GetFrequency() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in double CModedOscillator::GetFrequency() const : CModedOscillator ill-defined. Return 0.0." << endl;

  return this->frequency;
}



/* ########################################################################################
######   &GetIndex() const                                                           ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the space-time index of the oscillator                       ######
###########################################################################################
######   description:                                                                ######
######   Returns the index of the compact space (or of the gauge degrees of freedom) ######
######   in which this oscillator points, being 0,1,2,3 or 0,..,25.                  ######
######################################################################################## */
const unsigned &CModedOscillator::GetIndex() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in unsigned CModedOscillator::GetIndex() const : CModedOscillator ill-defined. Return 0." << endl;

  return this->index;
}



/* ########################################################################################
######   GetTransformationProperty(const CTwistVector &Twist) const                  ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Twist     : a CTwistVector object                                        ######
######   output:                                                                     ######
######   return value : transformation phase with respect to the rotation specified  ######
######                  by "Twist"                                                   ######
###########################################################################################
######   description:                                                                ######
######   Returns the eigenvalue of this oscillator under a rotation by "Twist".      ######
######################################################################################## */
double CModedOscillator::GetTransformationProperty(const CTwistVector &Twist) const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in double CModedOscillator::GetTransformationProperty(...) const : CModedOscillator ill-defined. Return 0." << endl;
  
  if (this->index < 4)
  {
    if (this->indexbar)
      return -Twist[this->index];
    else
      return Twist[this->index];
  }
  else
    return 0.0;
}



/* ########################################################################################
######   &GetType() const                                                            ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : either "LeftMover" or "RightMover"                           ######
###########################################################################################
######   description:                                                                ######
######   Returns either "LeftMover" or "RightMover" depending on whether this        ######
######   oscillator acts on the left- or on the right-mover of the string.           ######
######################################################################################## */
const MoversType &CModedOscillator::GetType() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in MoversType CModedOscillator::GetType() const : CModedOscillator ill-defined. Return NOT_DEF_MOVTYPE." << endl;

  return this->type;
}



/* ########################################################################################
######   &GetNumberOperator() const                                                  ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the number operator                                          ######
###########################################################################################
######   description:                                                                ######
######   Returns this oscillators' contribution to the number operator.              ######
######################################################################################## */
const double &CModedOscillator::GetNumberOperator() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in double CModedOscillator::GetNumberOperator() const : CModedOscillator ill-defined. Return 0.0." << endl;

  return this->NumberOperator;
}



/* ########################################################################################
######   &GetLadderOperatorType() const                                              ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : either "CreationOperator" or "AnnihilationOperator"          ######
###########################################################################################
######   description:                                                                ######
######   Returns either "CreationOperator" or "AnnihilationOperator" depending on    ######
######   whether this oscillators' frequency is negative or positive.                ######
######################################################################################## */
const LadderOperatorType &CModedOscillator::GetLadderOperatorType() const
{
  if (this->Oscillator_CheckStatus != CheckedAndGood)
    cout << "\n  Warning in LadderOperatorType CModedOscillator::GetLadderOperatorType() const : CModedOscillator ill-defined. Return NOT_DEF_LAD." << endl;
  
  return this->LadderOperator;
}
  
