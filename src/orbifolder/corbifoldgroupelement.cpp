
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;
using std::exit;

#include "corbifoldgroupelement.h"


/* ########################################################################################
######   COrbifoldGroupElement()                                                     ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a COrbifoldGroupElement object. No content is       ######
######   specified.                                                                  ######
######################################################################################## */
COrbifoldGroupElement::COrbifoldGroupElement()
{
}

/* ########################################################################################
######   COrbifoldGroupElement(SelfDualLattice Lattice)                              ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifoldGroupElement object. The shift is of lattice type ######
######   "Lattice".                                                                  ######
######################################################################################## */
COrbifoldGroupElement::COrbifoldGroupElement(SelfDualLattice Lattice)
  : Shift(Lattice)
{
}


/* ########################################################################################
######   COrbifoldGroupElement(SelfDualLattice Lattice, ...)                         ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement : a CSpaceGroupElement object containing (k,l,n_alpha)         ######
######   2) Shift     : the corresponding shift w V_0 + n V_1 + k V_2 + n_alpha W_alpha      ######
######   3) Twist     : the corresponding twist w v_0 + n v_1 + k v_2                        ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifoldGroupElement object where all content is          ######
######   specified.                                                                  ######
######################################################################################## */
COrbifoldGroupElement::COrbifoldGroupElement(const CSpaceGroupElement &SGElement, const CShiftVector &Shift, const CTwistVector &Twist)
 : SGElement(SGElement), Shift(Shift), Twist(Twist)
{
}



/* ########################################################################################
######   ~COrbifoldGroupElement()                                                    ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a COrbifoldGroupElement object.                      ######
######################################################################################## */
COrbifoldGroupElement::~COrbifoldGroupElement()
{
}
