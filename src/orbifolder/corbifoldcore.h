#ifndef CORBIFOLDCORE_H
#define CORBIFOLDCORE_H

#include "corbifoldgroup.h"
#include "csector.h"

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

//! COrbifoldCore.
/*!
COrbifoldCore creates the part of the orbifold that is independent of the choice of shifts and Wilson lines, i.e. the right-moving part of the string.
 */

class COrbifoldCore {
public:
// member functions
  COrbifoldCore(const COrbifoldGroup &OrbifoldGroup, vector<CSector> &Spectrum);
  ~COrbifoldCore();

  const CheckStatus &GetOrbifoldCore_CheckStatus() const {return this->OrbifoldCore_CheckStatus;};

private:
  CheckStatus OrbifoldCore_CheckStatus;
};

#endif
