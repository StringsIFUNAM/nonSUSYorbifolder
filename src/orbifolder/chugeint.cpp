#include "chugeint.h"
#include <iostream>
#include <math.h>

#include <boost/config.hpp>
#include <boost/rational.hpp>
#include "globalfunctions.h"
using boost::rational;


using std::ostream;
using std::flush;
using std::cout;
using std::endl;


using std::vector;

namespace Heap_CHugeInt
{
  CHugeInt ZERO;
  CHugeInt result;
  
  unsigned i = 0;
  unsigned j = 0;
  
  int int_i = 0;
  
  vector<unsigned char> step;
  vector<vector<unsigned char> > steps;
  
  unsigned entry   = 0;
  unsigned carry   = 0;
  unsigned sum     = 0;
  bool          b_carry  = 0;
  char          c_sum    = 0;
  unsigned char uc_carry = 0;
  unsigned char uc_sum   = 0;
  
  CHugeInt tmp;
  vector<CHugeInt> Test_cxb(10, tmp);
  unsigned char Elements_in_Test_cxb(2);
  unsigned c   = 0;
  unsigned pos = 0;
}



/* ##########################################################################
######   operator<<                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
######   checked: 16.04.2009                                           ######
########################################################################## */
ostream& operator << (ostream &os, const CHugeInt &a)
{
  if (!a.sign)
    os << "-";

  for (Heap_CHugeInt::int_i = a.Size-1; Heap_CHugeInt::int_i >= 0; --Heap_CHugeInt::int_i)
    os << (unsigned)a.Digits[Heap_CHugeInt::int_i];

  return os;
}


/* ##########################################################################
######   rational<int> D2Rat(const double x)                           ######
######   Ist die routine falsch? Saul war's!                           ######
######   Version: 18.7.2006                                            ######
########################################################################## */
rational<CHugeInt> D2RatHugeInt(double x)
{
  double tmp=0.0;

  for(int i=1; i<10000; ++i)
  {
    tmp = x * i;
    if(is_integer(tmp))
      return rational<CHugeInt>((CHugeInt)((int)roundf(tmp)),(CHugeInt)i);
  }

  return rational<CHugeInt>((CHugeInt)((int)roundf(10000*x)),(CHugeInt)10000);

/*  double fx = floor(x);
  rational<CHugeInt> result = CHugeInt((int)fx);
  x -= fx;

  const double prec = 0.00001;
  if (fabs(x) < prec)
    return result;

  double tmp = 0.0;
  
  for(int i = 2; i < 20000; ++i)
  {
    tmp = x * i;
    if(fabs(roundf(tmp) - tmp) < prec)
      return result + rational<CHugeInt>(CHugeInt((int)roundf(tmp)), CHugeInt(i));
  }
  cout << "Warning in rational<CHugeInt> D2CHugeInt(const double &x): no conversion from double " << x << " to CHugeInt found. Return 0." << endl;
  return CHugeInt(0);*/
}



double RatHugeInt2D(const rational<CHugeInt> &r)
{
  return ((double)(r.numerator().ToLongLongInt()))/((double)(r.denominator().ToLongLongInt()));
}



rational<int> RatHugeInt2Rat(const rational<CHugeInt> &r)
{
  return rational<int>((int)r.numerator().ToLongLongInt(),(int)r.denominator().ToLongLongInt());
}

/* ##########################################################################
######   operator+                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::operator+(const CHugeInt &b) const
{
  // if a = 0  and  b anything
  if (this->zero)
    return b;

  // if b = 0  and  a != 0
  if (b.zero)
  {
    Heap_CHugeInt::result.sign   = this->sign;
    Heap_CHugeInt::result.zero   = false;  // a != 0
    Heap_CHugeInt::result.Size   = this->Size;
    Heap_CHugeInt::result.Digits = this->Digits;

    return Heap_CHugeInt::result;
  }
  // now  a != 0  and  b != 0

  // a = -b
  if ((this->sign == !b.sign) && (this->Digits == b.Digits))
    return Heap_CHugeInt::ZERO;

  // result != 0
  Heap_CHugeInt::result.zero = false;

  if (this->sign == b.sign)
  {
    Heap_CHugeInt::result.sign   = this->sign;
    Heap_CHugeInt::result.Size   = this->Size;
    Heap_CHugeInt::result.Digits = this->Digits;
    Heap_CHugeInt::result.StandardAddition(b.Digits, b.Size);

    return Heap_CHugeInt::result;
  }

  // a > 0  and  b < 0
  if (this->sign)
  {
    // a < |b|
    if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    {
      // result < 0, e.g. 3 + (-5) = -(5 - 3)
      Heap_CHugeInt::result.sign   = false;
      Heap_CHugeInt::result.Size   = b.Size;
      Heap_CHugeInt::result.Digits = b.Digits;
      Heap_CHugeInt::result.StandardSubtraction(this->Digits, this->Size);

      return Heap_CHugeInt::result;
    }
    // a > |b|  because  a != -b
    else
    {
      // result > 0, e.g. 5 + (-3) = 5 - 3
      Heap_CHugeInt::result.sign   = true;
      Heap_CHugeInt::result.Size   = this->Size;
      Heap_CHugeInt::result.Digits = this->Digits;
      Heap_CHugeInt::result.StandardSubtraction(b.Digits, b.Size);

      return Heap_CHugeInt::result;
    }
  }
  // a < 0  and  b > 0
  else
  {
    // |a| < b
    if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    {
      // result > 0, e.g. -3 + 5 = 5 - 3
      Heap_CHugeInt::result.sign   = true;
      Heap_CHugeInt::result.Size   = b.Size;
      Heap_CHugeInt::result.Digits = b.Digits;
      Heap_CHugeInt::result.StandardSubtraction(this->Digits, this->Size);

      return Heap_CHugeInt::result;
    }
    // |a| > b  because  a != -b
    else
    {
      // result < 0, e.g. -5 + 3 = -(5 - 3)
      Heap_CHugeInt::result.sign   = false;
      Heap_CHugeInt::result.Size   = this->Size;
      Heap_CHugeInt::result.Digits = this->Digits;
      Heap_CHugeInt::result.StandardSubtraction(b.Digits, b.Size);

      return Heap_CHugeInt::result;
    }
  }
}




/* ##########################################################################
######   operator-                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::operator-(const CHugeInt &b) const
{
  // if b = 0  and  a anything
  if (b.zero)
  {
    Heap_CHugeInt::result.sign   = this->sign;
    Heap_CHugeInt::result.zero   = this->zero;
    Heap_CHugeInt::result.Size   = this->Size;
    Heap_CHugeInt::result.Digits = this->Digits;

    return Heap_CHugeInt::result;
  }

  // if a = 0  and  b != 0
  if (this->zero)
  {
    Heap_CHugeInt::result.sign   = !b.sign;
    Heap_CHugeInt::result.zero   = false;
    Heap_CHugeInt::result.Size   = b.Size;
    Heap_CHugeInt::result.Digits = b.Digits;

    return Heap_CHugeInt::result;
  }
  // now  a != 0  and  b != 0

  // a = b
  if ((this->sign == b.sign) && (this->Digits == b.Digits))
    return Heap_CHugeInt::ZERO;

  // result != 0
  Heap_CHugeInt::result.zero = false;
  
  // a > 0
  if (this->sign)
  {
    // b > 0
    if (b.sign)
    {
      // a < b
      if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
      {
        // result < 0, e.g. 3 - 5 = -(5 - 3)
        Heap_CHugeInt::result.sign   = false;
        Heap_CHugeInt::result.Size   = b.Size;
        Heap_CHugeInt::result.Digits = b.Digits;
        Heap_CHugeInt::result.StandardSubtraction(this->Digits, this->Size);

        return Heap_CHugeInt::result;
      }
      // a > b  and  a != b
      else
      {
        // result > 0, e.g. 5 - 3
        Heap_CHugeInt::result.sign   = true;
        Heap_CHugeInt::result.Size   = this->Size;
        Heap_CHugeInt::result.Digits = this->Digits;
        Heap_CHugeInt::result.StandardSubtraction(b.Digits, b.Size);

        return Heap_CHugeInt::result;
      }
    }
    // b < 0
    else
    {
      // result > 0, e.g. 5 - (-3) = 5 + 3
      Heap_CHugeInt::result.sign   = true;
      Heap_CHugeInt::result.Size   = this->Size;
      Heap_CHugeInt::result.Digits = this->Digits;
      Heap_CHugeInt::result.StandardAddition(b.Digits, b.Size);

      return Heap_CHugeInt::result;
    }
  }
  // a < 0
  else
  {
    // b > 0
    if (b.sign)
    {
      // result < 0, e.g. (-5) - 3 = -(5 + 3)
      Heap_CHugeInt::result.sign   = false;
      Heap_CHugeInt::result.Size   = this->Size;
      Heap_CHugeInt::result.Digits = this->Digits;
      Heap_CHugeInt::result.StandardAddition(b.Digits, b.Size);

      return Heap_CHugeInt::result;
    }
    // b < 0
    else
    {
      // |a| < |b|
      if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
      {
        // result > 0, e.g. (-3) - (-5) = 5 - 3
        Heap_CHugeInt::result.sign   = true;
        Heap_CHugeInt::result.Size   = b.Size;
        Heap_CHugeInt::result.Digits = b.Digits;
        Heap_CHugeInt::result.StandardSubtraction(this->Digits, this->Size);

        return Heap_CHugeInt::result;
      }
      // |a| > |b|  and  a != b
      else
      {
        // result < 0, e.g. (-5) - (-3) = -(5 - 3)
        Heap_CHugeInt::result.sign   = false;
        Heap_CHugeInt::result.Size   = this->Size;
        Heap_CHugeInt::result.Digits = this->Digits;
        Heap_CHugeInt::result.StandardSubtraction(b.Digits, b.Size);

        return Heap_CHugeInt::result;
      }
    }
  }
}



/* ##########################################################################
######   operator*                                                     ######
######                                                                 ######
######   Version: 15.09.2011                                           ######
########################################################################## */
CHugeInt CHugeInt::operator*(const CHugeInt &b) const
{
  // if a = 0 or b = 0
  if (this->zero || b.zero)
    return Heap_CHugeInt::ZERO;

  // now  a != 0  and  b != 0
  Heap_CHugeInt::result.sign = (this->sign == b.sign);
  Heap_CHugeInt::result.zero = false;

  if ((b.Size == 1) && (b.Digits[0] == 1))
  {
    Heap_CHugeInt::result.Size   = this->Size;
    Heap_CHugeInt::result.Digits = this->Digits;

    return Heap_CHugeInt::result;
  }

  if ((this->Size == 1) && (this->Digits[0] == 1))
  {
    Heap_CHugeInt::result.Size   = b.Size;
    Heap_CHugeInt::result.Digits = b.Digits;

    return Heap_CHugeInt::result;
  }

  Heap_CHugeInt::steps.clear();
  Heap_CHugeInt::steps.reserve(b.Size);
  
  Heap_CHugeInt::entry = 0;
  Heap_CHugeInt::carry = 0;
  Heap_CHugeInt::sum   = 0;
  
  for (Heap_CHugeInt::i = 0; Heap_CHugeInt::i < b.Size; ++Heap_CHugeInt::i)
  {
    Heap_CHugeInt::entry = b.Digits[Heap_CHugeInt::i];
    if (Heap_CHugeInt::entry != 0)  
    {
      Heap_CHugeInt::step.assign(Heap_CHugeInt::i + this->Size,0);
      
      for (Heap_CHugeInt::j = 0; Heap_CHugeInt::j < this->Size; ++Heap_CHugeInt::j)
      {
        Heap_CHugeInt::sum = Heap_CHugeInt::carry + (Heap_CHugeInt::entry * this->Digits[Heap_CHugeInt::j]);
    
        if (Heap_CHugeInt::sum > 9)
        {
          Heap_CHugeInt::step[Heap_CHugeInt::i + Heap_CHugeInt::j] = Heap_CHugeInt::sum % 10;
          Heap_CHugeInt::carry = Heap_CHugeInt::sum/10;
        }
        else
        {
          Heap_CHugeInt::step[Heap_CHugeInt::i + Heap_CHugeInt::j] = Heap_CHugeInt::sum;
          Heap_CHugeInt::carry = 0;
        }
      }
      while(Heap_CHugeInt::carry != 0)
      {
        Heap_CHugeInt::step.push_back(Heap_CHugeInt::carry % 10);
        Heap_CHugeInt::carry /= 10;
      }
      Heap_CHugeInt::steps.push_back(Heap_CHugeInt::step);
    }
  }
  Heap_CHugeInt::result.Digits.reserve((this->Size + 1) * (b.Size + 1));
  Heap_CHugeInt::result.AddAllPositive(Heap_CHugeInt::steps);
  
  return Heap_CHugeInt::result;
}



/* ##########################################################################
######   operator/                                                     ######
######                                                                 ######
######   Version: 15.09.2011                                           ######
########################################################################## */
CHugeInt CHugeInt::operator/(const CHugeInt &b) const
{
  if (b.zero)
  {
    cout << "Error in CHugeInt CHugeInt::operator/(const CHugeInt &b) const. Divide by zero: b = " << b << ". Return 0." << endl;
    return Heap_CHugeInt::ZERO;
  }

  // if a = 0
  if (this->zero)
    return Heap_CHugeInt::ZERO;

  // now a != 0

  // if b = \pm 1 then return \pm a
  if ((b.Size == 1) && (b.Digits[0] == 1))
  {
    if (b.sign)
      return *this;
    else
    {
      Heap_CHugeInt::result.sign   = !this->sign;
      Heap_CHugeInt::result.Size   = this->Size;
      Heap_CHugeInt::result.Digits = this->Digits;
      Heap_CHugeInt::result.zero   = false;

      return Heap_CHugeInt::result;
    }
  }

  // if |a| < |b| then return a/b = 0
  if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    return Heap_CHugeInt::ZERO;

  Heap_CHugeInt::result.sign = (this->sign == b.sign);
  Heap_CHugeInt::result.zero = false;

  // if |a| = |b| then return a/b = \pm 1
  if (this->Digits == b.Digits)
  {
    Heap_CHugeInt::result.Digits.assign(1,1);
    Heap_CHugeInt::result.Size = 1;

    return Heap_CHugeInt::result;
  }

  Heap_CHugeInt::result.Digits.clear();

  Heap_CHugeInt::tmp.Digits.clear();
  Heap_CHugeInt::tmp.sign = true;
  Heap_CHugeInt::tmp.zero = false;
  
  Heap_CHugeInt::Test_cxb[1] = b;
  Heap_CHugeInt::Elements_in_Test_cxb = 2;
  
  Heap_CHugeInt::pos = this->Size;
  while (Heap_CHugeInt::pos > 0)
  {
    --Heap_CHugeInt::pos;
    Heap_CHugeInt::tmp.Digits.insert(Heap_CHugeInt::tmp.Digits.begin(), this->Digits[Heap_CHugeInt::pos]);
    
    while ((Heap_CHugeInt::tmp.Digits.size() > 1) && (Heap_CHugeInt::tmp.Digits[Heap_CHugeInt::tmp.Digits.size()-1] == 0))
      Heap_CHugeInt::tmp.Digits.pop_back();
    
    Heap_CHugeInt::tmp.Size = Heap_CHugeInt::tmp.Digits.size();
    
    Heap_CHugeInt::c = 0;
    while ((Heap_CHugeInt::c < 10) && this->StandardLess(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::tmp.Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size, Heap_CHugeInt::tmp.Size, true))
    {
      ++Heap_CHugeInt::c;
      if ((Heap_CHugeInt::c == Heap_CHugeInt::Elements_in_Test_cxb) && (Heap_CHugeInt::c < 10))
      {
        ++Heap_CHugeInt::Elements_in_Test_cxb;
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c] = Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c-1];
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].StandardAddition(b.Digits, b.Size);
      }
    }
    Heap_CHugeInt::result.Digits.insert(Heap_CHugeInt::result.Digits.begin(), --Heap_CHugeInt::c);
    
    if (!Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].zero)
      Heap_CHugeInt::tmp.StandardSubtraction(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size);
  }
  
  while ((Heap_CHugeInt::result.Digits.size() > 1) && (Heap_CHugeInt::result.Digits[Heap_CHugeInt::result.Digits.size()-1] == 0))
    Heap_CHugeInt::result.Digits.pop_back();
  
  Heap_CHugeInt::result.Size = Heap_CHugeInt::result.Digits.size();

  return Heap_CHugeInt::result;
}


/* ##########################################################################
######   operator%                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::operator%(const CHugeInt &b) const
{
  if (b.zero)
  {
    cout << "Error in CHugeInt CHugeInt::operator%(const CHugeInt &b) const. Divide by zero: b = " << b << ". Return 0." << endl;
    return Heap_CHugeInt::ZERO;
  }

  // if b = \pm 1  or  |a| = |b|
  if (((b.Size == 1) && (b.Digits[0] == 1)) || (this->Digits == b.Digits))
    return Heap_CHugeInt::ZERO;

  // if a = 0  or  |a| < |b| then return a
  if (this->zero || this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    return *this;

  Heap_CHugeInt::tmp.Digits.clear();
  Heap_CHugeInt::tmp.sign = this->sign;
  
  Heap_CHugeInt::Test_cxb[1] = b;
  Heap_CHugeInt::Elements_in_Test_cxb = 2;
  
  Heap_CHugeInt::pos = this->Size;
  while (Heap_CHugeInt::pos > 0)
  {
    --Heap_CHugeInt::pos;
    Heap_CHugeInt::tmp.Digits.insert(Heap_CHugeInt::tmp.Digits.begin(), this->Digits[Heap_CHugeInt::pos]);
    
    while ((Heap_CHugeInt::tmp.Digits.size() > 1) && (Heap_CHugeInt::tmp.Digits[Heap_CHugeInt::tmp.Digits.size()-1] == 0))
      Heap_CHugeInt::tmp.Digits.pop_back();
    
    Heap_CHugeInt::tmp.Size = Heap_CHugeInt::tmp.Digits.size();
    
    Heap_CHugeInt::c = 0;
    while ((Heap_CHugeInt::c < 10) && this->StandardLess(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::tmp.Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size, Heap_CHugeInt::tmp.Size, true))
    {
      ++Heap_CHugeInt::c;
      if ((Heap_CHugeInt::c == Heap_CHugeInt::Elements_in_Test_cxb) && (Heap_CHugeInt::c < 10))
      {
        ++Heap_CHugeInt::Elements_in_Test_cxb;
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c] = Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c-1];
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].StandardAddition(b.Digits, b.Size);
      }
    }
    if (!Heap_CHugeInt::Test_cxb[--Heap_CHugeInt::c].zero)
      Heap_CHugeInt::tmp.StandardSubtraction(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size);
  }
  Heap_CHugeInt::tmp.Size = Heap_CHugeInt::tmp.Digits.size();
  if ((Heap_CHugeInt::tmp.Size == 1) && (Heap_CHugeInt::tmp.Digits[0] == 0))
    Heap_CHugeInt::tmp.zero = true;
  else
    Heap_CHugeInt::tmp.zero = false;
  
  return Heap_CHugeInt::tmp;
}



/* ##########################################################################
######   operator+=                                                    ######
######                                                                 ######
######   Version: 24.01.2007                                           ######
########################################################################## */
void CHugeInt::operator+=(const CHugeInt &b)
{
  // if b = 0  and  a anything
  if (b.zero)
    return;

  // if a = 0  and  b != 0
  if (this->zero)
  {
    this->sign   = b.sign;
    this->zero   = false; // b != 0
    this->Size   = b.Size;
    this->Digits = b.Digits;

    return;
  }
  // now  a != 0  and  b != 0

  // a = -b
  if ((this->sign == !b.sign) && (this->Digits == b.Digits))
  {
    this->sign = true;
    this->zero = true;
    this->Size = 1;
    this->Digits.assign(1,0);

    return;
  }

  // result != 0
  this->zero = false;
  
  if (this->sign == b.sign)
  {
    this->StandardAddition(b.Digits, b.Size);

    return;
  }
  
  // a > 0  and  b < 0
  if (this->sign)
  {
    // a < |b|
    if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    {
      // result < 0, e.g. 3 + (-5) = -(5 - 3)
      vector<unsigned char> DigitsA = this->Digits;
      this->StandardSubtraction(b.Digits, DigitsA);
      this->sign = false;

      return;
    }
    // a > |b|  because  a != -b
    else
    {
      // result > 0, e.g. 5 + (-3) = 5 - 3
      this->StandardSubtraction(b.Digits, b.Size);
      //this->sign is true already
      
      return;
    }
  }
  // a < 0  and  b > 0
  else
  {
    // |a| < b
    if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    {
      // result > 0, e.g. -3 + 5 = 5 - 3
      vector<unsigned char> DigitsA = this->Digits;
      this->StandardSubtraction(b.Digits, DigitsA);
      this->sign = true;

      return;
    }
    // |a| > b  because  a != -b
    else
    {
      // result < 0, e.g. -5 + 3 = -(5 - 3)
      this->StandardSubtraction(b.Digits, b.Size);
      //this->sign is false already

      return;
    }
  }
}



/* ##########################################################################
######   operator-=                                                    ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
########################################################################## */
void CHugeInt::operator-=(const CHugeInt &b)
{
  // if b = 0  and  a anything
  if (b.zero)
    return;

  // if a = 0  and  b != 0
  if (this->zero)
  {
    this->sign   = !b.sign;
    this->zero   = false;
    this->Size   = b.Size;
    this->Digits = b.Digits;

    return;
  }
  // now  a != 0  and  b != 0

  // a = b
  if ((this->sign == b.sign) && (this->Digits == b.Digits))
  {
    this->sign = true;
    this->zero = true;
    this->Size = 1;
    this->Digits.assign(1,0);

    return;
  }

  // result != 0
  this->zero = false;
  
  // a > 0
  if (this->sign)
  {
    // b > 0
    if (b.sign)
    {
      // a < b
      if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
      {
        // result < 0, e.g. 3 - 5 = -(5 - 3)
        vector<unsigned char> DigitsA = this->Digits;
        this->StandardSubtraction(b.Digits, DigitsA);
        this->sign = false;

        return;
      }
      // a > b  and  a != b
      else
      {
        // result > 0, e.g. 5 - 3
        this->StandardSubtraction(b.Digits, b.Size);
        //this->sign is true already

        return;
      }
    }
    // b < 0
    else
    {
      // result > 0, e.g. 5 - (-3) = 5 + 3
      this->StandardAddition(b.Digits, b.Size);
      //this->sign is true already

      return;
    }
  }
  // a < 0
  else
  {
    // b > 0
    if (b.sign)
    {
      // result < 0, e.g. (-5) - 3 = -(5 + 3)
      this->StandardAddition(b.Digits, b.Size);
      //this->sign is false already

      return;
    }
    // b < 0
    else
    {
      // |a| < |b|
      if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
      {
        // result > 0, e.g. (-3) - (-5) = 5 - 3
        vector<unsigned char> DigitsA = this->Digits;
        this->StandardSubtraction(b.Digits, DigitsA);
        this->sign = true;

        return;
      }
      // |a| > |b|  and  a != b
      else
      {
        // result < 0, e.g. (-5) - (-3) = -(5 - 3)
        this->StandardSubtraction(b.Digits, b.Size);
        //this->sign is false already

        return;
      }
    }
  }
}



/* ##########################################################################
######   operator*=                                                    ######
######                                                                 ######
######   Version: 15.09.2011                                           ######
########################################################################## */
void CHugeInt::operator*=(const CHugeInt &b)
{
  // if a = 0
  if (this->zero)
    return;

  // if b = 0
  if (b.zero)
  {
    this->sign = true;
    this->zero = true;
    this->Size = 1;
    this->Digits.assign(1,0);

    return;
  }
  // now  a != 0  and  b != 0

  this->sign = (this->sign == b.sign);
  this->zero = false;

  // if b = \pm 1
  if ((b.Size == 1) && (b.Digits[0] == 1))
    return;

  // if a = \pm 1
  if ((this->Size == 1) && (this->Digits[0] == 1))
  {
    this->Size   = b.Size;
    this->Digits = b.Digits;

    return;
  }

  Heap_CHugeInt::steps.clear();
  Heap_CHugeInt::steps.reserve(b.Size);
  
  Heap_CHugeInt::entry = 0;
  Heap_CHugeInt::carry = 0;
  Heap_CHugeInt::sum  = 0;
  
  for (Heap_CHugeInt::i = 0; Heap_CHugeInt::i < b.Size; ++Heap_CHugeInt::i)
  {
    Heap_CHugeInt::entry = b.Digits[Heap_CHugeInt::i];
    if (Heap_CHugeInt::entry != 0)  
    {
      Heap_CHugeInt::step.assign(Heap_CHugeInt::i + this->Size,0);
      
      for (Heap_CHugeInt::j = 0; Heap_CHugeInt::j < this->Size; ++Heap_CHugeInt::j)
      {
        Heap_CHugeInt::sum = Heap_CHugeInt::carry + (Heap_CHugeInt::entry * this->Digits[Heap_CHugeInt::j]);
    
        if (Heap_CHugeInt::sum > 9)
        {
          Heap_CHugeInt::step[Heap_CHugeInt::i+Heap_CHugeInt::j] = Heap_CHugeInt::sum % 10;
          Heap_CHugeInt::carry = Heap_CHugeInt::sum/10;
        }
        else
        {
          Heap_CHugeInt::step[Heap_CHugeInt::i+Heap_CHugeInt::j] = Heap_CHugeInt::sum;
          Heap_CHugeInt::carry = 0;
        }
      }
      while(Heap_CHugeInt::carry != 0)
      {
        Heap_CHugeInt::step.push_back(Heap_CHugeInt::carry % 10);
        Heap_CHugeInt::carry /= 10;
      }
      Heap_CHugeInt::steps.push_back(Heap_CHugeInt::step);
    }
  }
  this->Digits.reserve((this->Size + 1) * (b.Size + 1));
  this->AddAllPositive(Heap_CHugeInt::steps);
}



/* ##########################################################################
######   operator/=                                                    ######
######                                                                 ######
######   Version: 15.09.2011                                           ######
########################################################################## */
void CHugeInt::operator/=(const CHugeInt &b)
{
  if (b.zero)
  {
    cout << "Error in void CHugeInt::operator/=(const CHugeInt &b). Divide by zero: b = " << b << ". Return." << endl;
    return;
  }

  // if a = 0
  if (this->zero)
  {
    this->Digits.assign(1,0);
    this->sign = true;
    this->zero = true;
    this->Size = 1;

    return;
  }

  // now a != 0

  // if b = \pm 1 then return \pm a
  if ((b.Size == 1) && (b.Digits[0] == 1))
  {
    if (b.sign)
      return;
    else
    {
      this->sign = !this->sign;
      return;
    }
  }

  // if |a| < |b| then return a/b = 0
  if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
  {
    this->Digits.assign(1,0);
    this->sign = true;
    this->zero = true;
    this->Size = 1;

    return;
  }

  this->sign = (this->sign == b.sign);
  this->zero = false;

  // if |a| = |b| then return a/b = \pm 1
  if (this->Digits == b.Digits)
  {
    this->Digits.assign(1,1);
    this->Size = 1;

    return;
  }
  
  Heap_CHugeInt::tmp.Digits.clear();
  Heap_CHugeInt::tmp.sign = true;
  Heap_CHugeInt::tmp.zero = false;
  
  Heap_CHugeInt::Test_cxb[1] = b;
  Heap_CHugeInt::Elements_in_Test_cxb = 2;
  
  Heap_CHugeInt::result = Heap_CHugeInt::ZERO;
  Heap_CHugeInt::result.Digits.clear();
  
  Heap_CHugeInt::pos = this->Size;
  while (Heap_CHugeInt::pos > 0)
  {
    --Heap_CHugeInt::pos;
    Heap_CHugeInt::tmp.Digits.insert(Heap_CHugeInt::tmp.Digits.begin(), this->Digits[Heap_CHugeInt::pos]);
    
    while ((Heap_CHugeInt::tmp.Digits.size() > 1) && (Heap_CHugeInt::tmp.Digits[Heap_CHugeInt::tmp.Digits.size()-1] == 0))
      Heap_CHugeInt::tmp.Digits.pop_back();
    
    Heap_CHugeInt::tmp.Size = Heap_CHugeInt::tmp.Digits.size();
    
    Heap_CHugeInt::c = 0;
    while ((Heap_CHugeInt::c < 10) && this->StandardLess(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::tmp.Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size, Heap_CHugeInt::tmp.Size, true))
    {
      ++Heap_CHugeInt::c;
      if ((Heap_CHugeInt::c == Heap_CHugeInt::Elements_in_Test_cxb) && (Heap_CHugeInt::c < 10))
      {
        ++Heap_CHugeInt::Elements_in_Test_cxb;
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c] = Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c-1];
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].StandardAddition(b.Digits, b.Size);
      }
    }
    Heap_CHugeInt::result.Digits.insert(Heap_CHugeInt::result.Digits.begin(), --Heap_CHugeInt::c);
    
    if (!Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].zero)
      Heap_CHugeInt::tmp.StandardSubtraction(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size);
  }
    
  while ((Heap_CHugeInt::result.Digits.size() > 1) && (Heap_CHugeInt::result.Digits[Heap_CHugeInt::result.Digits.size()-1] == 0))
    Heap_CHugeInt::result.Digits.pop_back();
  
  this->Digits = Heap_CHugeInt::result.Digits;
  this->Size   = this->Digits.size();
}



/* ##########################################################################
######   operator%=                                                    ######
######                                                                 ######
######   Version: 15.09.2011                                           ######
########################################################################## */
void CHugeInt::operator%=(const CHugeInt &b)
{
  if (b.zero)
  {
    cout << "Error in void CHugeInt::operator%=(const CHugeInt &b). Divide by zero: b = " << b << ". Return." << endl;
    return;
  }

  // if a = 0
  if (this->zero)
    return;

  // now a != 0

  Heap_CHugeInt::b_carry = (b.Size == 1);
  // if b = 0  or  |a| = |b|
  if ((Heap_CHugeInt::b_carry && (b.Digits[0] == 1)) || (this->Digits == b.Digits))
  {
    this->Digits.assign(1,0);
    this->sign = true;
    this->zero = true;
    this->Size = 1;

    return;
  }

  // if b = 2
  if (Heap_CHugeInt::b_carry && (b.Digits[0] == 2))
  {
    // if a is even
    if ((this->Digits[0] % 2) == 0)
    {
      this->Digits.assign(1,0);
      this->sign = true;
      this->zero = true;
      this->Size = 1;

      return;
    }
    // if a is odd
    else
    {
      // keep the sign of a so that, e.g. for a = -27 we have -27 % 2 = -1
      this->Digits.assign(1,1);
      this->zero = false;
      this->Size = 1;

      return;
    }
  }

  // if |a| < |b| then return a
  if (this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false))
    return;

  Heap_CHugeInt::tmp.Digits.clear();
  
  Heap_CHugeInt::Test_cxb[1] = b;
  Heap_CHugeInt::Elements_in_Test_cxb = 2;
  
  Heap_CHugeInt::pos = this->Size;
  while (Heap_CHugeInt::pos > 0)
  {
    --Heap_CHugeInt::pos;
    Heap_CHugeInt::tmp.Digits.insert(Heap_CHugeInt::tmp.Digits.begin(), this->Digits[Heap_CHugeInt::pos]);
    
    while ((Heap_CHugeInt::tmp.Digits.size() > 1) && (Heap_CHugeInt::tmp.Digits[Heap_CHugeInt::tmp.Digits.size()-1] == 0))
      Heap_CHugeInt::tmp.Digits.pop_back();
    
    Heap_CHugeInt::tmp.Size = Heap_CHugeInt::tmp.Digits.size();
    
    Heap_CHugeInt::c = 0;
    while ((Heap_CHugeInt::c < 10) && this->StandardLess(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::tmp.Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size, Heap_CHugeInt::tmp.Size, true))
    {
      ++Heap_CHugeInt::c;
      if ((Heap_CHugeInt::c == Heap_CHugeInt::Elements_in_Test_cxb) && (Heap_CHugeInt::c < 10))
      {
        ++Heap_CHugeInt::Elements_in_Test_cxb;
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c] = Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c-1];
        Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].StandardAddition(b.Digits, b.Size);
      }
    }
    if (!Heap_CHugeInt::Test_cxb[--Heap_CHugeInt::c].zero)
      Heap_CHugeInt::tmp.StandardSubtraction(Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Digits, Heap_CHugeInt::Test_cxb[Heap_CHugeInt::c].Size);
  }
  
  if ((Heap_CHugeInt::tmp.Size == 1) && (Heap_CHugeInt::tmp.Digits[0] == 0))
  {
    this->sign = true;
    this->zero = true;
  }
  else
    this->zero = false;
    
  this->Digits = Heap_CHugeInt::tmp.Digits;
  this->Size   = this->Digits.size();
}



void CHugeInt::operator++()
{
  CHugeInt One(1);
  this->operator+=(One);
}



void CHugeInt::operator--()
{
  CHugeInt One(1);
  this->operator-=(One);
}



/* ##########################################################################
######   operator+                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::operator+() const
{
  Heap_CHugeInt::result.Digits = this->Digits;
  Heap_CHugeInt::result.sign   = this->sign;
  Heap_CHugeInt::result.zero   = this->zero;
  Heap_CHugeInt::result.Size   = this->Size;

  return Heap_CHugeInt::result;
}


/* ##########################################################################
######   operator-                                                     ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::operator-() const
{
  Heap_CHugeInt::result.Digits = this->Digits;
  Heap_CHugeInt::result.zero   = this->zero;
  Heap_CHugeInt::result.Size   = this->Size;
  if (this->zero)
    Heap_CHugeInt::result.sign = true;
  else
    Heap_CHugeInt::result.sign = !this->sign;

  return Heap_CHugeInt::result;
}



/* ##########################################################################
######   operator!                                                     ######
######                                                                 ######
######   Version: 11.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
bool CHugeInt::operator!() const
{
  if (this->zero)
    return true;
  else
    return false;
}



/* ##########################################################################
######   operator=                                                     ######
######                                                                 ######
######   Version: 11.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
void CHugeInt::operator=(const CHugeInt &b)
{
  this->Digits = b.Digits;
  this->sign   = b.sign;
  this->zero   = b.zero;
  this->Size   = b.Size;
}



/* ##########################################################################
######   operator=                                                     ######
######                                                                 ######
######   Version: 06.02.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
void CHugeInt::operator=(const int &a)
{
  if (a == 0)
  {
    this->sign = true;
    this->zero = true;
    this->Size = 1;
    this->Digits.assign(1,0);

    return;
  }

  // now a != 0

  this->Digits.clear();
  this->sign = (a > 0);
  this->zero = false;

  unsigned div = (unsigned)a;

  this->Digits.push_back(div % 10);
  div /= 10;
  while(div != 0)
  {
    this->Digits.push_back(div % 10);
    div /= 10;
  }
  this->Size = this->Digits.size();
}


/* ##########################################################################
######   operator=                                                     ######
######                                                                 ######
######   Version: 06.02.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
void CHugeInt::operator=(const unsigned &a)
{
  this->sign = true;

  if (a == 0)
  {
    this->zero = true;
    this->Size = 1;
    this->Digits.assign(1,0);

    return;
  }

  // now a != 0

  this->Digits.clear();
  this->zero = false;
  unsigned div = a;

  this->Digits.push_back(div % 10);
  div /= 10;
  while(div != 0)
  {
    this->Digits.push_back(div % 10);
    div /= 10;
  }
  this->Size = this->Digits.size();
}


    
/* ##########################################################################
######   operator==                                                    ######
######                                                                 ######
######   Version: 10.01.2008                                           ######
########################################################################## */
bool CHugeInt::operator==(const CHugeInt &b) const
{
  if ((this->sign != b.sign) || (this->Size != b.Size) || (this->Digits != b.Digits))
    return false;

  return true;
}



/* ##########################################################################
######   operator!=                                                    ######
######                                                                 ######
######   Version: 10.01.2008                                           ######
########################################################################## */
bool CHugeInt::operator!=(const CHugeInt &b) const
{
  if ((this->sign != b.sign) || (this->Size != b.Size) || (this->Digits != b.Digits))
    return true;

  return false;
}



/* ##########################################################################
######   operator>=                                                    ######
######                                                                 ######
######   Version: 19.12.2007                                           ######
########################################################################## */
bool CHugeInt::operator>=(const CHugeInt &b) const
{
  if (this->sign == b.sign)
  {
    if (this->sign)
      return this->StandardLess(b.Digits, this->Digits, b.Size, this->Size, true);
    else
      return this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, true);
  }
  else
  {
    if (this->sign)
      return true;
    else
      return false;
  }
}



/* ##########################################################################
######   operator<=                                                    ######
######                                                                 ######
######   Version: 19.12.2007                                           ######
########################################################################## */
bool CHugeInt::operator<=(const CHugeInt &b) const
{
  if (this->sign == b.sign)
  {
    if (this->sign)
      return this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, true);
    else
      return this->StandardLess(b.Digits, this->Digits, b.Size, this->Size, true);
  }
  else
  {
    if (this->sign)
      return false;
    else
      return true;
  }
}



/* ##########################################################################
######   operator>                                                     ######
######                                                                 ######
######   Version: 19.12.2007                                           ######
########################################################################## */
bool CHugeInt::operator>(const CHugeInt &b) const
{
  if (this->sign == b.sign)
  {
    if (this->sign)
      return this->StandardLess(b.Digits, this->Digits, b.Size, this->Size, false);
    else
      return this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false);
  }
  else
  {
    if (this->sign)
      return true;
    else
      return false;
  }
}



/* ##########################################################################
######   operator<                                                     ######
######                                                                 ######
######   Version: 10.01.2008                                           ######
########################################################################## */
bool CHugeInt::operator<(const CHugeInt &b) const
{
  if (this->sign == b.sign)
  {
    if (this->sign)
      return this->StandardLess(this->Digits, b.Digits, this->Size, b.Size, false);
    else
      return this->StandardLess(b.Digits, this->Digits, b.Size, this->Size, false);
  }
  else
  {
    if (this->sign)
      return false;
    else
      return true;
  }
}


/* ##########################################################################
######   AddAllPositive(...)                                           ######
######                                                                 ######
######   changes: this->Size                                           ######
######            this->Digits                                         ######
######                                                                 ######
######   Version: 23.01.2008                                           ######
########################################################################## */
void CHugeInt::AddAllPositive(const vector<vector<unsigned char> > &Digits)
{
  const size_t s1 = Digits.size();
  if (s1 == 0)
  {
    this->Size = 1;
    this->Digits.assign(1,0);
    return;
  }

  this->Digits.clear();

  vector<bool> DigitsAlive(s1, true);
  size_t Number_Of_Alive = s1;

  Heap_CHugeInt::sum   = 0;
  Heap_CHugeInt::carry = 0;

  Heap_CHugeInt::i = 0;
  while (true)
  {
    Heap_CHugeInt::sum = Heap_CHugeInt::carry;
    for (Heap_CHugeInt::j = 0; Heap_CHugeInt::j < s1; ++Heap_CHugeInt::j)
    {
      if (DigitsAlive[Heap_CHugeInt::j])
      {
        if (Heap_CHugeInt::i < Digits[Heap_CHugeInt::j].size())
          Heap_CHugeInt::sum += Digits[Heap_CHugeInt::j][Heap_CHugeInt::i];
        else
        {
          --Number_Of_Alive;
          DigitsAlive[Heap_CHugeInt::j] = false;
        }
      }
    }
    if (Number_Of_Alive == 0)
    {
      while(Heap_CHugeInt::carry != 0)
      {
        this->Digits.push_back(Heap_CHugeInt::carry % 10);
        Heap_CHugeInt::carry /= 10;
      }
      this->Size = this->Digits.size();
      return;
    }

    if (Heap_CHugeInt::sum > 9)
    {
      this->Digits.push_back(Heap_CHugeInt::sum % 10);
      Heap_CHugeInt::carry = Heap_CHugeInt::sum/10;
    }
    else
    {
      this->Digits.push_back(Heap_CHugeInt::sum);
      Heap_CHugeInt::carry = 0;
    }
    ++Heap_CHugeInt::i;
  }
  cout << "Error in void CHugeInt::AddAllPositive(const vector<vector<unsigned char> > &Digits)" << endl;
  exit(1);
}



/* ##########################################################################
######   abs                                                           ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt CHugeInt::abs() const
{
  Heap_CHugeInt::result.Digits = this->Digits;
  Heap_CHugeInt::result.sign   = true;
  Heap_CHugeInt::result.zero   = this->zero;
  Heap_CHugeInt::result.Size   = this->Size;
  
  return Heap_CHugeInt::result;
}



/* ##########################################################################
######   StandardSubtraction(...)  only when a >= b                    ######
######                                                                 ######
######   changes: this->Size                                           ######
######            this->Digits                                         ######
######                                                                 ######
######   Version: 09.01.2008                                           ######
########################################################################## */
void CHugeInt::StandardSubtraction(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, bool check)
{
  if (check && this->StandardLess(DigitsA, DigitsB, DigitsA.size(), DigitsB.size(), false))
  {
    cout << "Error in void CHugeInt::StandardSubtraction(...). (a < b)" << endl;
    exit(1);
  }
  
  const size_t s1 = DigitsA.size();
  const size_t s2 = DigitsB.size();
  
  this->Digits.clear();  // does not change the capacity of the vector
  
  unsigned i = 0;
  
  char sum = 0;
  char carry = 0;
  for (i = 0; i < s2; ++i)
  {
    sum = DigitsA[i] - DigitsB[i] + carry;
    if (sum < 0)
    {
      this->Digits.push_back(sum+10);
      carry = -1;
    }
    else
    {
      this->Digits.push_back(sum);
      carry = 0;
    }
  }
  for (i = s2; i < s1; ++i)
  {
    sum = DigitsA[i] + carry;
    if (sum < 0)
    {
      this->Digits.push_back(sum+10);
      carry = -1;
    }
    else
    {
      this->Digits.push_back(sum);
      carry = 0;
    }
  }
  
  if ((this->Digits.size() == 1) && (this->Digits[0] == 0))
  {
    this->Size = 1;
    return;
  }
  
  while ((this->Digits.size() > 1) && (this->Digits[this->Digits.size()-1] == 0))
    this->Digits.pop_back();

  this->Size = this->Digits.size();
}



/* ##########################################################################
######   StandardAddition(...)   when a,b >= 0                         ######
######                                                                 ######
######   changes: this->Size                                           ######
######            this->Digits                                         ######
######                                                                 ######
######   Version: 22.01.2008                                           ######
########################################################################## */
void CHugeInt::StandardAddition(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, const size_t &SizeA, const size_t &SizeB)
{
  size_t minsize = SizeA;
  bool s1_smaller_s2 = true;
  
  if (SizeB < SizeA)
  {
    minsize = SizeB;
    s1_smaller_s2 = false;
  }
  
  this->Digits.clear();  // does not change the capacity of the vector
  
  unsigned i = 0;
  
  unsigned char sum = 0;
  bool carry = 0;
  
  for (i = 0; i < minsize; ++i)
  {
    sum = DigitsA[i] + DigitsB[i] + carry;
    if (sum > 9)
    {
      this->Digits.push_back(sum-10);
      carry = 1;
    }
    else
    {
      this->Digits.push_back(sum);
      carry = 0;
    }
  }
  
  if (SizeA == SizeB)
  {
    if (carry == 1)
      this->Digits.push_back(1);
      
    this->Size = this->Digits.size();
    return;
  }
    
  if (s1_smaller_s2)
  {
    for (i = SizeA; i < SizeB; ++i)
    {
      sum = DigitsB[i] + carry;
      if (sum > 9)
      {
        this->Digits.push_back(sum-10);
        carry = 1;
      }
      else
      {
        this->Digits.push_back(sum);
        carry = 0;
      }
    }
  }
  else
  {
    for (i = SizeB; i < SizeA; ++i)
    {
      sum = DigitsA[i] + carry;
      if (sum > 9)
      {
        this->Digits.push_back(sum-10);
        carry = 1;
      }
      else
      {
        this->Digits.push_back(sum);
        carry = 0;
      }
    }
  }
  if (carry == 1)
    this->Digits.push_back(1);
    
  this->Size = this->Digits.size();
}



/* ##########################################################################
######   StandardSubtraction(...)  only when a >= b                    ######
######                                                                 ######
######   changes: this->Size                                           ######
######            this->Digits                                         ######
######                                                                 ######
######   Version: 25.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
void CHugeInt::StandardSubtraction(const vector<unsigned char> &DigitsB, const size_t &SizeB)
{
  Heap_CHugeInt::c_sum   = 0;
  Heap_CHugeInt::b_carry = 0;
  
  for (Heap_CHugeInt::i = 0; Heap_CHugeInt::i < SizeB; ++Heap_CHugeInt::i)
  {
    Heap_CHugeInt::c_sum = this->Digits[Heap_CHugeInt::i] - DigitsB[Heap_CHugeInt::i] - Heap_CHugeInt::b_carry;
    if (Heap_CHugeInt::c_sum < 0)
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::c_sum+10;
      Heap_CHugeInt::b_carry = 1;
    }
    else
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::c_sum;
      Heap_CHugeInt::b_carry = 0;
    }
  }
  for (Heap_CHugeInt::i = SizeB; Heap_CHugeInt::i < this->Size; ++Heap_CHugeInt::i)
  {
    Heap_CHugeInt::c_sum = this->Digits[Heap_CHugeInt::i] - Heap_CHugeInt::b_carry;
    if (Heap_CHugeInt::c_sum < 0)
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::c_sum+10;
      Heap_CHugeInt::b_carry = 1;
    }
    else
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::c_sum;
      Heap_CHugeInt::b_carry = 0;
    }
  }
  
  for (Heap_CHugeInt::i = this->Size-1; Heap_CHugeInt::i > 0; --Heap_CHugeInt::i)
  {
    if (this->Digits[Heap_CHugeInt::i] == 0)
      this->Digits.pop_back();
    else
    {
      this->Size = Heap_CHugeInt::i + 1;
      return;
    }
  }
  this->Size = 1;
}



/* ##########################################################################
######   StandardAddition(...)   when a,b >= 0                         ######
######                                                                 ######
######   changes: this->Size                                           ######
######            this->Digits                                         ######
######                                                                 ######
######   Version: 25.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
void CHugeInt::StandardAddition(const vector<unsigned char> &DigitsB, const size_t &SizeB)
{
  size_t minsize = this->Size;
  bool SA_smaller_SB = true;

  if (SizeB < this->Size)
  {
    minsize = SizeB;
    SA_smaller_SB = false;
  }

  Heap_CHugeInt::uc_sum  = 0;
  Heap_CHugeInt::b_carry = 0;

  for (Heap_CHugeInt::i = 0; Heap_CHugeInt::i < minsize; ++Heap_CHugeInt::i)
  {
    Heap_CHugeInt::uc_sum = this->Digits[Heap_CHugeInt::i] + DigitsB[Heap_CHugeInt::i] + Heap_CHugeInt::b_carry;
    if (Heap_CHugeInt::uc_sum > 9)
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum-10;
      Heap_CHugeInt::b_carry = 1;
    }
    else
    {
      this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum;
      Heap_CHugeInt::b_carry = 0;
    }
  }

  if (this->Size == SizeB)
  {
    if (Heap_CHugeInt::b_carry == 1)
    {
      ++this->Size;
      this->Digits.push_back(1);
    }
    return;
  }

  if (SA_smaller_SB)
  {
    this->Digits.resize(SizeB,0);
    for (Heap_CHugeInt::i = this->Size; Heap_CHugeInt::i < SizeB; ++Heap_CHugeInt::i)
    {
      Heap_CHugeInt::uc_sum = DigitsB[Heap_CHugeInt::i] + Heap_CHugeInt::b_carry;
      if (Heap_CHugeInt::uc_sum > 9)
      {
        this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum-10;
        Heap_CHugeInt::b_carry = 1;
      }
      else
      {
        this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum;
        Heap_CHugeInt::b_carry = 0;
      }
    }
    this->Size = SizeB;
  }
  else
  {
    for (Heap_CHugeInt::i = SizeB; Heap_CHugeInt::i < this->Size; ++Heap_CHugeInt::i)
    {
      Heap_CHugeInt::uc_sum = this->Digits[Heap_CHugeInt::i] + Heap_CHugeInt::b_carry;
      if (Heap_CHugeInt::uc_sum > 9)
      {
        this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum-10;
        Heap_CHugeInt::b_carry = 1;
      }
      else
      {
        this->Digits[Heap_CHugeInt::i] = Heap_CHugeInt::uc_sum;
        Heap_CHugeInt::b_carry = 0;
      }
    }
  }
  if (Heap_CHugeInt::b_carry == 1)
  {
    ++this->Size;
    this->Digits.push_back(1);
  }
}



/* ##########################################################################
######   StandardLess   A < B   when A,B >= 0                          ######
######                                                                 ######
######   Version: 24.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
bool CHugeInt::StandardLess(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, const size_t &SizeA, const size_t &SizeB, bool TrueIfEqual) const
{
  if (SizeA < SizeB)
    return true;
  
  if (SizeA == SizeB)
  {
    for (Heap_CHugeInt::int_i = SizeA-1; Heap_CHugeInt::int_i >= 0; --Heap_CHugeInt::int_i)
    {
      if (DigitsA[Heap_CHugeInt::int_i] < DigitsB[Heap_CHugeInt::int_i])
        return true;
      else
      {
        if (DigitsA[Heap_CHugeInt::int_i] != DigitsB[Heap_CHugeInt::int_i])
          return false;
      }
    }
    // both numbers are equal
    return TrueIfEqual;
  }
  return false;
}



/* ##########################################################################
######   long long int CHugeInt::ToLongLongInt() const                 ######
######                                                                 ######
######   Version: 22.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
long long int CHugeInt::ToLongLongInt() const
{
  long long int result = 0;
  long long int b = 1;

  for (unsigned i = 0; i < this->Size; ++i)
  {
    result += (this->Digits[i] * b);
    b *= 10;
  }

  if (!this->sign)
    result *= -1;

  return result;
}


long long int ToLongLongInt(unsigned a)
{
  return a;
}
  
long long int ToLongLongInt(CHugeInt a)
{
  return a.ToLongLongInt();
}

/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 09.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt::CHugeInt()
{
  this->sign = true;
  this->zero = true;
  this->Size = 1;
  this->Digits.push_back(0);
}



/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 22.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt::CHugeInt(int a)
{
  if (a == 0)
  {
    this->sign = true;
    this->zero = true;
    this->Size = 1;
    this->Digits.push_back(0);
    return;
  }
  
  this->sign = (a > 0);  // a != 0
  this->zero = false;
  
  if (!this->sign)  // a < 0
    a *= -1;
  
  this->Digits.push_back(a % 10);
  a /= 10;
  while(a != 0)
  {
    this->Digits.push_back(a % 10);
    a /= 10;
  }
  this->Size = this->Digits.size();
}



/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 11.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt::CHugeInt(unsigned a)
{
  this->sign = true;

  if (a == 0)
  {
    this->zero = true;
    this->Size = 1;
    this->Digits.push_back(0);
    return;
  }
  this->zero = false;
  
  this->Digits.push_back(a % 10);
  a /= 10;
  while(a != 0)
  {
    this->Digits.push_back(a % 10);
    a /= 10;
  }
  this->Size = this->Digits.size();
}


/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 11.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt::CHugeInt(unsigned long long a)
{
  this->sign = true;

  if (a == 0)
  {
    this->zero = true;
    this->Size = 1;
    this->Digits.push_back(0);
    return;
  }
  this->zero = false;
  
  this->Digits.push_back(a % 10);
  a /= 10;
  while(a != 0)
  {
    this->Digits.push_back(a % 10);
    a /= 10;
  }
  this->Size = this->Digits.size();
}


/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 23.01.2008                                           ######
######   checked: 06.02.2008                                           ######
########################################################################## */
CHugeInt::CHugeInt(const CHugeInt &tmp)
 : Digits(tmp.Digits), sign(tmp.sign), zero(tmp.zero), Size(tmp.Size)
{
}


/* ##########################################################################
######   CHugeInt()                                                    ######
######                                                                 ######
######   Version: 28.02.2008                                           ######
######                                                                 ######
########################################################################## */
CHugeInt::CHugeInt(const string &tmp)
{
  int to = 0;
  string Minus = "-";

  if (tmp[0] == Minus[0])
  {
    this->sign = false;
    to = 1;
  }
  else
    this->sign = true;

  for (int i = tmp.size()-1; i >= to; --i)
    this->Digits.push_back( atoi(tmp.substr(i,1).c_str()) );

  this->Size = this->Digits.size();
  if ((this->Size == 1) && (this->Digits[0] == 0))
    this->zero = true;
  else
    this->zero = false;
}



/* ##########################################################################
######   ~CHugeInt()                                                   ######
######                                                                 ######
######   Version: 18.12.2007                                           ######
########################################################################## */
CHugeInt::~CHugeInt()
{
}

bool CHugeInt::Check() const
{
  if ((this->Size > 1) && (this->Digits[this->Size-1] == 0))
    return false;
  else
    return true;
}
