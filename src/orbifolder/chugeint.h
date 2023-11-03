#ifndef CHUGEINT_H
#define CHUGEINT_H

#include <vector>
#include <iostream>  
#include <string>
#include <limits>
#include <cmath>

#include <boost/config.hpp>
#include <boost/rational.hpp>

using std::string;
using std::ostream;
using std::vector;
using boost::rational;


class CHugeInt{
public:
    CHugeInt();
    CHugeInt(int a);
    CHugeInt(unsigned a);
    CHugeInt(unsigned long long a);
    CHugeInt(const CHugeInt &tmp);
    CHugeInt(const string &tmp);

    ~CHugeInt();

    CHugeInt operator+(const CHugeInt &b) const;
    CHugeInt operator-(const CHugeInt &b) const;
    CHugeInt operator*(const CHugeInt &b) const;
    CHugeInt operator/(const CHugeInt &b) const;
    CHugeInt operator%(const CHugeInt &b) const;
    
    void operator+=(const CHugeInt &b);
    void operator-=(const CHugeInt &b);
    void operator*=(const CHugeInt &b);
    void operator/=(const CHugeInt &b);
    void operator%=(const CHugeInt &b);

    void operator++();
    void operator--();
    
    CHugeInt operator+() const;
    CHugeInt operator-() const;
    
    bool operator!() const;
    
    void operator=(const CHugeInt &b);
    void operator=(const int &a);
    void operator=(const unsigned &a);
    
    bool operator==(const CHugeInt &b) const;
    bool operator!=(const CHugeInt &b) const;
    bool operator>=(const CHugeInt &b) const;
    bool operator<=(const CHugeInt &b) const;
    bool operator>(const CHugeInt &b) const;
    bool operator<(const CHugeInt &b) const;
    
    CHugeInt abs() const;
    long long int ToLongLongInt() const;
    
    void AddAllPositive(const vector<vector<unsigned char> > &Digits);
    
    void StandardAddition(const vector<unsigned char> &DigitsB, const size_t &SizeB);
    void StandardSubtraction(const vector<unsigned char> &DigitsB, const size_t &SizeB);
    
    void StandardAddition(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, const size_t &SizeA, const size_t &SizeB);
    void StandardSubtraction(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, bool check = false);
    
    bool StandardLess(const vector<unsigned char> &DigitsA, const vector<unsigned char> &DigitsB, const size_t &SizeA, const size_t &SizeB, bool TrueIfEqual) const;

    bool Check() const;

    vector<unsigned char> Digits;
    bool sign;
    bool zero;
    size_t Size;
};

ostream& operator << (ostream &os, const CHugeInt &a);
rational<CHugeInt> D2RatHugeInt(double x);
double RatHugeInt2D(const rational<CHugeInt> &r);
rational<int> RatHugeInt2Rat(const rational<CHugeInt> &r);
long long int ToLongLongInt(unsigned a);
long long int ToLongLongInt(CHugeInt a);


namespace std
{
  template<>
  class numeric_limits<CHugeInt>
  {
    public:
      static const bool is_specialized = true;
      static CHugeInt min() {return CHugeInt(-2099999999);};
      static CHugeInt max() {return CHugeInt(2099999999);};
      static const int  digits = 0;
      static const int  digits10 = 0;
      static const bool is_signed = true;
      static const bool is_integer = true;
      static const bool is_exact = true;
      static const int radix = 0;
      static CHugeInt epsilon() throw();
      static CHugeInt round_error() throw();
      
      static const int  min_exponent = 0;
      static const int  min_exponent10 = 0;
      static const int  max_exponent = 0;
      static const int  max_exponent10 = 0;
      
      static const bool has_infinity = false;
      static const bool has_quiet_NaN = false;
      static const bool has_signaling_NaN = false;
      //static const float_denorm_style has_denorm = denorm absent;
      static const bool has_denorm_loss = false;
      static CHugeInt infinity() throw();
      static CHugeInt quiet_NaN() throw();
      static CHugeInt signaling_NaN() throw();
      static CHugeInt denorm_min() throw();
      
      static const bool is_iec559 = false;
      static const bool is_bounded = false;
      static const bool is_modulo = false;
      
      static const bool traps = false;
      static const bool tinyness_before = false;
      //static const float_round_style round_style = round_toward_zero;
      
  };
}

#endif
