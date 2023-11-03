#include "clinalg.h"
#include "linalg.hpp"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "globalfunctions.h"

using std::cout;
using std::endl;
using std::exit;

template <class T> CLinAlg<T>::CLinAlg()
 : One(T(1),T(1)), Zero(T(0),T(1)) //: One(1,1), Zero(0,1)
{
}


template <class T> CLinAlg<T>::~CLinAlg()
{
}



/* ###########################################################################
######   bool criterion_number_of_zeros(...)                            ######
########################################################################### */
template <class T> bool criterion_number_of_zeros(const vector<rational<T> > &v1, const vector<rational<T> > &v2)
{
  const size_t s1 = v1.size();
  if (s1 != v2.size())
  {
    cout << "\n  Warning in bool criterion_number_of_zeros(...). Exit." << endl;
    return false;
  }

  unsigned zeros_v1 = 0;
  unsigned zeros_v2 = 0;

  for (unsigned i = 0; i < s1; ++i)
  {
    if (v1[i] == T(0))
      ++zeros_v1;
    if (v2[i] == T(0))
      ++zeros_v2;
  }
  return (zeros_v1 > zeros_v2);
}


/* ########################################################################################
######   void CreateSolutionsFromCone(...)                                           ######
######                                                                               ######
######   Version: 14.10.2008                                                         ######
######################################################################################## */
/*void CreateSolutionsFromCone(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<vector<unsigned> > &AllNumbers, const vector<vector<rational<CHugeInt> > > &UndefBasisOfKernel, const vector<vector<rational<CHugeInt> > > &Cone)
{
  const size_t number_of_Digits = currentNumber.size();

  rational<CHugeInt> Zero;
  // runs through the positions of the current number
  for (unsigned i = 0; i < currentMaxDigit; ++i)
  {
    currentNumber[currentPosition] = i;

    unsigned nextPosition = currentPosition + 1;

    if (nextPosition < number_of_Digits)
      CreateSolutionsFromCone(currentNumber, nextPosition, MaxDigits[nextPosition], MaxDigits, AllNumbers, UndefBasisOfKernel, Cone);
    else
    {
      const size_t q1 = Cone.size();
      const size_t k1 = Cone[0].size();
      const size_t s2 = UndefBasisOfKernel[0].size();

      vector<rational<CHugeInt> > TestVector(k1, CHugeInt(0));
      vector<rational<CHugeInt> > Solution(s2, CHugeInt(0));

      unsigned j = 0;
      unsigned k = 0;

      for (j = 0; j < q1; ++j)
      {
        const vector<rational<CHugeInt> > &ConeVector = Cone[j];
        for (k = 0; k < k1; ++k)
          TestVector[k] = TestVector[k] + (currentNumber[j] * ConeVector[k] / (CHugeInt(10)));
      }

      bool all_zero = true;
      for (j = 0; j < s2; ++j)
      {
        for (k = 0; k < k1; ++k)
          Solution[j] += UndefBasisOfKernel[k][j] * TestVector[k];
        if (Solution[j] > Zero)
          all_zero = false;
        if (Solution[j] < Zero)
        {
          cout << "\n  Warning in void CreateSolutionsFromCone(...) of type CHugeInt: solution is wrong. Return." << endl;
          return;
        }
      }
      if (!all_zero)
      {
        const size_t s0 = AllNumbers.size();

        bool solution_unknown = true;
        if (s0 != 0)
        {
          unsigned l = 0;

          vector<unsigned> currentNumber2(s0, 0);
          vector<unsigned> MaxDigits2(s0, 2);
          vector<vector<unsigned> > AllNumbers2;
          RecursiveCounting(currentNumber2, 0, MaxDigits2[0], MaxDigits2, AllNumbers2);

          const size_t s1 = AllNumbers2.size();

          vector<unsigned> tmp(number_of_Digits, 0);
          for (j = 0; solution_unknown && (j < s1); ++j)
          {
            tmp.assign(number_of_Digits, 0);

            const vector<unsigned> &AllNumber2 = AllNumbers2[j];

            for (k = 0; k < s0; ++k)
            {
              for (l = 0; l < number_of_Digits; ++l)
                tmp[l] += AllNumber2[k] * AllNumbers[k][l];
            }
            if (currentNumber == tmp)
              solution_unknown = false;
          }
        }
        if (solution_unknown)
        {
          cout << "save Solution!" << endl;
          AllNumbers.push_back(currentNumber);
        }
      }
    }
  }
}*/


/*
void CreateSolutionsFromCone(vector<unsigned> &currentNumber, unsigned currentPosition, unsigned currentMaxDigit, vector<unsigned> &MaxDigits, vector<vector<unsigned> > &AllNumbers, const vector<vector<rational<int> > > &UndefBasisOfKernel, const vector<vector<rational<int> > > &Cone)
{
  const size_t number_of_Digits = currentNumber.size();

  // runs through the positions of the current number
  for (unsigned i = 0; i < currentMaxDigit; ++i)
  {
    currentNumber[currentPosition] = i;

    unsigned nextPosition = currentPosition + 1;

    if (nextPosition < number_of_Digits)
      CreateSolutionsFromCone(currentNumber, nextPosition, MaxDigits[nextPosition], MaxDigits, AllNumbers, UndefBasisOfKernel, Cone);
    else
    {
      const size_t q1 = Cone.size();
      const size_t k1 = Cone[0].size();
      const size_t s2 = UndefBasisOfKernel[0].size();

      vector<rational<int> > TestVector(k1, 0);
      vector<rational<int> > Solution(s2, 0);

      unsigned j = 0;
      unsigned k = 0;

      for (j = 0; j < q1; ++j)
      {
        const vector<rational<int> > &ConeVector = Cone[j];
        for (k = 0; k < k1; ++k)
          TestVector[k] = TestVector[k] + (currentNumber[j] * ConeVector[k] / 10);
      }

      bool all_zero = true;
      for (j = 0; j < s2; ++j)
      {
        for (k = 0; k < k1; ++k)
          Solution[j] += UndefBasisOfKernel[k][j] * TestVector[k];
        if (Solution[j] > 0)
          all_zero = false;
        if (Solution[j] < 0)
        {
          cout << "\n  Warning in void CreateSolutionsFromCone(...) of type int: solution is wrong!" << endl;
          exit(1);
        }
      }
      if (!all_zero)
      {
        const size_t s0 = AllNumbers.size();

        stable_sort(AllNumbers.begin(), AllNumbers.end(), criterion4);

        bool solution_unknown = true;
        if (s0 != 0)
        {
          unsigned l = 0;

          vector<unsigned> currentNumber2(s0, 0);
          vector<unsigned> MaxDigits2(s0, 2);
          vector<vector<unsigned> > AllNumbers2;
          RecursiveCounting(currentNumber2, 0, MaxDigits2[0], MaxDigits2, AllNumbers2);

          const size_t s1 = AllNumbers2.size();

          vector<unsigned> tmp(number_of_Digits, 0);
          for (j = 0; solution_unknown && (j < s1); ++j)
          {
            tmp.assign(number_of_Digits, 0);

            const vector<unsigned> &AllNumber2 = AllNumbers2[j];

            for (k = 0; k < s0; ++k)
            {
              for (l = 0; l < number_of_Digits; ++l)
                tmp[l] += AllNumber2[k] * AllNumbers[k][l];
            }
            if (currentNumber == tmp)
              solution_unknown = false;
          }
        }
        if (solution_unknown)
        {
          cout << "save Solution! ";
          for (unsigned q1 = 0; q1 < currentNumber.size(); ++q1)
            cout << currentNumber[q1] << " ";
          cout << endl;
          AllNumbers.push_back(currentNumber);
        }
      }
    }
  }
}

*/

/* ##########################################################################
######   bool CLinAlg::FindKernel(...)                                 ######
######                                                                 ######
######   Version: 22.09.2011                                           ######
########################################################################## */
template <class T> bool CLinAlg<T>::FindKernel(const vector<vector<rational<T> > > &M, vector<vector<rational<T> > > &BasisOfKernel, bool TestResult)
{
  /*cout << "Matrix M\n";
  for (unsigned i = 0 ; i < M.size(); ++i)
  {
    for (unsigned j = 0 ; j < M[i].size(); ++j)
      cout << M[i][j] << "  ";
    cout << "\n";
  }
  cout << endl;*/

  vector<vector<rational<T> > > Matrix = M;

  size_t s1 = Matrix.size();
  if (s1 == 0)
  {
    cout << "\n  Warning in bool CLinAlg::FindKernel(...): Matrix has 0 rows. Return false." << endl;
    return false;
  }
  size_t s2 = Matrix[0].size();

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  // begin: make M quadratic if neccessary
  if (s1 != s2)
  {
    rational<T> sum = this->Zero;
    vector<rational<T> >          s2Line(s2, this->Zero);
    vector<vector<rational<T> > > QuadraticM(s2, s2Line);

    for (i = 0; i < s2; ++i)
    {
      for (j = 0; j < s2; ++j)
      {
        sum = this->Zero;
        for (k = 0; k < s1; ++k)
          sum += M[k][i] * M[k][j];
        QuadraticM[i][j] = sum;
      }
    }
    Matrix = QuadraticM;

    s1 = Matrix.size();
    s2 = Matrix[0].size();
  }
  // end: make M quadratic if neccessary
  if (!this->ReducedRowEchelonForm(Matrix))
  {
    cout << "\n  Warning in bool CLinAlg::FindKernel(...) : ReducedRowEchelonForm failed. Return false." << endl;
    return false;
  }

  /*cout << "Matrix\n";
  for (unsigned i = 0 ; i < Matrix.size(); ++i)
  {
    for (unsigned j = 0 ; j < Matrix[i].size(); ++j)
      cout << Matrix[i][j] << "  ";
    cout << "\n";
  }
  cout << endl;*/

  if (BasisOfKernel.size() != 0)
  {
    cout << "\n  Warning in bool CLinAlg::FindKernel(...): BasisOfKernel is now cleared!" << endl;
    BasisOfKernel.clear();
  }

  // begin: find the positions of the 1
  vector<unsigned> Pos_of_entries;

  bool Found_one = false;

  // run through the rows of the matrix
  for (i = 0; i < s1; ++i)
  {
    Found_one = false;
    // run through the columns of the matrix starting on the diagonal
    for (j = i; !Found_one && (j < s2); ++j)
    {
      // search for the first 1
      if (Matrix[i][j] == this->One)
      {
        Pos_of_entries.push_back(j);
        Found_one = true;
      }
    }
  }
  // end: find the positions of the 1

  const size_t s3 = Pos_of_entries.size();

  // begin: find the non-trivial columns
  bool column_nontrivial = false;

  // run through the columns of the matrix
  for (i = 0; i < s2; ++i)
  {
    column_nontrivial = false;

    if (find(Pos_of_entries.begin(),Pos_of_entries.end(),i) == Pos_of_entries.end())
    {
      vector<rational<T> > BasisVector(s2, this->Zero);

      unsigned line = 0;
      for (j = 0; j < s3; ++j)
      {
        BasisVector[Pos_of_entries[j]] = -Matrix[line][i];
        ++line;
      }
      BasisVector[i] = this->One;

      BasisOfKernel.push_back(BasisVector);
    }
  }
  // end: find the non-trivial columns


  // begin: test the dimension of the kernel
  const size_t s4 = BasisOfKernel.size();

  /*vector<rational<int> > row(s2,0);
  vector<vector<rational<int> > > RatIntMatrix(s1, row);
  for (i = 0; i < s1; ++i)
  {
    for (j = 0; j < s2; ++j)
      RatIntMatrix[i][j] = rational<int>(ToLongLongInt(Matrix[i][j].numerator()),ToLongLongInt(Matrix[i][j].denominator()));
  }
    cout << "G" << endl;
  const vector<vector<rational<int> > > Basis = findBasis<rational<int> >(RatIntMatrix);
    cout << "H" << endl;
    
  if (Basis.size() + s4 != s2)
  {
    if (findBasis<rational<int> >(RatIntMatrix).size() + s4 != s2)
    {
      cout << "\n  Warning in bool CLinAlg::FindKernel(...): dim(im M) + dim(ker M) != dim(M):\ndim(im M)  = " << Basis.size() << "\ndim(ker M) = " << s4 << "\ndim(M)     = " << s2 << "\nQuadratic Matrix:\n";
      for (i = 0 ; i < RatIntMatrix.size(); ++i)
      {
        for (j = 0 ; j < RatIntMatrix[i].size(); ++j)
          cout << RatIntMatrix[i][j] << "  ";
        cout << "\n";
      }
      cout << endl;
      exit(1);
    }
  }*/
  // end: test the dimension of the kernel

  if (!TestResult)
    return true;

  // begin: test the kernel
  s1 = M.size();
  s2 = M[0].size();
  vector<rational<T> > NullVector(s2, this->Zero);

  for (i = 0; i < s4; ++i)
  {
    const vector<rational<T> > &Basis = BasisOfKernel.at(i);

    if (Basis.size() != s2)
    {
      cout << "\n  Warning in bool CLinAlg::FindKernel(...). Return false." << endl;
      return false;
    }

    for (j = 0; j < s1; ++j)
    {
      rational<T> sum = this->Zero;
      for (k = 0; k < s2; ++k)
        sum += (M.at(j).at(k) * Basis.at(k));

      if (sum != this->Zero)
      {
        cout << "\n  Warning in bool CLinAlg::FindKernel(...): Solution is not in the kernel:" << endl;
        for (k = 0; k < s2; ++k) 
          cout << Basis[k] << "  ";
        cout << " Return false." << endl;
        return false;
      }
    }
  }
  // end: test the kernel

  return true;
}



/* ########################################################################################
######   bool CLinAlg::FindPositiveIntegerKernel(...)                                ######
######                                                                               ######
######   Version: 22.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) M : matrix M whose kernel (positive and integer) is looked for           ######
######   output:                                                                     ######
######   1) BasisOfKernel : a basis of kernel vectors with non-negative entries      ######
######################################################################################## */
template <class T> bool CLinAlg<T>::FindPositiveIntegerKernel(const vector<vector<rational<T> > > &M, vector<vector<rational<T> > > &BasisOfKernel)
{
  if (BasisOfKernel.size() != 0)
  {
    cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): BasisOfKernel was not empty. Now cleared." << endl;
    BasisOfKernel.clear();
  }

  // begin: find kernel
  vector<vector<rational<T> > > UndefBasisOfKernel;

  //cout << "Now starting FindKernel..." << endl;

  if (!this->FindKernel(M, UndefBasisOfKernel, true))
  {
    cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): cannot find the kernel of the matrix. Return false." << endl;
    return false;
  }
  // end: find kernel

  const size_t k1 = UndefBasisOfKernel.size();
  //cout << "kernel is spanned by " << k1 << " vectors." << endl;

  if (k1 == 0)
    return true;

  const size_t s1 = M.size();
  const size_t s2 = M[0].size();

  if (s2 != UndefBasisOfKernel[0].size())
  {
    cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): check the dimensions. Return false." << endl;
    return false;
  }

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  /*cout << "UndefBasisOfKernel: \n";
  for (unsigned i = 0; i < k1; ++i)
  {
    for (unsigned k = 0; k < s2; ++k)
      cout << UndefBasisOfKernel[i][k] << " ";
    cout << endl;
  }*/

  // simplify the inequalities (by making linear combinations with positive coefficients)
  this->SimplifyInequality(UndefBasisOfKernel);

  /*cout << "UndefBasisOfKernel (simplified): \n";
  for (unsigned i = 0; i < k1; ++i)
  {
    for (unsigned k = 0; k < s2; ++k)
      cout << UndefBasisOfKernel[i][k] << " ";
    cout << endl;
  }*/

  // begin: transpose the matrix of basis vectors
  //        UndefBasisOfKernel is the matrix with i-th row b_i,
  //        where b_i is from the kernel (M b_i = 0)
  vector<rational<T> >          k1Line(k1, T(0));
  vector<vector<rational<T> > > UndefBasisOfKernelTranspose(s2, k1Line);
  for (i = 0; i < k1; ++i)
  {
    for (j = 0; j < s2; ++j)
      UndefBasisOfKernelTranspose[j][i] = UndefBasisOfKernel[i][j];
  }
  // end: transpose the matrix of basis vectors

  // begin: find solutions x to the inequality
  //        sum_i b_i x_i > (0,..,0)^T
  //        where the vectors x span a cone

  //cout << "Now starting SolveSystemOfLinearInequalities..." << endl;

  vector<vector<rational<T> > > Cone;
  this->SolveSystemOfLinearInequalities(UndefBasisOfKernelTranspose, Cone);
  // end: find solutions to the inequality

  const size_t q1 = Cone.size();
  if (q1 == 0)
  {
    //cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): cone is empty!" << endl;
    return false;
  }

  //cout << "Cone is spanned by " << q1 << " vectors." << endl;

  rational<T> max(T(1)); //rational<T> max(1);
  rational<T> sp(T(0)); //rational<T> sp(0);
  bool all_zero = true;

  // run through the basis vectors of the cone
  for (i = 0; i < q1; ++i)
  {
    const vector<rational<T> > &ConeVector = Cone[i];

    all_zero = true;
    // begin: create positive kernel-vector (called Solution) using the cone vector
    vector<rational<T> > Solution(s2, T(0));
    for (j = 0; j < s2; ++j)
    {
      for (k = 0; k < k1; ++k)
        Solution[j] += UndefBasisOfKernel[k][j] * ConeVector[k];
      if (Solution[j] > T(0))
        all_zero = false;
      if (Solution[j] < T(0))
      {
        cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): solution is wrong. Return false." << endl;
        return false;
      }
    }
    // end: create positive kernel-vector (called Solution) using the cone vector

    if (!all_zero)
    {
      // begin: make the solution integer
      max = T(2);
      while (max != T(1))
      {
        max = T(1);
        for (j = 0; j < s2; ++j)
        {
          if (Solution[j].denominator() > max)
            max = Solution[j].denominator();
        }
        if (max != T(1))
        {
          for (j = 0; j < s2; ++j)
            Solution[j] *= max;
        }
      }
      // end: make the solution integer

      // begin: check the solution
      for (j = 0; j < s1; ++j)
      {
        sp = T(0);
        for (k = 0; k < s2; ++k)
          sp += M[j][k] * Solution[k];

        if (sp != T(0))
        {
          cout << "\n  Warning in bool CLinAlg::FindPositiveIntegerKernel(...): solution is not in the kernel. Return false." << endl;
          return false;
        }
      }
      // end: check the solution
      BasisOfKernel.push_back(Solution);
    }
  }
  // end: check the solution

  //cout << "the positive kernel is spanned by " << BasisOfKernel.size() << " vectors." << endl;

  return true;
}



/* ##########################################################################
######   bool CLinAlg::GetParticularSolution(...)                      ######
######                                                                 ######
######   Version: 02.10.2008                                           ######
########################################################################## */
template <class T> bool CLinAlg<T>::GetParticularSolution(const vector<vector<rational<T> > > &A, const vector<rational<T> > &b, vector<rational<T> > &ParticularSolution)
{
  unsigned i, j;

  const size_t s1 = A.size();
  const size_t s2 = A[0].size();

  // begin: print A
  /*cout << "{";
  for (i = 0 ; i < s1; ++i)
  { 
    cout << "{";
    for (j = 0 ; j < s2-1; ++j)
      cout << A[i][j].numerator().ToLongLongInt() << "/" << A[i][j].denominator().ToLongLongInt() << ", ";

    cout << A[i][s2-1].numerator().ToLongLongInt() << "/" << A[i][s2-1].denominator().ToLongLongInt() << "},\n ";
  }
  cout << "}" << endl;*/
  // end: print A

  // begin: check if particular solution exists
  {
    doubleMatrix tmp_A(s1), tmp_Ab(s1);
    double tmp = 0.0;
    for(i = 0; i < s1; ++i)
    {
      const vector<rational<T> > &Ai = A[i];
      const rational<T>          &bi = b[i];
      for(j = 0; j < s2; ++j)
      {
        const rational<T> &Aij = Ai[j];
        tmp = ((double)ToLongLongInt(Aij.numerator()))/((double)ToLongLongInt(Aij.denominator()));
        tmp_A[i].push_back(tmp);
        tmp_Ab[i].push_back(tmp);
      }
      tmp_Ab[i].push_back(((double)ToLongLongInt(bi.numerator()))/((double)ToLongLongInt(bi.denominator())));
    }

    if(findBasis<double>(tmp_A).size() != findBasis<double>(tmp_Ab).size())
    {
      //cout << "Warning: There is no particular solution!" << endl;
      return false;
    }
  }
  // end: check if particular solution exists

  // add the column b to the matrix A to get the matrix (A|b)
  vector<vector<rational<T> > > Ab = A;
  for (i = 0; i < s1; ++i)
    Ab[i].push_back(b[i]);

  if (!this->ReducedRowEchelonForm(Ab))
  {
    cout << "\n  Warning in bool CLinAlg::GetParticularSolution(...): No particular solution found, but it should exist!" << endl;
    return false;
  }

  ParticularSolution.resize(s2, T(0));

  //cout <<"Particular Solution:\n";
  unsigned row = 0;
  for (i = 0; i < s2; ++i)
  {
    if (row < Ab.size())
    {
      if (Ab[row][i] == T(1))  //if (Ab[row][i] == 1)
      {
        ParticularSolution[i] = Ab[row][s2];
        ++row;
      }
      else
        ParticularSolution[i] = T(0);
    }
    //cout << ParticularSolution[i] << ", ";
  }
  //cout << endl;

  // begin: check A x = b
  vector<rational<T> > Ax;
  rational<T> tmp(T(0)); //rational<T> tmp(0);
  for (i = 0 ; i < s1; ++i)
  {
    tmp = this->Zero; //tmp = rational<T>(0);
    for (j = 0 ; j < s2; ++j)
      tmp += A[i][j] * ParticularSolution[j];

    if (tmp != b[i])
    {
      cout << "\n  Warning in bool CLinAlg::GetParticularSolution(...): solution is wrong!" << endl;
      return false;
    }
  }
  // end: check A x = b
  return true;
}


/* ########################################################################################
######   CLinAlg<T>::GramSchmidt(...)                                                ######
######                                                                               ######
######   Version: 12.04.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
template <class T> bool CLinAlg<T>::GramSchmidt(vector<vector<rational<T> > > Basis, vector<vector<rational<T> > > &OrthogonalBasis, bool TestResult)
{
  const size_t s1 = Basis.size();

  if (OrthogonalBasis.size() != 0)
  {
    cout << "\n  Warning in CLinAlg<T>::GramSchmidt(...) : OrthogonalBasis was not empty. Now cleared." << endl;
    OrthogonalBasis.clear();
  }

  if (s1 < 2)
  {
    OrthogonalBasis.insert(OrthogonalBasis.end(), Basis.begin(), Basis.end());
    return true;
  }

  stable_sort(Basis.begin(), Basis.end(), criterion_number_of_zeros<T>);

  const size_t d = Basis[0].size();
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  rational<T> tmp = T(0);
  vector<rational<T> > Length_of_OrthogonalBasisVectors;
  vector<rational<T> > NewVector;
  for (i = 0; i < s1; ++i)
  {
    const vector<rational<T> > &BasisVector = Basis[i];
    NewVector = BasisVector;
    for (j = 0; j < i; ++j)
    {
      const vector<rational<T> > &OrthogonalBasisVector = OrthogonalBasis[j];

      tmp = T(0);
      for (k = 0; k < d; ++k)
        tmp += BasisVector[k] * OrthogonalBasisVector[k];

      for (k = 0; k < d; ++k)
        NewVector[k] -= (OrthogonalBasisVector[k] * tmp)/Length_of_OrthogonalBasisVectors[j];
    }
    tmp = T(0);
    for (k = 0; k < d; ++k)
      tmp += NewVector[k] * NewVector[k];

    Length_of_OrthogonalBasisVectors.push_back(tmp);
    OrthogonalBasis.push_back(NewVector);
  }

  if (!TestResult)
    return true;

  for (i = 0; i < s1; ++i)
  {
    const vector<rational<T> > &Vector1 = OrthogonalBasis[i];

    for (j = i+1; j < s1; ++j)
    {
      const vector<rational<T> > &Vector2 = OrthogonalBasis[j];

      tmp = T(0);
      for (k = 0; k < d; ++k)
        tmp += Vector1[k] * Vector2[k];
      if (tmp != T(0))
      {
        cout << "\n  Warning in CLinAlg<T>::GramSchmidt(...) : Gram-Schmidt failed to produce orthogonal basis. Return false." << endl;
        return false;
      }
    }
  }
  return true;
}


/* ##########################################################################
######   bool CLinAlg::ReducedRowEchelonForm(...)                      ######
######                                                                 ######
######   Version: 15.10.2008                                           ######
########################################################################## */
template <class T> bool CLinAlg<T>::ReducedRowEchelonForm(vector<vector<rational<T> > > &M)
{
  const size_t s1 = M.size();

  if (s1 < 2)
  {
    cout << "\n  Warning in bool CLinAlg::ReducedRowEchelonForm(...): Matrix is too small." << endl;
    return false;
  }

  const size_t s2 = M[0].size();

  unsigned i = 0;
  unsigned j = 0;

  rational<T> Entry = this->Zero; //rational<T> Entry(0,1);
  rational<T> tmp   = this->Zero; //rational<T> tmp(0,1);

  unsigned lead = 0;

  // run through the rows of the matrix
  for (unsigned row = 0; row < s1; ++row)
  {
    i = row;
    // begin: search for a row with a non-zero entry at the lead position
    while (M[i][lead] == this->Zero)
    {
      if (++i >= s1)
      {
        i = row;
        ++lead;
        if(lead >= s2)
          return true;
      }
    }
    // end: search for a row with a non-zero entry at the lead position


    // exchange the rows
    if (i != row)
      swap(M[i], M[row]);

    // begin: bring the first non-zero entry of the current row to 1
    rational<T> EntryInverse = rational<T>(M[row][lead].denominator(), M[row][lead].numerator());
    if (EntryInverse != this->One)
    {
      for (j = 0; j < s2; ++j)
        M[row][j] *= EntryInverse;
    }
    // end: bring the first non-zero entry of the current row to 1

    // begin: use the 1 of the current row to get zeros in the other rows
    const vector<rational<T> > &M_row = M[row];
    for (i = 0; i < s1; ++i)
    {
      if (i != row)
      {
        Entry = M[i][lead];
        if (Entry.denominator() < 0)
        {
          cout << "\n  Warning in template <class T> bool CLinAlg<T>::ReducedRowEchelonForm(...): the denominator is smaller zero. Return false." << endl;
          return false;
        }
        if (Entry != this->Zero)
        {
          for (j = 0; j < s2; ++j)
            M[i][j] = M[i][j] - (Entry * M_row[j]);
        }
      }
    }
    // end: use the 1 of the current row to get zeros in the other rows
    ++lead;

    if (lead >= s2)
      return true;
  }
  return true;
}



/* ########################################################################################
######   FindOrdersOfWilsonLines(const vector<vector<rational<T> > > &input, ...)    ######
######                                                                               ######
######   Version: 08.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input : a matrix of relations between Wilson lines                       ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
template <class T> bool CLinAlg<T>::FindOrdersOfWilsonLines(const vector<vector<rational<T> > > &input, vector<vector<rational<T> > > &output)
{
  output = input;
  const size_t s1 = output.size();

  if (s1 < 2)
  {
    cout << "\n  Warning in bool CLinAlg::FindOrdersOfWilsonLines(...): Matrix is too small." << endl;
    return false;
  }

  const size_t s2 = output[0].size();

  unsigned i = 0;
  unsigned j = 0;

  rational<T> Entry = this->Zero; //rational<T> Entry(0,1);
  rational<T> tmp   = this->Zero;  //rational<T> tmp(0,1);

  unsigned lead = 0;
  rational<T> min(T(9999),T(1)); //rational<T> min(9999,1);
  int line_of_min = -1;
  bool go_on = true;
  
  // run through the rows of the matrix
  for (unsigned row = 0; row < s1; ++row)
  {
    line_of_min = -1;
    min = rational<T>(T(9999),T(1)); //min = rational<T>(9999,1);

    i = row;
    // begin: search for a row with the smallest non-zero entry at the lead position
    go_on = true;
    while (go_on)
    {
      if (output[i][lead] != this->Zero)
      {
        if (output[i][lead] < this->Zero) //if (output[i][lead] < 0)
        {
          for (j = 0; j < s2; ++j)
            output[i][j] *= -this->One;
        }
        if (output[i][lead] < min)
        {
          line_of_min = i;
          min = output[i][lead];
        }
      }
      if (++i >= s1)
      {
        if (line_of_min == -1)
        {
          i = row;
          ++lead;
          if(lead >= s2)
            return true;
        }
        else
          go_on = false;
      }
    }
    // end: search for a row with the smallest non-zero entry at the lead position

    if (line_of_min != row)
      swap(output[line_of_min], output[row]);

    rational<T> F1 = output[row][lead];

    // begin: use the entry of the current row to get zeros in the other rows
    const vector<rational<T> > &M_row = output[row];
    for (i = 0; i < s1; ++i)
    {
      if (i != row)
      {
        Entry = output[i][lead];
        if (Entry.denominator() < 0)
        {
          cout << "\n  Warning in template <class T> bool CLinAlg<T>::FindOrdersOfWilsonLines(...): the denominator is smaller zero. Return false." << endl;
          return false;
        }
        if (Entry != this->Zero)
        {
          for (j = 0; j < s2; ++j)
            output[i][j] = (F1 * output[i][j]) - (Entry * M_row[j]);
        }
      }
    }
    // end: use the entry of the current row to get zeros in the other rows
    ++lead;

    if (lead >= s2)
      return true;
  }
  return true;
}



/* ########################################################################################
######   bool CLinAlg<T>::ScaleToIntegerVector(...)                                  ######
######                                                                               ######
######   Version: 14.04.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
template <class T> bool CLinAlg<T>::ScaleToIntegerVector(vector<rational<T> > &Vector)
{
  unsigned i = 0;
  unsigned j = 0;

  T den = T(0);
  const size_t s1 = Vector.size();
  for (i = 0; i < s1; ++i)
  {
    if (Vector[i].denominator() != this->One)
    {
      den = Vector[i].denominator();
      for (j = 0; j < s1; ++j)
        Vector[j] *= den;
    }
  }
  return true;
}



/* ########################################################################################
######   bool CLinAlg::SimplifyInequality(...)                                       ######
######                                                                               ######
######   Version: 02.04.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) M : matrix M is simplified by adding rows such that as many entries      ######
######          as possible become positive                                          ######
######################################################################################## */
template <class T> bool CLinAlg<T>::SimplifyInequality(vector<vector<rational<T> > > &M)
{
  const size_t s1 = M.size();

  if (s1 < 2)
  {
    cout << "\n  Warning in bool CLinAlg::SimplifyInequality(...): Matrix is too small." << endl;
    return false;
  }

  const size_t s2 = M[0].size();

  unsigned i = 0;
  unsigned j = 0;

  rational<T> Entry = this->Zero; //rational<T> Entry(0,1);
  rational<T> tmp   = this->Zero; //rational<T> tmp(0,1);

  unsigned row = 0;
  unsigned lead = 0;
  bool positive_lead_not_found = true;

  // run through the columns
  while (lead < s2)
  {
    positive_lead_not_found = true;

    // begin: find the first row with a positive entry at the lead-position
    for (row = 0; positive_lead_not_found && (row < s1); ++row)
    {
      if (M[row][lead] > this->Zero)
        positive_lead_not_found = false;
    }
    --row;
    // end: find the first row with a positive entry at the lead-position

    // if a positive entry has been found
    if (!positive_lead_not_found)
    {
      // take the row that contains the positive entry at the lead-position
      const vector<rational<T> > &M_row = M[row];

      // run through the other rows of the matrix
      for (i = 0; i < s1; ++i)
      {
        // if the current row has a negative entry at the lead-position
        if ((i != row) && (M[i][lead] < this->Zero))
        {
          // make the negative entry positive by adding the row that contains the positive entry at the lead-position
          Entry = M[i][lead]/M[row][lead];

          for (j = 0; j < s2; ++j)
              M[i][j] -= Entry * M_row[j];
        }
      }
    }
    ++lead;
  }
  return true;
}



/* ##########################################################################
######   bool CLinAlg::SolveSystemOfLinearInequalities(...)            ######
######                                                                 ######
######   find solutions for Ax >= 0                                    ######
######   Version: 02.04.2008                                           ######
########################################################################## */
template <class T> void CLinAlg<T>::SolveSystemOfLinearInequalities(const vector<vector<rational<T> > > &A, vector<vector<rational<T> > > &Cone)
{
  /*
  const size_t m = A.size();
  const size_t n = A[0].size();

  vector<vector<rational<T> > > S = A;
  vector<vector<rational<T> > > U_current;
  vector<vector<rational<T> > > V_current;
  vector<vector<rational<T> > > S_current;
  vector<vector<rational<T> > > V1_current;
  vector<vector<rational<T> > > V2_current;

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned j2 = 0;
  unsigned l = 0;
  unsigned s = 0;
  rational<T> tmp1 = T(0);
  rational<T> tmp2 = T(0);

  rational<T> lvs = T(0);
  rational<T> lvk = T(0);

  // begin: U_current contains standard basis of R^n
  vector<rational<T> > empty(n,T(0));
  for (i = 0; i < n; ++i)
  {
    vector<rational<T> > tmp = empty;
    tmp[i] = T(1);
    U_current.push_back(tmp);
  }
  // end: U_current contains standard basis of R^n

  size_t u1 = 0;
  size_t v1 = 0;
  size_t s1 = 0;
  size_t s2 = 0;
  bool stop = true;

  i = 1;
  while(true)
  {
    cout << "i=" << i << ": U_current\n";
    size_t u1 = U_current.size();
    for (unsigned i1 = 0; i1 < u1; ++i1)
    {
      const vector<rational<T> > &u_i = U_current[i1];
      const size_t u2 = u_i.size();
      for (unsigned j1 = 0; j1 < u2; ++j1)
        cout << u_i[j1].numerator() << "/" << u_i[j1].denominator() << " ";
      cout << endl;
    }
    cout << endl;
    const vector<rational<T> > &S_i = S[i];

    u1 = U_current.size();

    // begin: find u from U_current with zero scalar product
    stop = false;
    vector<rational<T> > u;

    for (j = 0; !stop && (j < u1); ++j)
    {
      const vector<rational<T> > &u_j = U_current[j];
      tmp1 = T(0);
      for (k = 0; k < n; ++k)
        tmp1 += S_i[k] * u_j[k];
      if (tmp1 != T(0))
      {
        u = u_j;
        stop = true;
      }
    }
    // end: find u from U_current with zero scalar product

    if (stop)
    {
      vector<vector<rational<T> > > new_U_current;
      for (j = 0; j < u1; ++j)
      {
        const vector<rational<T> > &u_j = U_current[j];

        vector<rational<T> > new_u(n,T(0));

        tmp1 = T(0);
        tmp2 = T(0);
        for (k = 0; k < n; ++k)
        {
          tmp1 += S_i[k] * u[k];
          tmp2 += S_i[k] * u_j[k];
        }
        for (k = 0; k < n; ++k)
          new_u[k] = tmp1 * u_j[k] - tmp2 * u[k];

        new_U_current.push_back(new_u);
      }
      U_current = new_U_current;

      vector<vector<rational<T> > > new_V_current;
      vector<rational<T> > n_vector = -tmp1 * u;
      for (k = 0; k < n; ++k)
        n_vector[k] = n_vector[k]/abs(tmp1);
      new_V_current.push_back(n_vector);

      v1 = V_current.size();
      for (j = 0; j < v1; ++j)
      {
        const vector<rational<T> > &v_j = V_current[j];

        vector<rational<T> > new_v(n,T(0));

        tmp1 = T(0);
        tmp2 = T(0);
        for (k = 0; k < n; ++k)
        {
          tmp1 += S_i[k] * n_vector[k];
          tmp2 += S_i[k] * v_j[k];
        }
        for (k = 0; k < n; ++k)
          new_v[k] = (n_vector[k] * tmp2) - (v_j[k] * tmp1);

        new_V_current.push_back(new_v);
      }
      V_current = new_V_current;
    }
    else
    {
      // begin: create V1_current
      V1_current.clear();
      v1 = V_current.size();
      for (j = 0; j < v1; ++j)
      {
        const vector<rational<T> > &v_j = V_current[j];

        tmp1 = T(0);
        for (k = 0; k < n; ++k)
          tmp1 += S_i[k] * v_j[k];

        if (tmp1 <= T(0))
          V1_current.push_back(v_j);
      }
      // end: create V1_current
      V2_current.clear();

      // choose v_k
      for (k = 0; k < v1; ++k)
      {
        const vector<rational<T> > &v_k = V_current[k];

        lvk = T(0);
        for (l = 0; l < n; ++l)
          lvk += S_i[l] * v_k[l];

        // with scalar product negative
        if (lvk < T(0))
        {
          // choose v_s
          for (s = 0; s < v1; ++s)
          {
            const vector<rational<T> > &v_s = V_current[s];

            lvs = T(0);
            for (l = 0; l < n; ++l)
              lvs += S_i[l] * v_s[l];

            // with scalar product positive
            if (lvs > T(0))
            {
              // begin: create S^*
              vector<vector<rational<T> > > S_star;

              s1 = S_current.size();
              for (j = 0; j < s1; ++j)
              {
                const vector<rational<T> > &S_j = S_current[j];

                tmp1 = T(0);
                tmp2 = T(0);
                for (l = 0; l < n; ++l)
                {
                  tmp1 += S_j[l] * v_k[l];
                  tmp2 += S_j[l] * v_s[l];
                }
                if ((tmp1 == T(0)) && (tmp2 == T(0)))
                  S_star.push_back(S_j);
              }
              // end: create S^*

              // if S^* not empty
              s2 = S_star.size();
              if (s2 != 0)
              {
                // run through V_current, but do not take v_k or v_s
                for (j = 0; j < v1; ++j)
                {
                  if ((j != k) && (j != s))
                  {
                    const vector<rational<T> > &v_j = V_current[j];

                    // run through S^*
                    for (j2 = 0; j2 < s2; ++j2)
                    {
                      const vector<rational<T> > &S_star_j2 = S_star[j2];

                      // compute scalar product S_star_j2 * v_j
                      tmp1 = T(0);
                      for (l = 0; l < n; ++l)
                        tmp1 += S_star_j2[l] * v_j[l];

                      // if non-zero
                      if (tmp1 != T(0))
                      {
                        // add new vector to V2_current
                        vector<rational<T> > new_v2_vector(n,T(0));
                        for (l = 0; l < n; ++l)
                          new_v2_vector[l] = -lvs * v_k[l] + lvk * v_s[l];

                        V2_current.push_back(new_v2_vector);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    V_current = V1_current;
    V_current.insert(V_current.end(), V2_current.begin(), V2_current.end());
    S_current.push_back(S_i);

    ++i;
    if (i == m)
      break;
  }
  cout << "U_current\n";
  u1 = U_current.size();
  for (i = 0; i < u1; ++i)
  {
    const vector<rational<T> > &u_i = U_current[i];
    const size_t u2 = u_i.size();
    for (j = 0; j < u2; ++j)
      cout << u_i[j].numerator() << "/" << u_i[j].denominator() << " ";
    cout << endl;
  }
  cout << endl;

  cout << "V_current\n";
  u1 = V_current.size();
  for (i = 0; i < u1; ++i)
  {
    const vector<rational<T> > &v_i = V_current[i];
    const size_t v2 = v_i.size();
    for (j = 0; j < v2; ++j)
      cout << v_i[j].numerator() << "/" << v_i[j].denominator() << " ";
    cout << endl;
  }
  cout << endl;
*/

  if (Cone.size() != 0)
  {
    cout << "\n  Warning in void CLinAlg::SolveSystemOfLinearInequalities(...): Cone was not empty. Now cleared." << endl;
    Cone.clear();
  }
  const size_t s1 = A.size();

  if (s1 == 0)
  {
    cout << "\n  Warning in void CLinAlg::SolveSystemOfLinearInequalities(...): Matrix A is empty (1)." << endl;
    return;
  }
  const size_t s2 = A[0].size();

  if (s2 == 0)
  {
    cout << "\n  Warning in void CLinAlg::SolveSystemOfLinearInequalities(...): Matrix A is empty (2)." << endl;
    return;
  }

  //cout << "SolveSystemOfLinearInequalities: A.size() = " << s1 << endl;
  //cout << "SolveSystemOfLinearInequalities: A[0].size() = " << s2 << endl;

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  unsigned m = 0;

  if (s2 == 1)
  {
    bool all_positive = true;
    bool all_negative = true;
    for (i = 0; i < s1; ++i)
    {
      if (all_positive && (A[0][i] < T(0)))
        all_positive = false;
      if (all_negative && (A[0][i] > T(0)))
        all_negative = false;
    }
    if (!all_positive && !all_negative)
      return;

    vector<rational<T> > ConeVector(1,T(0));
    if (all_positive)
      ConeVector[0] = T(1);
    if (all_negative)
      ConeVector[0] = T(-1);
    Cone.push_back(ConeVector);
    return;
  }
  vector<vector<rational<T> > > A2 = A;
  /*
  cout << "A: \n";
  for (unsigned i = 0; i < A2.size(); ++i)
  {
    for (unsigned k = 0; k < A2[i].size(); ++k)
      cout << A2[i][k] << " ";
    cout << endl;
  }
  */
  size_t c1 = 0;

  const vector<rational<T> >    s2Line(s2, T(0));
  vector<vector<rational<T> > > EmptyMatrix(s2, s2Line);

  vector<vector<rational<T> > > MatrixA_i;
  vector<vector<rational<T> > > BasisOfH_i;
  vector<vector<rational<T> > > MatrixA_ij;
  vector<vector<rational<T> > > BasisOfH_ij;

  vector<vector<rational<T> > > BasisOfIntersections;

  vector<vector<rational<T> > > AllBasisOfH_i;

  rational<T> sp(T(0)); //rational<T> sp(0);
  bool check_again = false;
  bool use_vector  = true;

  // run through the rows of A
  for (i = 0; i < s1; ++i)
  {
    const vector<rational<T> > &A_i = A2[i];

    // begin: solve A_i x = 0
    //        the space spanned by x is named H_i,
    //        then z = (x + epsilon A_i) is a solution of A_i z > 0 for epsilon > 0
    BasisOfH_i.clear();

    // begin: make A_i quadratic
    MatrixA_i = EmptyMatrix;
    for (j = 0; j < s2; ++j)
    {
      for (k = 0; k < s2; ++k)
        MatrixA_i[j][k] = A_i[j] * A_i[k];
    }
    // end: make A_i quadratic

    if (!this->FindKernel(MatrixA_i, BasisOfH_i, true))
    {
      cout << "\n  Warning in void CLinAlg::SolveSystemOfLinearInequalities(...): Kernel of A_i not found." << endl;
      return;
    }
    AllBasisOfH_i.insert(AllBasisOfH_i.end(), BasisOfH_i.begin(), BasisOfH_i.end());
    // end: solve A_i x = 0

    // begin: compute H_ij, the intersection of H_i and H_j for j < i
    for (j = 0; j < i; ++j)
    {
      const vector<rational<T> > &A_j = A2[j];

      // begin: check whether A_i and A_j are linear dependent
      bool linear_dep = false;
      if ((A_i == A_j) || (A_i == s2Line))
        linear_dep = true;
      else
      {
        // assume that they are linear dependent
        linear_dep = true;

        bool first_found = false;
        rational<T> first_A_iOverA_j = T(0);

        for (k = 0; linear_dep && (k < s2); ++k)
        {
          if (((A_i[k] == T(0)) && (A_j[k] != T(0))) || ((A_i[k] != T(0)) && (A_j[k] == T(0))))
            linear_dep = false;
          else
          if ((A_i[k] != T(0)) && (A_j[k] != T(0)))
          {
            if (!first_found)
            {
              first_A_iOverA_j = A_i[k]/A_j[k];
              first_found = true;
            }
            else
            {
              if (first_A_iOverA_j != A_i[k]/A_j[k])
                linear_dep = false;
            }
          }
        }
      }
      // end: check whether A_i and A_j are linear dependent

      if (!linear_dep)
      {
        BasisOfH_ij.clear();

        // begin: make A_ij quadratic
        MatrixA_ij = EmptyMatrix;
        for (k = 0; k < s2; ++k)
        {
          for (l = 0; l < s2; ++l)
            MatrixA_ij[k][l] = (A_i[k] * A_i[l]) + (A_j[k] * A_j[l]);
        }
        // end: make A_ij quadratic

        if (!this->FindKernel(MatrixA_ij, BasisOfH_ij, true))
        {
          cout << "\n  Warning in void CLinAlg::SolveSystemOfLinearInequalities(...): Kernel of A_ij not found." << endl;
          return;
        }
        c1 = BasisOfH_ij.size();

        // begin: check whether the basis vectors of H_ij are compatible with A
        for (k = 0; k < c1; ++k)
        {
          vector<rational<T> > &Basisvector = BasisOfH_ij[k];

          check_again = false;
          for (l = 0; !check_again && (l < s1); ++l)
          {
            const vector<rational<T> > &A_l = A2[l];

            sp = T(0);
            for (m = 0; m < s2; ++m)
              sp += Basisvector[m] * A_l[m];

            if (sp < T(0))
            {
              check_again = true;
              for (m = 0; m < s2; ++m)
                Basisvector[m] *= T(-1); //Basisvector[m] *= (-1);
            }
          }
          use_vector = true;
          if (check_again)
          {
            for (l = 0; use_vector && (l < s1); ++l)
            {
              const vector<rational<T> > &A_l = A2[l];

              sp = T(0);
              for (m = 0; m < s2; ++m)
                sp += Basisvector[m] * A_l[m];

              if (sp < T(0))
                use_vector = false;
            }
          }
          if (use_vector)
          {
            // begin: if the vector is not known, add it
            bool new_vector = true;
            const size_t t1 = BasisOfIntersections.size();

            for (l = 0; new_vector && (l < t1); ++l)
            {
              if (BasisOfIntersections[l] == Basisvector)
                new_vector = false;
            }
            if (new_vector)
            {
              BasisOfIntersections.push_back(Basisvector);
              //cout << "good basis vector: ";
              //for (unsigned q = 0; q < Basisvector.size(); ++q)
              //  cout << Basisvector[q] << " ";
              //cout << endl;
            }
            // end: if the vector is not known, add it
          }
        }
        // end: check whether the basis vectors of H_ij are compatible with A
      }
    }
    // end: compute H_ij, the intersection of H_i and H_j for j < i
  }

  //size_t q1 = BasisOfIntersections.size();

  bool orthogonal = true;
  for (i = 0; i < AllBasisOfH_i.size(); ++i)
  {
    vector<rational<T> > &BasisVectorH_i = AllBasisOfH_i[i];

    orthogonal = true;
    for (j = 0; orthogonal && (j < BasisOfIntersections.size()); ++j)
    {
      const vector<rational<T> > &B = BasisOfIntersections[j];

      sp = T(0);
      for (k = 0; k < s2; ++k)
        sp += BasisVectorH_i[k] * B[k];

      if (sp != T(0))
        orthogonal = false;
    }
    if (orthogonal)
    {
      check_again = false;
      for (l = 0; !check_again && (l < s1); ++l)
      {
        const vector<rational<T> > &A_l = A2[l];

        sp = T(0);
        for (m = 0; m < s2; ++m)
          sp += BasisVectorH_i[m] * A_l[m];

        if (sp < T(0))
        {
          check_again = true;
          for (m = 0; m < s2; ++m)
            BasisVectorH_i[m] *= T(-1); //BasisVectorH_i[m] *= (-1);
        }
      }
      use_vector = true;
      if (check_again)
      {
        for (l = 0; use_vector && (l < s1); ++l)
        {
          const vector<rational<T> > &A_l = A2[l];

          sp = T(0);
          for (m = 0; m < s2; ++m)
            sp += BasisVectorH_i[m] * A_l[m];

          if (sp < T(0))
            use_vector = false;
        }
      }
      if (use_vector)
      {
        /*cout << "orthogonal row added ";
        for (unsigned q = 0; q < BasisVectorH_i.size(); ++q)
          cout << BasisVectorH_i[q] << " ";
        cout << endl;*/
        BasisOfIntersections.push_back(BasisVectorH_i);
      }
    }
  }
  // end: add rows of A that are orthogonal to all previous solutions

  Cone = BasisOfIntersections;
}


template class CLinAlg<long long int>;
template class CLinAlg<int>;
template class CLinAlg<CHugeInt>;

