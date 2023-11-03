
#include "cmassmatrix.h"

using std::endl;
using std::setw;


/* ########################################################################################
######   CMassMatrix::CMassMatrix()                                                  ######
######                                                                               ######
######   Version: 29.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CMassMatrix object. Creates an empty mass matrix. ######
######################################################################################## */
CMassMatrix::CMassMatrix()
{
  this->Label_Row    = "";
  this->Label_Column = "";
}



/* ########################################################################################
######   CMassMatrix::CMassMatrix(...)                                               ######
######                                                                               ######
######   Version: 04.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum        : the vev-config where to look for the fields and couplings######
######   2) Label_Row     : field label A                                            ######
######   3) Label_Column  : field label B, create mass matrix from A_i B_j M_ij      ######
######   4) AutoTranspose : transpose the matrix if too many columns                 ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CMassMatrix object. Creates the mass matrix of the fields  ######
######   "Label_Row" and "Label_Column" from "Vacuum". If "AutoTranspose" transpose  ######
######   the matrix automatically if it has too many columns.                        ######
######################################################################################## */
CMassMatrix::CMassMatrix(const SConfig &Vacuum, const string &Label_Row, const string &Label_Column, bool AutoTranspose)
{
  this->FieldIndices_Row.clear();
  this->FieldIndices_Column.clear();

  const vector<CField> &Fields = Vacuum.Fields;
  const size_t f1 = Fields.size();
  
  unsigned i = 0;
  unsigned j = 0;

  // begin: find the field indices of "Label_Row"- and "Label_Column"-states
  if (Label_Row == Label_Column)
  {
    for(i = 0; i < f1; ++i)
    {
      if (Fields[i].Labels[Vacuum.use_Labels] == Label_Row)
        this->FieldIndices_Row.push_back(i);
    }
    // begin: sort
    bool not_inserted = true;
    vector<unsigned> tmp;
    const size_t s1 = this->FieldIndices_Row.size();
    size_t s2 = 0;
    for (i = 0; i < s1; ++i)
    {
      s2 = tmp.size();
      not_inserted = true;
      unsigned newNumber = Fields[this->FieldIndices_Row[i]].Numbers[Vacuum.use_Labels];

      for (j = 0; not_inserted && (j < s2); ++j)
      {
        if (newNumber < Fields[tmp[j]].Numbers[Vacuum.use_Labels])
        {
          not_inserted = false;
          tmp.insert(tmp.begin()+j, this->FieldIndices_Row[i]);
        }
      }
      if (not_inserted)
        tmp.push_back(this->FieldIndices_Row[i]);
    }
    this->FieldIndices_Row = tmp;
    // end: sort

    this->FieldIndices_Column = this->FieldIndices_Row;
  }
  else
  {
    for(i = 0; i < f1; ++i)
    {
      if (Fields[i].Labels[Vacuum.use_Labels] == Label_Row)
        this->FieldIndices_Row.push_back(i);
      else
      if (Fields[i].Labels[Vacuum.use_Labels] == Label_Column)
        this->FieldIndices_Column.push_back(i);
    }

    // begin: sort rows
    bool not_inserted = true;
    vector<unsigned> tmp;
    size_t s1 = this->FieldIndices_Row.size();
    size_t s2 = 0;
    for (i = 0; i < s1; ++i)
    {
      s2 = tmp.size();
      not_inserted = true;
      unsigned newNumber = Fields[this->FieldIndices_Row[i]].Numbers[Vacuum.use_Labels];

      for (j = 0; not_inserted && (j < s2); ++j)
      {
        if (newNumber < Fields[tmp[j]].Numbers[Vacuum.use_Labels])
        {
          not_inserted = false;
          tmp.insert(tmp.begin()+j, this->FieldIndices_Row[i]);
        }
      }
      if (not_inserted)
        tmp.push_back(this->FieldIndices_Row[i]);
    }
    this->FieldIndices_Row = tmp;
    // end: sort rows

    // begin: sort columns
    tmp.clear();
    s1 = this->FieldIndices_Column.size();
    s2 = 0;
    for (i = 0; i < s1; ++i)
    {
      s2 = tmp.size();
      not_inserted = true;
      unsigned newNumber = Fields[this->FieldIndices_Column[i]].Numbers[Vacuum.use_Labels];

      for (j = 0; not_inserted && (j < s2); ++j)
      {
        if (newNumber < Fields[tmp[j]].Numbers[Vacuum.use_Labels])
        {
          not_inserted = false;
          tmp.insert(tmp.begin()+j, this->FieldIndices_Column[i]);
        }
      }
      if (not_inserted)
        tmp.push_back(this->FieldIndices_Column[i]);
    }
    this->FieldIndices_Column = tmp;
    // end: sort columns
  }
  // end: find the field indices of "Label_Row"- and "Label_Column"-states

  const size_t number_rows    = this->FieldIndices_Row.size();
  const size_t number_columns = this->FieldIndices_Column.size();

  if ((number_rows == 0) || (number_columns == 0))
  {
    this->FieldIndices_Row.clear();
    this->FieldIndices_Column.clear();
    this->Matrix.clear();
    //this->vev_Matrix.clear();

    return;
  }

  // if too many columns, transpose the matrix
  if (AutoTranspose && (number_rows < number_columns) && (number_columns >= 10))
  {
    this->Label_Row    = Label_Column;
    this->Label_Column = Label_Row;

    this->FieldIndices_Row.swap(this->FieldIndices_Column);
  }
  else
  {
    this->Label_Row    = Label_Row;
    this->Label_Column = Label_Column;
  }

  this->Update(Vacuum);
}



/* ########################################################################################
######   bool CMassMatrix::IsEmpty() const                                           ######
######                                                                               ######
######   Version: 05.08.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is the matrix empty?                                         ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this mass matrix is empty.                                   ######
######################################################################################## */
bool CMassMatrix::IsEmpty() const
{
  // s1 = number of rows
  // s2 = number of columns
  const size_t s1 = this->FieldIndices_Row.size();
  const size_t s2 = this->FieldIndices_Column.size();

  unsigned j = 0;

  for (unsigned i = 0; i < s1; ++i)
  {
    const vector<vector<YukawaCoupling> > &Matrix_row = this->Matrix[i];

    for (j = 0; j < s2; ++j)
    {
      if (Matrix_row[j].size() != 0)
        return false;
    }
  }

  return true;
}



/* ########################################################################################
######   bool CMassMatrix::IsRankMaximal(...) const                                  ######
######                                                                               ######
######   Version: 29.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CMassMatrix::IsRankMaximal(unsigned &row_rank, unsigned &column_rank) const
{
  /*
  // s1 = number of (blown-up) rows
  // s2 = number of (blown-up) columns
  const size_t s1 = this->vev_Matrix.size();

  if (s1 == 0)
  {
    row_rank    = 0;
    column_rank = 0;

    return false;
  }

  const size_t s2 = this->vev_Matrix[0].size();

  row_rank    = findBasis<double>(this->vev_Matrix).size();
  column_rank = findBasis<double>(transpose(this->vev_Matrix)).size();

  if (s1 == s2)
  {
    if ((row_rank < s1) || (column_rank < s1))
      return false;
  }
  else
  {
    if (s1 < s2)
    {
      if ((row_rank < s1) || (column_rank < s1))
        return false;
    }
    else
    {
      if ((row_rank < s2) || (column_rank < s2))
        return false;
    }
  }
*/
  return true;
}



/* ########################################################################################
######   void CMassMatrix::SetRandom_vev()                                           ######
######                                                                               ######
######   Version: 16.08.2006                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
void CMassMatrix::SetRandom_vev()
{
  /*const size_t v1 = this->vev_Matrix.size();
  if (v1 == 0)
    return;

  const size_t v2 = this->vev_Matrix[0].size();

  unsigned i = 0;
  unsigned j = 0;

  for (i = 0; i < v1; ++i)
  {
    vector<double> &line = this->vev_Matrix[i];

    for (j = 0; j < v2; ++j)
    {
      if (line[j] != 0)
        line[j] = -5.576453 * ((double)(rand() + 1))/((double)(RAND_MAX + 1));
    }
  }*/
}



/* ########################################################################################
######   Update(const SConfig &Vacuum)                                               ######
######                                                                               ######
######   Version: 04.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum : the vev-config where to look for the fields and couplings       ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Update this mass matrix using the fields and couplings from "Vacuum".       ######
######################################################################################## */
void CMassMatrix::Update(const SConfig &Vacuum)
{
  // count the number of "Label_Row" and "Label_Column"
  const size_t number_of_columns = FieldIndices_Column.size();
  const size_t number_of_rows    = FieldIndices_Row.size();

  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  size_t t1 = 0;
  size_t t2 = 0;
  
  // begin: clear the matrix
  this->Matrix.clear();
  this->Matrix.resize(number_of_rows);

  for (i = 0; i < number_of_rows; ++i)
    this->Matrix[i].resize(number_of_columns);
  // end: clear the matrix

  vector<string> Labels;
  string tmp_Label1 = "";
  string tmp_Label2 = "";

  vector<YukawaCoupling> Couplings;
  CAnalyseModel Analyse;

  // run through the rows
  for (i = 0; i < number_of_rows; ++i)
  {
    std::ostringstream os1;
    os1 << Vacuum.Fields[this->FieldIndices_Row[i]].Numbers[Vacuum.use_Labels];
    tmp_Label1 = this->Label_Row + "_" + os1.str();
    
    // run through the columns
    for (j = 0; j < number_of_columns; ++j)
    {
      std::ostringstream os2;
      os2 << Vacuum.Fields[this->FieldIndices_Column[j]].Numbers[Vacuum.use_Labels];
      tmp_Label2 = this->Label_Column + "_" + os2.str();
      
      Labels.clear();
      Labels.push_back(tmp_Label1);
      Labels.push_back(tmp_Label2);

      Couplings.clear();
      Analyse.FindCouplings(Vacuum, Labels, Couplings, true);

      t1 = Couplings.size();
      for (k = 0; k < t1; ++k)
      {
        vector<unsigned> &FieldIndices = Couplings[k].FieldIndices;

        t2 = FieldIndices.size();
        for (l = 0; l < t2; ++l)
        {
          if (Vacuum.Fields[FieldIndices[l]].Labels[Vacuum.use_Labels] == this->Label_Row)
          {
            FieldIndices.erase(FieldIndices.begin() + l);
            break;
          }
        }
        t2 = FieldIndices.size();
        for (l = 0; l < t2; ++l)
        {
          if (Vacuum.Fields[FieldIndices[l]].Labels[Vacuum.use_Labels] == this->Label_Column)
          {
            FieldIndices.erase(FieldIndices.begin() + l);
            break;
          }
        }
      }
      this->Matrix[i][j] = Couplings;
    }
  }
}



/* ########################################################################################
######   CMassMatrix::~CMassMatrix()                                                 ######
######                                                                               ######
######   Version: 13.6.2006                                                          ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CMassMatrix object.                                ######
######################################################################################## */
CMassMatrix::~CMassMatrix()
{
}
