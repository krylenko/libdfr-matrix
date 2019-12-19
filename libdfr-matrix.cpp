 /******************* simple matrix/vector library ************************/
/*                          daniel ford                                  */
/*                           april 2011                                  */
/*************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include "libdfr-matrix.h"

/************************************************************************/
/*																		*/
/*						MATRIX FUNCTIONS								*/
/*																		*/	
/*************************************************************************/

// default constructor
Matrix::Matrix()
    : rows_(1), cols_(1)		// set dimensions
{

  // allocate memory
  int k=0;
  element_ = new double* [rows_];
  for (k = 0; k < rows_; k++)
      element_[k] = new double [cols_];

  // initialize matrix to zero
  zeros();

}

// normal constructor
Matrix::Matrix(const int& rows, const int& cols)
	: rows_(rows), cols_(cols)		// set dimensions
{

  // allocate memory
  int k=0;
  element_ = new double* [rows];
  for (k = 0; k < rows; k++)
      element_[k] = new double [cols];

  // initialize matrix to zero
  zeros();
  
}

// range fill constructor
Matrix::Matrix(const double& start, const double& end, const double& interval)
	: cols_(1)
{

  double length = (end-start)/interval;
  int L = (int)length;
  
  if( (L < length) && (L!=1) )
	rows_ = L+2;
  if( L == 1)
    rows_ = L;
  else
	rows_ = L+1;
  
  // allocate memory
  element_ = new double* [rows_];
  for (int k = 0; k < rows_; k++)
      element_[k] = new double [cols_];

  // fill matrix using specified interval
  for(int m=0; m < rows_; m++)
	element_[m][0] = start+interval*m;
	
}

// copy constructor
Matrix::Matrix(const Matrix& in)
{

	// allocate memory
	int i,j,k;
	element_ = new double* [in.rows()];
	for (k = 0; k < in.rows(); k++)
		element_[k] = new double [in.cols()];

	// make a copy of main array 
	for (i=0; i<in.rows(); i++)
	{
	  for(j=0;j<in.cols();j++)
		element_[i][j] = in.element_[i][j];
	}	  
	
	rows_ = in.rows();
	cols_ = in.cols();
	
}

// overloaded assignment operator
Matrix& Matrix::operator=(const Matrix& a)
{

	// check for assignment to itself
	if(this != &a)
	{

		// delete old memory
		delete [] element_;
		// allocate new memory
		int i,j,k;
		element_ = new double* [a.rows()];
		for (k = 0; k < a.rows(); k++)
			element_[k] = new double [a.cols()];

		// make a copy of main array 
		for (i=0; i<a.rows(); i++)
		{
		  for(j=0;j<a.cols();j++)
			element_[i][j] = a.element_[i][j];
		}	  
		
		rows_ = a.rows();
		cols_ = a.cols();
	}
	else printf("failed");
	return *this;
	
}

///////////////////////////////////////////////////////////////////////////////
//                         MATRIX GENERATION                                 //
///////////////////////////////////////////////////////////////////////////////

/***************************** IDENTITY MATRIX **************************/

// create identity matrix (overwrites any existing contents)
void Matrix::I()
{

  if( rows_ != cols_ ){
    cout << "matrix not square: Matrix::I" << endl;
    exit(-1);
  }

  int i,j;
	  
  //fill matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
      if( i == j )
          element_[i][j] = 1;
      else element_[i][j] = 0;
      }
  }
}

/***************************** ZERO MATRIX **************************/

// fill matrix with zeros
void Matrix::zeros()
{
  int i,j;

  //fill matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
      element_[i][j] = 0;
  }
}

/***************************** FILLED MATRIX **************************/

// generate matrix filled with specified value
void Matrix::fill(double scalar)
{

  int i,j;

  //fill matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
		element_[i][j] = scalar;
  }

}

// populate diagonal of square matrix with specified value
void Matrix::fillDiag(double scalar)
{

  int i,j;

  // populate diagonal of matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++) {
          if (i == j) {
            element_[i][j] = scalar;
          }
      }
  }

}

/***************************** RANDOM MATRIX **************************/

// generate random matrix
void Matrix::random(double max)
{
  int i,j;
  
  // seed random number generator
  srand(time(NULL));

  //fill matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
      element_[i][j] = max*(double)rand()/RAND_MAX;
  }

}

///////////////////////////////////////////////////////////////////////////////
//                         MATRIX OPERATIONS                                 //
///////////////////////////////////////////////////////////////////////////////

/************************************** TRANSPOSE MATRIX ********************/

// transpose matrix (replaces matrix with its transpose)
Matrix Matrix::T() const
{

  Matrix trans(cols_,rows_);
  int i,j;
    
  // transpose matrix
  for (i=0; i<cols_; i++)
  {
      for(j=0;j<rows_;j++)
		trans[i][j] = element_[j][i];
  }
  
  return trans;
  
}

/********************************* INVERT MATRIX **************************/
//    NOTE: As of 10.10.2010, returns bad results when a diagonal term of
//    the U matrix is zero (b/c further progress requires dividing by
//    these terms)

//invert square matrix using LU decomposition
Matrix Matrix::inv()
{

  Matrix invers(*this);

  int size=0;
  if( rows_ == cols_)
	size = rows_;
  else{
    cout << "matrix not square: Matrix::inv" << endl;
    exit(-1);
  }

  int i,j,k,n;
  double temp=0;

  // set L=identity, U=input
  Matrix L(size,size);
  L.I();
  Matrix U(invers);

  // multiplication factor storage array
  Matrix Umult(size,size);

  // C to solve LZ=C
  Matrix C(size,size);
  C.I();

  // Z to solve LZ=C
  Matrix Z(size,size);
  
  // swap index
  Matrix swapidx(1,size);

  // compute U and L matrices
  // good non-swapping code for backup
  n=1;
  for(j=0;j<size;j++){
     for(i=j+1;i<size;i++){
        Umult[i][j] = U[i][j]/U[j][j];
                //update matrix using multiplication factors
                for(k=j;k<size;k++){
                    U[i][k] -= (U[j][k]*Umult[i][j]);
                }
                //create L from multiplication factors
                L[i][j]=Umult[i][j];
     }
  n++;
  }

  //solve for Z
  for(i=0;i<size;i++){          //iterate for each column
     for(j=0;j<size;j++){  		//iterate for each column entry
        temp=0;
        for(k=0;k<size;k++){   //add elements for each column entry
           if( k != i ){
           temp += Z[k][j]*L[i][k];
           }
        }
     Z[i][j] = C[i][j]-temp;
     }
  }


  //use Z to calculate columns of inverse
  for(i=0;i<size;i++){          //iterate for each column
     for(j=(size-1);j>=0;j--){  //iterate for each column entry
        temp=0;
        for(k=j;k<size;k++){   //add elements for each column entry
           if( k != j ){
           temp += U[j][k]*invers[k][i];
           }
        }
        invers[j][i] = (Z[j][i]-temp)/U[j][j];
     }
  }

  return invers;

}

// pseudoinverse of matrix of any shape
Matrix Matrix::pinv()
{

	Matrix in(*this);
	return ((in.T()*in).inv())*in.T();

}

// reshape matrix
Matrix Matrix::reshape(const int& rows, const int& cols)
{

  if(rows*cols != rows_*cols_)
  {
    cout << "size mismatch: Matrix::reshape" << endl;
    exit(-1);
  }

  Matrix reshp(rows,cols);
  int i,j,n,p;
  
  // reshape matrix
  n=0, p=0;
  for(i=0; i<rows_; i++)
  {

	for (j=0; j<cols_; j++)
	{
          reshp[n][p] = element_[i][j];
		  p++;

		  if( cols == 1 )
		  {
			p=0;
			n++;
		  }
		  
		  else if( (cols_ > cols) && (cols != 1) )
		  {
			if( p == cols )
			{
				p=0;
				
				if( rows_ > rows)
				{
					if( (i+1) % rows == 0 )
						n++;
				}
				else
				{
					if( (((rows/rows_)-1) - n) != 0 )
						n++;
				} 	
			}
		  }
		  else
		  {
			if( p == cols )
			{
				p=0;
				
				if( rows_ > rows)
				{
					if( (i+1) >= (rows_/rows)*(n+1) )
						n++;
				}
				else
				{
					if( (((rows/rows_)-1) - n) != 0 )
						n++;
				} 	
			}
		  } 
    }
  }
  
  return reshp;
  
}

/****************************** MATRIX DOT PRODUCT ************************/

Matrix Matrix::operator*(const Matrix& input2)
{

  Matrix temp(*this);
  Matrix dot(temp.rows(),input2.cols());

  int i,j,m;

  //exit if there's a size mismatch
  if( temp.cols() != input2.rows() )
  {
	  cout << "size mismatch: Matrix::*" << endl;
      exit(-1);  
  }

  else
  {
	  //take dot product of matrices
	  for (i=0;i<temp.rows();i++)
	  {
		  for(j=0;j<input2.cols();j++)
		  {
				for(m=0;m<temp.cols();m++)
				{
					dot[i][j] += temp[i][m]*input2[m][j];
				}
		  }
	  }
  }

  return dot;	

}

/****************************** SLICE MATRIX ************************/

Matrix Matrix::slice(	const int& rowStart, const int& rowEnd, 
						const int& colStart, const int& colEnd)
{

  Matrix temp(*this);
  Matrix slice((rowEnd-rowStart)+1,(colEnd-colStart)+1);
  int n=0, p=0;

  //extract slice
  for(int j=colStart;j<=colEnd;j++){
    n = 0;
    for (int i=rowStart; i<=rowEnd; i++){
      slice[n][p] = temp[i][j];
		  n++;
    }
    p++;
  }

  return slice;

}

/****************************** MATRIX * SCALAR ************************/

Matrix Matrix::operator*(const double& b)
{

  Matrix temp(*this);
  int i,j;

  //multiply matrices
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		temp[i][j] = element_[i][j]*b;
      }
  }
  
  return temp;
}

/****************************** ELEMENT-WISE MULTIPLY ************************/
Matrix Matrix::emult(const Matrix& b)
{
  
  if( (b.rows() != rows_) || (b.cols() != 1) || (cols_ != 1) ){
    cout << "size mismatch: Matrix::emult" << endl;
    exit(-1);
  }

  Matrix temp(*this);
  
  for(int i=0;i<rows_;i++)
    temp[i][0] *= b[i][0];
    
  return temp;
  
}

/****************************** MATRIX / SCALAR ************************/

Matrix Matrix::operator/(const double& b)
{

  if(b == 0){
    cout << "divide by zero: Matrix::/" << endl;
    exit(-1);
  }

  Matrix temp(*this);
  int i,j;

  //multiply matrices
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		temp[i][j] = element_[i][j]/b;
      }
  }
  
  return temp;
}

/****************************** MATRIX + MATRIX ************************/

Matrix Matrix::operator+(const Matrix& b)
{

  Matrix temp(*this);

  if( (temp.rows()!=b.rows()) || (temp.cols()!=b.cols()) ){
    cout << "size mismatch: Matrix::+" << endl;
    exit(-1);
  }

  int i,j;

  //add matrices
  for (i=0; i<temp.rows_; i++)
  {
      for(j=0;j<temp.cols_;j++)
	  {
		temp[i][j] = element_[i][j]+b[i][j];
      }
  }

  return temp;
}

/****************************** MATRIX - MATRIX ************************/

Matrix Matrix::operator-(const Matrix& b) const
{

  Matrix temp(*this);

  if( (temp.rows()!=b.rows()) || (temp.cols()!=b.cols()) ){
    cout << "size mismatch: Matrix::-" << endl;
    exit(-1);
  }

  int i,j;

  //add matrices
  for (i=0; i<temp.rows_; i++)
  {
      for(j=0;j<temp.cols_;j++)
	  {
		temp[i][j] = element_[i][j]-b[i][j];
      }
  }

  return temp;
}

/****************************** MATRIX + SCALAR ************************/

Matrix Matrix::operator+(const double& b)
{

  Matrix temp(*this);

  int i,j;

  //add matrices
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		temp[i][j] = element_[i][j]+b;
      }
  }

  return temp;
}

/****************************** MATRIX - SCALAR ************************/

Matrix Matrix::operator-(const double& b)
{

  Matrix temp(*this);

  int i,j;

  //add matrices
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		temp[i][j] = element_[i][j]-b;
      }
  }

  return temp;
}

/****************************** COMPUTE SUM ****************************/

double Matrix::sum()
{

  int i,j;
  double sum = 0.0;

  //sum matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		sum += element_[i][j];
      }
  }

  return sum;
}

/*********************** COMPUTE SUM ALONG AN AXIS **************************/

double Matrix::sum(const int& axis, const int& idx)
{

  double sum = 0.0;

  //sum along given row
  if( axis == ROWS){
    if( idx > rows_-1 ){
      cout << "index overflow: Matrix::sum" << endl;
      exit(-1);
    }
    for (int i=0; i<cols_; ++i)
      sum += element_[idx][i];
  }

  //sum along given column
  if( axis == COLS){
    if( idx > cols_-1 ){
      cout << "index overflow: Matrix::sum" << endl;
      exit(-1);
    }
    for (int i=0; i<rows_; ++i)
      sum += element_[i][idx];
  }  
  
  return sum;
}

/****************************** COMPUTE AVERAGE ****************************/

double Matrix::mean()
{

  Matrix temp(*this);  
  return temp.sum()/(rows_*cols_);
  
}

/*********************** COMPUTE AVERAGE ALONG AN AXIS ************************/

double Matrix::mean(const int& axis, const int& idx)
{

  double mean = 0;
  Matrix temp(*this);
  if( axis == ROWS ){
    mean = temp.sum(ROWS,idx)/temp.cols();
  }
  if( axis == COLS ){
    mean = temp.sum(COLS,idx)/temp.rows();
  }
  return mean;
  
}

/****************************** VECTOR MAGNITUDE ****************************/
// expects matrix to be column vector
double Matrix::mag(){

  if( cols_ != 1 ){
    cout << "not column vector: Matrix::mag" << endl;
    exit(-1);
  }
  
  double mag = 0;
  for(int i=0;i<rows_;++i)
    mag += element_[i][0]*element_[i][0];
  
  return sqrt(mag);

}

/****************************** COMPUTE STD DEV ****************************/

double Matrix::stdDev()
{

  Matrix temp(*this); 

  int i,j;
  double std_dev = 0.0;
  double mean = temp.mean();
  
  //sum matrix
  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		std_dev += (mean-element_[i][j])*(mean-element_[i][j]);
      }
  }

  return std_dev;
}

/******************* STD DEV OF AN AXIS ****************************/

double Matrix::stdDev(const int& axis, const int& idx)
{

  Matrix temp(*this); 

  double std_dev = 0;
  
  // compute along given row
  if( axis == ROWS ){
    double mean = temp.mean(ROWS,idx);
    for (int i=0; i<cols_; i++)
      std_dev += (mean-element_[idx][i])*(mean-element_[idx][i]);
  }
  // compute along given column
  if( axis == COLS ){
    double mean = temp.mean(COLS,idx);
    for (int i=0; i<rows_; i++)
      std_dev += (mean-element_[i][idx])*(mean-element_[i][idx]);
  }  
  
  return std_dev;
}

/****************************** COMPUTE DETERMINANT ****************************/

// only up to 2x2!

double Matrix::det()
{

  // make sure matrix is 2x2
  if( rows_ > 2 || cols_ > 2 )
  {
    cout << "wrong size: Matrix::det" << endl;
    exit(-1);
  } 

  // return single element if matrix is 1x1
  if( rows_ == 1 || cols_ == 1 )
    return element_[0][0];
  
  double deter = 0.0;
  
  //determinant = ad-cb
  deter = element_[0][0]*element_[1][1] -
		  element_[0][1]*element_[1][0];

  return deter;
}

/****************************** NORMALIZE MATRIX ************************/
Matrix Matrix::norm()
{

  Matrix temp(*this);
  double tmpsum = temp.sum();
  if( tmpsum == 0 ){
    cout << "divide by zero: Matrix::norm" << endl;
    exit(-1);
  }
  for(int i=0; i<temp.rows(); i++)
  {
      for(int j=0;j<temp.cols();j++)
      {
        temp[i][j] = temp[i][j]/tmpsum;
      }
  }

  return temp;
  
}

/****************************** NORMALIZE COLUMNS/ROWS ************************/
Matrix Matrix::norm(const int& axis)
{

  Matrix temp(*this);
  
  if( axis == COLS ){
  
    for(int j=0; j<temp.cols(); j++){
      double tmpsum = temp.sum(COLS,j);
      if( tmpsum == 0 ){
        cout << "divide by zero: Matrix::norm" << endl;
        exit(-1);
      }
      for(int i=0;i<temp.rows();i++){
          temp[i][j] = temp[i][j]/tmpsum;
      }
    }
  
  }
  
  if( axis == ROWS ){
  
    for(int i=0; i<temp.rows(); i++){
      double tmpsum = temp.sum(ROWS,i);
      if( tmpsum == 0 ){
        cout << "divide by zero: Matrix::norm" << endl;
        exit(-1);
      }
      for(int j=0;j<temp.rows();j++){
          temp[i][j] = temp[i][j]/tmpsum;
      }
    }
  
  }

  return temp;
  
}

/*************************** ELEMENTWISE NATURAL LOG **********************/
Matrix Matrix::ln(){

  Matrix temp(*this);
  for (int i=0; i<temp.rows(); i++)
  {
      for(int j=0;j<temp.cols();j++)
      {
        temp[i][j] = log(temp[i][j]);
      }
  }

  return temp;
  
}

/****************************** MATRIX EXPONENTIAL ************************/

// uses Taylor series expansion, based on Matlab expm routine
// typically requires at least 10 iterations for reasonable precision
Matrix Matrix::exp(const int& iterations)
{

	Matrix A(*this);
	Matrix E(A);
	E.zeros();
	Matrix F(E);
	F.I();
	
	for(int k=1;k<iterations;++k)
	{
		E = E + F;
		F = (A*F)/k;
	}

	return E;
}

/****************************** MATRIX POWER ************************/

// raise matrix to integer power - naive implementation
Matrix pow(const int& power){
// to be implemented
}

/****************************** COLUMN MAX ************************/

Matrix Matrix::max(){

  Matrix maxtmp(cols_,1);
  
  for(int j=0; j<cols_; j++)
  {
      double maxS = element_[0][j];
      for(int i=0;i<rows_;i++)
      {
        if( element_[i][j] > maxS )
          maxS = element_[i][j];
      }
      maxtmp[j][0] = maxS;
  }

  return maxtmp;  

}

/****************************** ROW MAX ************************/
/*
Matrix Matrix::max(){

  Matrix max(rows_,1);
  
  for(int j=0; j<rows_; j++)
  {
      double maxS = element_[0][j];
      for(int i=0;i<cols_;i++)
      {
        if( element_[j][i] > maxS )
          maxS = element_[j][i];
      }
      max[j][0] = maxS;
  }

  return max;  

}
*/

/****************************** SCALAR MAX ************************/
double Matrix::maxS(){

  double max = element_[0][0];

  for(int i=0; i<rows_; i++)
  {
      for(int j=0;j<cols_;j++)
      {
        if( element_[i][j] > max )
          max = element_[i][j];
      }
  }

  return max;  

}

/****************************** ARGMAX ************************/
Matrix Matrix::argmax(){

  Matrix maxtmp(cols_,1);
  
  for(int j=0; j<cols_; j++)
  {
      double maxS = element_[0][j];
      maxtmp[j][0] = 0;
      for(int i=0;i<rows_;i++)
      {
        if( element_[i][j] > maxS ){
          maxS = element_[i][j];
          maxtmp[j][0] = i;
        }
      }
  }

  return maxtmp;

}

/****************************** DUPLICATE VECTOR INTO MATRIX ************************/
Matrix Matrix::dupe(const int& num){

  Matrix temp(*this);
  Matrix final(temp);

  if( rows_ == 1 ){
    for(int i=0;i<num;i++)
      final = final.conc(temp,ROWS);
  }
  else if( cols_ == 1 ){
    for(int i=0;i<num;i++)
      final = final.conc(temp,COLS);  
  }
  else {
    cout << "not a vector: Matrix::dupe" << endl;
    exit(-1);
  }

  return final;
  
}

/****************************** CONCATENATE MATRICES ************************/

Matrix Matrix::conc(const Matrix& b, const int& axis)
{

  int i,j,newrows,newcols;
  
  Matrix temp(*this);

  if( axis == COLS )
  {
    newrows = rows_;
    newcols = cols_+b.cols();
	
    if(rows_!=b.rows())
    {
      cout << "rows mismatch: Matrix::conc" << endl;
    }
  }  
  
  if( axis == ROWS )
  {
    newrows = rows_+b.rows();
    newcols = cols_;
	
    if(cols_!=b.cols())
    {
      cout << "columns mismatch: Matrix::conc" << endl;
    }
  }
  
  Matrix concat(newrows,newcols);
  
  if( axis == ROWS )
  {

      //concatenate matrices along rows
      for (i=0; i<(newrows); i++)
	  {
          for(j=0;j<cols_;j++)
		  {
			if( i<rows_ )
				concat[i][j] = temp[i][j];
			if( i>=rows_)
				concat[i][j] = b[i-rows_][j];
          }
      }
  }
  
  if( axis == COLS )
  {

      //concatenate matrices along columns
      for (i=0; i<rows_; i++)
	  {
          for(j=0;j<(newcols);j++)
		  {
			if( j<cols_ )
				concat[i][j] = temp[i][j];
			if( j>=cols_ )
				concat[i][j] = b[i][j-cols_];
          }
      }
  }

  return concat;

}

/****************************** COLLAPSE MATRIX ************************/

Matrix Matrix::collapse(const int& axis)
{

	const int COLS = 0;
	const int ROWS = 1;
	
	Matrix cllpse(1,rows_*cols_);
	int i,j,m = 0;
	
	if( axis == ROWS )
	{
		for(i=0;i<rows_;++i)
		{
			for(j=0;j<cols_;++j)
			{
				cllpse[0][m] = element_[i][j];
				++m;
			}
		}
	}
	
	if( axis == COLS )
	{
		for(j=0;j<cols_;++j)
		{
			for(i=0;i<rows_;++i)
			{
				cllpse[0][m] = element_[i][j];
				++m;
			}
		}
	}
	
	return cllpse;
}

///////////////////////////////////////////////////////////////////////////////
//                       INPUT AND OUTPUT                                    //
///////////////////////////////////////////////////////////////////////////////

/**************************** FILL MATRIX MANUALLY *************************/

//fill matrix with keyboard input
void Matrix::fillMan()
{

  int i,j;

  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
      printf("Enter element %d,%d: ",i,j);
      scanf("%lf", &element_[i][j]);
      fflush(NULL);
      }
  printf("\n");
  }

}

/**************************** FILL MATRIX FROM FILE *************************/

//fill n x m matrix from file
void Matrix::fromFile(const char filename[])
{

  ifstream inp;
  inp.open(filename, ios::in);

  int i,j;

  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
         inp >> element_[i][j];
      }
  }

  inp.close();
  
}

/**************************** PRINT MATRIX TO SCREEN **********************/

//print n x m matrix to screen
void Matrix::print(const int& blanks, const int& precision) const
{

  cout << setiosflags(ios::fixed) << setprecision(precision) << setiosflags(ios::left);
  int i,j;

  for (i=0; i<rows_; i++)
  {
    for(j=0;j<cols_;j++)
	{
		cout << setw(precision*2+2) << element_[i][j];
    }
	if( blanks == NOBLANKS )
		cout << "\n";
	else cout << "\n\n";
  }
  cout << "\n";
}

/**************************** PRINT MATRIX TO FILE **********************/
//print matrix to file
void Matrix::printFile(const char filename[], const int& blanks, const int& precision) const
{

  ofstream outp;
  outp.open(filename, ios::out);
  
  outp << setiosflags(ios::fixed) << setprecision(precision) << setiosflags(ios::left);
  int i,j;

  for (i=0; i<rows_; i++)
  {
      for(j=0;j<cols_;j++)
	  {
		outp << setw(precision*2+2) << element_[i][j];
      }
  if( blanks == NOBLANKS )
	outp << "\n";
  else outp << "\n\n";
  }

  outp.close();
  
}

/************************************************************************/
/*																		*/
/*						TRANSFORM FRAME FUNCTIONS						*/
/*																		*/	
/*************************************************************************/

// constructor
TransformFrame::TransformFrame(	const double& theta, const double& alpha,
								const double& A, const double& D)
	: 	Matrix(4,4),					// call Matrix constructor explicitly
		frame_(4,4)						// init frame matrix
{
	I();								// make frame_ an identity matrix
	updateLink(theta, alpha, A, D);		// fill the rest of the entries
}

// update entries in composed transform matrix
void TransformFrame::updateLink(	const double& theta, const double& alpha,
									const double& A, const double& D)
{

	theta_ = theta; alpha_ = alpha; A_ = A; D_ = D;

	element_[0][0] = cos(theta_);
	element_[0][1] = -sin(theta_)*cos(alpha_);
	element_[0][2] = sin(theta_)*sin(alpha_);
	element_[0][3] = A_*cos(theta_);
	
	element_[1][0] = sin(theta_);
	element_[1][1] = cos(theta_)*cos(alpha_);
	element_[1][2] = -cos(theta_)*sin(alpha_);
	element_[1][3] = A_*sin(theta_);
	
	element_[2][1] = sin(alpha_);
	element_[2][2] = cos(alpha_);
	element_[2][3] = D_;
	
}

Matrix TransformFrame::fwdTrans()
{

	Matrix coords(1,3);
	coords[0][0] = element_[0][3]; 
	coords[0][1] = element_[1][3];
	coords[0][2] = element_[2][3];
	
	return coords;
	
}

TransformFrame TransformFrame::operator*(const TransformFrame& input2)
{

  TransformFrame temp(*this);
  TransformFrame mult;

  int i,j,m;

	  // take dot product of matrices	
	  // shouldn't have to redefine this, but works for now
	  for (i=0;i<temp.rows();i++)
	  {
		  for(j=0;j<input2.cols();j++)
		  {
				for(m=0;m<temp.cols();m++)
				{
					mult[i][j] += temp[i][m]*input2[m][j];
				}
		  }
	  } 
  
  return mult;	

}										

