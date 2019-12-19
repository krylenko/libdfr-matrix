/******************* simple matrix/vector library ************************/
/*                          daniel ford                                  */
/*                           april 2011                                  */
/*************************************************************************/

#ifndef LIBDFR_MATRIX_H
#define LIBDFR_MATRIX_H

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

enum{NOBLANKS, BLANKS};	// adds/removes blank line in printing
enum{ROWS, COLS};		// for concatenation, etc.

// define base class for the matrix object
class Matrix
{
	public:
	
		// initialization
        Matrix();                               // default constructor, makes 1x1 matrix
        Matrix(const int& rows,				    // normal constructor
               const int& cols=1);              // defaults to column vector
        Matrix(const double& start,			    // range fill constructor
               const double& end,
               const double& interval);
        Matrix(const Matrix& in);				// deep copy constructor
        Matrix& operator=(const Matrix& a);     // overloaded assignment operator
		~Matrix(){delete [] element_;};			// destructor
		
		// dimension and matrix entry getters
		int rows() const {return rows_;}
		int cols() const {return cols_;}
		double* operator[](int i)const{return element_[i];}
		
		// matrix creation functions
		void zeros();							// fill matrix with zero
		void I();								// create identity matrix
												// (exits if matrix is not square)
		void fill(double scalar);				// fill matrix with scalar
        void fillDiag(const double scalar);     // set all diagonal elements to scalar
		void random(double max=1);				// fill matrix with random numbers
												// [0-max]
		void fillMan();							// enter matrix entries manually
		
		// matrix manipulation
		Matrix T() const;								// transpose matrix
		Matrix reshape(	const int& rows,		// reshape matrix
						const int& cols);	 
		Matrix inv();							// invert (exits if matrix not square)
		Matrix pinv();							// compute pseudoinverse
		Matrix slice(	const int& rowStart,	// slice matrix
						const int& rowEnd, 
						const int& colStart,
						const int& colEnd);
						
		Matrix conc(	const Matrix& b,		// concatenate matrices
						const int& axis);									

		Matrix collapse(const int& axis);		// collapse to single-row matrix
						
		// multiplication and division
		Matrix operator*(const Matrix& b);		// multiply two matrices (dot product)
		Matrix operator*(const double& b);		// multiply matrix by scalar
		Matrix operator/(const double& b);		// divide matrix by scalar
        Matrix emult(const Matrix& b);          // element-wise multiplication of two col vecs
    
		// addition and subtraction
		Matrix operator+(const Matrix& b);		// add a matrix
		Matrix operator+(const double& b);		// add a scalar
		Matrix operator-(const Matrix& b) const;		// subtract a matrix
		Matrix operator-(const double& b);		// subtract a scalar

		// statistics
        double sum();                                   // for whole matrix
		double mean();
		double stdDev();
        double sum(const int& axis, const int& idx);    // along a row or column
        double mean(const int& axis, const int& idx);
        double mag();                                   // magnitude of a vector
        double stdDev(const int& axis, const int& idx);
		
		// other functions
        double det();                                   // calculate determinant (only for 2x2!)
        Matrix norm();                                  // normalize all matrix entries
        Matrix norm(const int& axis=COLS);              // normalize columns or rows
        Matrix ln();                                    // elementwise natural log
        Matrix exp(const int& iterations=10);           // matrix exponential
        Matrix pow(const int& power);                   // naive matrix power function
    
        Matrix max();                 // returns column vector with max entry in each column
        double maxS();                // returns max value of entire matrix
        Matrix argmax();              // return index of max point(s)
		Matrix dupe(const int& num);  // duplicate a row or column vector as many times as indicated 
    
		// display/utility
		void fromFile(const char filename[]);		// load matrix from file
        void printFile(const char filename[],       // print to file
					   const int& blanks=NOBLANKS,
					   const int& precision=7) const;
        void print(const int& blanks=NOBLANKS,
				   const int& precision=7) const;	// print to screen
	
	protected:
		// data storage for the class
		double** element_;
		int rows_;
		int cols_;
	
};

// transform frame class
class TransformFrame: public Matrix
{
	private:

		// 4x4 rotation and translation matrix
		Matrix frame_;	
		
		// Denavit-Hartenberg link parameters
		double theta_;	// rotation about Z axis
		double alpha_;	// rotation about X axis
		double A_;		// link length along X axis
		double D_;		// link length along Z axis
		
	public:
	
		// constructor
		TransformFrame(	const double& theta=0.,	
						const double& alpha=0.,	
						const double& A=0.,		
						const double& D=0.);				

		// update frame parameters
		void updateLink( 	const double& theta,
							const double& alpha,
							const double& A, 
							const double& D	);		
							
		// return [ x y z ] coordinates of end point
		Matrix fwdTrans();
		
		// multiply two frames
		// right-hand multiply, so (Link1*Link2*Link3).fwdTrans()
		// will return the end point of Link3
		TransformFrame operator*(const TransformFrame& b);
	
};

#endif
