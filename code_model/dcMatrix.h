//
//  dcMatrix.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef dcMatrix_h
#define dcMatrix_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <cstdlib>


using namespace std;

class dcMatrix
{
public:
	unsigned int nbRows;	//Dimensions
	unsigned int nbCols;
    
	vector<double> val;	//Value for each element
    
	double & operator()(unsigned int,unsigned int);
	double & operator[](unsigned int);
    
    unsigned int getNbRows()     {return nbRows;}
    unsigned int getNbCols()     {return nbCols;}
    
    
    // Constructors
	
	//dcMatrix(){nbRows=0;nbCols=0;}
	dcMatrix(){}
	
	dcMatrix(unsigned int h,unsigned int l) 
	{
		nbRows=h; nbCols=l; 
		vector<double> tmp(l*h,0.0);
		val=tmp;
	}
	
	dcMatrix(unsigned int n) 
	{
		nbRows=n; nbCols=n;
		vector<double> tmp(n*n, 0.0);
		val=tmp;
	}
	
    dcMatrix(string pathFile);
	
	dcMatrix(vector<double> v);	// dcMatrix (n,1) from a single vector(n)
	
	
    
    void    resize(unsigned int n)  
	{
		nbRows=n;nbCols=n; 
		val.resize(n*n);
	}
    
	void    resize(unsigned int nRows, unsigned int nCols)  
	{
		nbRows=nRows; nbCols=nCols; 
		this->val.resize(nRows*nCols);
	}
    
	void    display();
        
	void    RandomInit();
    
    // Files operations
	void    FromFile(const char*);
    void    FromFile(string);
	void	FromFile_Rows(string fileName, unsigned int nrow);
	
    void    WriteToFile(string);
	void    WriteToFileCSV(string);
	void    WriteToFileCSV(string fileName, vector<string> header);
	
	
    // Operations on embedded vectors
    vector<double>  extractColumn(unsigned int j_col);
	vector<double>	extractRow(unsigned int i_row);
	
	void            addRowVector(vector<double> v);
	void            addRowVector(vector<unsigned long> v);
	
    void            addColVector(vector<double> v);
    
	void			removeRow(unsigned int i_row);	// removes row 'i_row' and resize dcMatrix
	void			removeCol(unsigned int j_col);	// removes column 'j_col' and resize dcMatrix

	
	
	// Extract the row #i of the matrix
	// "i" is calculated such that it is
	// the smallest element of column "j_col"
  	vector<double>	extractRow_cond_minElement(unsigned int j_col);
	
    // Operations on elements
    
    void    setAllValues(double value);
	void	setValueFromMatrix(dcMatrix M);
	
	double  sumAllElements();
	double  sumLine(unsigned int i);		// sum all elements of line #i
	double  sumColumn(unsigned int j);	// sum all elements of line #i
	
	// conditional sum 
	double	sumColumnIf(unsigned int colToSum, unsigned int colToTest,
						double lowerBound, double upperBound);

	// counts nb of elements which are lower<element<upper  
	unsigned int		countColumnIf(unsigned int colToTest,
						  double lowerBound, double upperBound);
    
    void    setColumnValues(unsigned int colNb, vector<double> v);
	void	setRowValues(unsigned int rowNb_start0, vector<double> v);
	void	setRowValues(unsigned int rowNb_start0, vector<unsigned long> v);
    
    dcMatrix  transpose();
    
    bool    isSymetric();

    double  determinant();
	
	dcMatrix	getMinor(unsigned int row, unsigned int col);
	
	dcMatrix	inverse();
    
    dcMatrix  Cholesky();
    
    double  getMinimumValue();
    double  getMaximumValue();
	
};

dcMatrix operator + (dcMatrix &A,dcMatrix &B);
dcMatrix operator - (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (double a,dcMatrix &A);

dcMatrix Id(unsigned int n);

dcMatrix power(dcMatrix A,unsigned int n);

dcMatrix cholesky(dcMatrix A);	//renvoi la dcMatrix triangul L tq://si A symetrique,carre //L*transpo(L)=A

double distance_Matrix(dcMatrix A, dcMatrix B, double power);	// Euclidian distance b/w two matrices

dcMatrix rowBind(dcMatrix A, dcMatrix B);

dcMatrix	transpo(dcMatrix A);        // FIX ME : to delete if not used elsewhere
double	Det(dcMatrix A);           // FIX ME : to delete if not used elsewhere
unsigned int		test_sym(dcMatrix A);       // FIX ME : to delete if not used elsewhere



#endif
