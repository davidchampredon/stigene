//
//  dcDataFrame.h
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-08-18.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#ifndef __LocalSTI__dcDataFrame__
#define __LocalSTI__dcDataFrame__

#include <iostream>
#include "dcTools.h"
#include "dcMatrix.h"




using namespace std;

class dcDataFrame
{

	/// Simple Data frame (similar to R)
	///
	/// For example:
	///				| _colname[0]	| _colname[1]
	///	--------------------------------------------
	/// _rowname[0]	|	1.23		|	4.56
	/// _rowname[1]	|	5.67		|	6.78
	/// _rowname[2]	|	8.79		|	3.84
	///	--------------------------------------------
	///
	
	vector<string>	_rowname;
	vector<string>	_colname;
	dcMatrix		_value;



public:
    
	
	// ======================
    // ==== CONSTRUCTORS ====
    // ======================
	
	
	dcDataFrame(){
		dcMatrix M(0,0);
		_value = M;
	}
	
	dcDataFrame(vector<double> v, string header){
		dcMatrix M(v);
		_value = M;
		_colname.push_back(header);
		
		for(unsigned long i=0; i<v.size(); i++) _rowname.push_back(to_string(i));
	}

	dcDataFrame(vector<unsigned long> v, string header){
		
		vector<double> x;
		for(int i=0; i<v.size(); i++) x.push_back((double)(v[i]));
		dcMatrix M(x);
		_value = M;
		_colname.push_back(header);
		
		for(unsigned long i=0; i<v.size(); i++) _rowname.push_back(to_string(i));
	}
	
	
	dcDataFrame(vector<string> rowname, dcMatrix M, vector<string> colName)
	{
		if(rowname.size()!=M.getNbRows())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_rowname' vector size and dcMatrix nb rows 'value' do not match"<<endl;
			exit(1);
		}
		
		if(colName.size()!=M.getNbCols())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_colname' (headers) vector size and dcMatrix cols 'value' do not match"<<endl;
			exit(1);
		}
		
		_rowname = rowname;
		_value = M;
		_colname = colName;
	}
	
	
	dcDataFrame(vector<string> rowname, dcMatrix M)
	{
		if(rowname.size()!=M.getNbRows())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_rowname' vector size and dcMatrix nb rows 'value' do not match"<<endl;
			exit(1);
		}
		
		_rowname = rowname;
		_value = M;
		
		vector<string> tmp(0);
		for (int j=0;j<_value.getNbCols(); j++)
		{
			tmp.push_back("C"+int2string(j));
		}

		_colname = tmp;
	}
	
	
	dcDataFrame(dcMatrix M)
	{
		/// Construct a dcDataFrame from a dcMatrix
		/// (row and column names are given defualt values)
		
		int ncol = M.getNbCols();
		int nrow = M.getNbRows();
		stopif(nrow==0 || ncol==0, "Cannot construct  a dcDataFrame from empty matrix");
		
		_value = M;
		
		// Default column names
		vector<string> tmp(0);
		for (int j=0;j<ncol; j++)
			tmp.push_back("C"+int2string(j));
		
		_colname = tmp;
		
		// Default row names
		tmp.clear();
		for (int j=0;j<nrow; j++)
			tmp.push_back("R"+int2string(j));
		
		_rowname = tmp;
	}

	
	dcDataFrame(string filename, bool headers)
	{
		//
		/// CONSTRUCT A DATA FRAME FROM A FILE
		/// FILE MAY HAVE HEADERS (=FIRST LINE DEFINING VALUES NAMES)
		/// 1ST COLUMN = VARIABLES NAMES
		
		
		ifstream fh(filename);
		
		// Retrieve first line
		vector<string> headersNames = getFirstLineHeaders(fh);
		unsigned long n = headersNames.size();
		
		if (headers)
		{
			// Assume first column is not relevant (because column of variable names)
			headersNames.erase(headersNames.begin());
			_colname = headersNames;
		}
		
		if (!headers)
		{
			_colname.clear();
			
			for (int i=0; i<n-1; i++)
			{
				_colname.push_back("V"+int2string(i));
			}
			
		}
		
		// Retrieve variable names = 1st column
		vector<string> tmp_varname;
		vectorFromCSVfile_string(tmp_varname, filename.c_str(), 1);
		
		// Assume first row is not relevant if headers
		if (headers) tmp_varname.erase(tmp_varname.begin());
		
		_rowname = tmp_varname;
		
		
		// Read the values
		
		// Read the matrix including the varnames and (potential) headers
		dcMatrix M;
		MatrixFromCSVfile(M, filename, n);
		
		// Remove 1st column = var names
		M.removeCol(0);
		
		// Remove 1st line (if headers)
		if (headers) M.removeRow(0);
		
		_value = M;
	}
	
	
	void clear(){
		_value.resize(0);
		_colname.clear();
		_rowname.clear();
	}
	
	bool is_empty();
	
	// ==== MANIPULATIONS ====
	
	void addrow(string varname, vector<double> values);
	void addrow(string varname, double value);
	void addcol(string colname, vector<double> values);
	
	// ==== SET FUNCTIONS ====
	
	void set_rowname(vector<string> x){_rowname=x;}
	void set_colname(vector<string> x){_colname=x;}
	
	
	// ==== RETRIEVE DATA ====
	
	vector<string>	get_rowname() {return _rowname;}
	vector<string>	get_colname() {return _colname;}
	dcMatrix		get_value() {return _value;}
	
	
	double	getValue(string varname, string valuename);
	double	getValue(string varname, unsigned int j);
	double	getValue(string varname);
	
	unsigned int get_nrows() {return _value.getNbRows();}
	unsigned int get_ncols() {return _value.getNbCols();}
	
	
	// ==== FILE MANIPULATION ====
	void	saveToCSV(string filename,bool headers);
	
	// ==== DISPLAY ====

	void display();
	
	
};

// =================================
// ===== OUTSIDE CLASS ======
// =================================

dcDataFrame rbind(dcDataFrame x, dcDataFrame y);

void thetest(dcDataFrame& d);



//void to_paste_in_main_cpp()
//{
//	dcDataFrame D;
//	
//	D.addrow("test", 0.2);
//	
//	//D.display();
//	
//	dcDataFrame DD;
//	vector<double> xx,yy;
//	xx.push_back(1.234);
//	xx.push_back(5.678);
//	
//	std::vector<double>::iterator it;
//	it = xx.begin();
//	it = xx.insert ( it , 99.9 );
//	displayVector(xx);
//
//}
//
#endif /* defined(__LocalSTI__dcDataFrame__) */
