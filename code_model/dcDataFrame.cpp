//
//  dcDataFrame.cpp
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-08-18.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//


#include "dcDataFrame.h"


bool dcDataFrame::is_empty(){
	return (_value.getNbRows()==0);
}


// ==== MANIPULATIONS ====


void dcDataFrame::addrow(string varname, vector<double> values)
{
	// Add a row to the data frame
	// Size of values MUST equal the nb of columns of dcMatrix _value
	
	if (values.size()!= _value.getNbCols() &&  _value.val.size()>0 )
	{
		cerr << "ERROR dcDataFrame [addrow]:";
		cerr << "'values' vector size and dcMatrix nb columns 'value' do not match"<<endl;
		exit(1);
	}
	
	// If first insertion in empty data frame
	// then create default headers
	if(_colname.size()==0){
		vector<string> tmp(0);
		for (int j=0;j<values.size(); j++)
			tmp.push_back("V"+int2string(j));
		_colname = tmp;
	}
	
	_rowname.push_back(varname);
	_value.addRowVector(values);
}

void dcDataFrame::addrow(string varname,double value)
{
	vector<double> tmp;
	tmp.push_back(value);
	addrow(varname,tmp);
}

void dcDataFrame::addcol(string colname, vector<double> values){
	
	if (values.size()!= _value.getNbRows() &&  _value.val.size()>0 ){
		cerr << "ERROR dcDataFrame [addcol]:";
		cerr << "'values' vector size and dcMatrix nb rows 'value' do not match"<<endl;
		exit(1);
	}
	
	// If first insertion in empty data frame
	// then create default headers
	if(_colname.size()==0) _colname.push_back("V1");

	_colname.push_back(colname);
	_value.addColVector(values);
}


// =======================
// ==== SET FUNCTIONS ====
// =======================





// =======================
// ==== RETRIEVE DATA ====
// =======================


double	dcDataFrame::getValue(string varname, unsigned int j)
{
	// Find row position of 'varname'
	unsigned int i;
	bool found = false;
	
	for (i=0; i<=_rowname.size() && !found; i++)
	{
		if (_rowname[i]==varname) found=true;
	}
	i--;
	
	if (i==_rowname.size())
	{
		cout << endl << " ERROR [dcDataFrame::getValue]: variable name '"<<varname<<"' does not exist!"<<endl;
		exit(1);
	}
	
	return _value(i,j);
}


double	dcDataFrame::getValue(string rowname, string colname)
{
	// Find column position of 'valuename'
	unsigned int j;
	bool found = false;
	
	for (j=0; j<=_colname.size() && !found; j++)
	{
		if (_colname[j]==colname) found=true;
	}
	j--;
	
	if (j==_colname.size())
	{
		cout << endl << " ERROR [dcDataFrame::getValue]: value name '"<<colname<<"' does not exist!"<<endl;
		exit(1);
	}
	
	return getValue(rowname, j);
}





// =======================
// ==== DISPLAY ====
// =======================


void dcDataFrame::display()
{
	unsigned int n = _rowname.size();
	unsigned int ncol = _value.getNbCols();
	
	cout << endl << "dcDataFrame Content:"<<endl;
	
	string sep = "\t\t|\t";
	
	if (ncol==0) cout << "EMPTY!"<<endl;
	
	// Print headers
	if (ncol>0)
	{
		cout << "_rowname" << sep;
		
		for (int i=0; i<ncol; i++)
			cout << _colname[i] << sep;
		
		cout<<endl;
	}
	
	if (ncol==1)
	{
		coutline(20);
		
		for (int i=0; i<n; i++)
		{
			cout << _rowname[i] << sep <<  _value[i] <<endl;
		}
		
		coutline(20);
	}
	
	
	if (ncol>1)
	{
		coutline(80);
		
		for (int i=0; i<n; i++)
		{
			cout << _rowname[i] << sep;
			
			for (int j=0;j<ncol; j++)
				cout <<  _value(i,j)<<"\t";
			
			cout << endl;
		}
		
		coutline(80);
	}
	
}




void dcDataFrame::saveToCSV(string filename, bool col_headers)
{
	/// Save a dcDataFrame into a CSV file
	/// Option to write column headers or not.
	
	unsigned int N = _rowname.size();
	unsigned int M = _colname.size();
	
	ofstream f(filename);
	
	if (col_headers)
	{
		f << "row_name,";
		for(unsigned int i=0;i<M;i++)
		{
			f<<_colname[i];
			if(i< M-1) f<<",";
		}
		f << endl;
	}
	
	for (unsigned int i=0;i<N;i++)
	{
		f << _rowname[i] << ",";
		for(unsigned int j=0;j<M;j++)
		{
			f<<_value(i,j);
			if(j< M-1) f<<",";
		}
		f << endl;
	}
	
}


// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////


dcDataFrame rbind(dcDataFrame x, dcDataFrame y)
{
	/// Binds 2 dcDataFrames along rows (same number of columns)
	
	int nx = x.get_colname().size();
	int ny = y.get_colname().size();
	
	string errmsg = "cannot row bind if column number different ("+int2string(nx)+" != "+int2string(ny) + ")";
	stopif( nx!=ny ,errmsg);
	
	dcDataFrame z;
	z=x;
	
	for(int i=0;i<y.get_rowname().size();i++)
		z.addrow(y.get_rowname()[i], y.get_value().extractRow(i));
	
	return z;
}


// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////

void thetest(dcDataFrame& d)
{
	vector<double> xx;
	xx.push_back(1.234);
	xx.push_back(5.678);
	
	d.addrow("kiki", xx);
}






