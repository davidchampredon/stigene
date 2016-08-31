//
//  dcTools.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 . All rights reserved.

// Update : 2013-07-18

#ifndef dcTools_h
#define dcTools_h

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <sys/time.h>

// My Libraries

#include "dcMatrix.h"
#include "RV.h"



using namespace std;



void stopif(bool condition, string error_msg,
			int error_code=1, const char ff[]=__FUNCTION__);






int factorial(int i);
long combination(int n, int k);

template <class T> T dc_max(T a, T b) {
	T result;
	result = (a>b)? a : b;
	return (result);
}

template <class T> T dc_min(T a, T b) {
	T result;
	result = (a<b)? a : b;
	return (result);
}

void	coutline(unsigned int n);

double uniforme ();     // Simulates a random variable Uniform law on [0;1]


// ===================================================
// ================ FILES MANIPULATION ===============
// ===================================================

unsigned int nbLinesFile(string pathFile);

void delete_out_files(string path_out);

void vectorFromFile(vector<long>& res, const char * theFileName);
void vectorFromFile(vector<double>& res, const char * theFileName);
void vectorFromFile(vector<int>& res, string theFileName);

void vectorFromCSVfile(vector<double>& res, const char * theFileName, int column);
void vectorFromCSVfile_string(vector<string>& res, const char * theFileName, int column);

void MatrixFromCSVfile(dcMatrix& M,string filename, unsigned long ncol);


template <class T> void vectorToFile(vector<T> v, string fileName)
{
	ofstream theFile(fileName.c_str());
    
    for (int i=0; i<v.size(); i++) {
        theFile << v[i] << endl;
    }
}


template <class T> void vectorToCSVFile_Row(vector<T> x, string fileName)
{
	ofstream theFile(fileName.c_str(), ios::app);
    unsigned long n = x.size();
    
    for (int i=0; i<n-1; i++) {
        theFile << x[i] << ",";
    }
    
    theFile << x[n-1] << endl;
}

//template <class T> T 
double getParameterFromFile(string paramName, string fileName);
string getParameterFromFile_string(string paramName, string fileName);


vector<string>	getFirstLineHeaders( istream& ins );
unsigned int	getNumberColumns(string fileName);

dcMatrix			distributionFromFile(string filename);


typedef vector <double> record_t;
typedef vector <record_t> data_t;

istream& operator >> ( istream& ins, record_t& record );
istream& operator >> ( istream& ins, data_t& data );

void write_headers_if_emptyFile(string filename, string headers);



// ===================================================
// ================ VECTOR OPERATIONS ================
// ===================================================





// -------------------------------


// Note: Templates MUST be declared in the header file (not in .cpp)

template <class T> void displayVector(vector<T> v)
{
	cout << endl<< "(size="<<v.size()<<")"<<endl<<"[";
	for (int i=0; i<v.size()-1; i++) 
	{
		cout << v[i] << "; ";
		if ((i+1)%10==0) cout<<endl;
	}
	cout << v[v.size()-1];
	cout<< "]" << endl;	
}


template <class T> void displayVector(vector< vector<T> > v)
{
	cout << endl<< "(size="<<v.size()<<")"<<endl<<"[";
	for (int i=0; i<v.size(); i++)
	{
		displayVector(v[i]);
	}
	
}


vector<double> multVector(double a, vector<double> x);


vector<double> addVector(vector<double> a, vector<double> b);
vector<double> substractVector(vector<double> a, vector<double> b);
vector<double> divideVector(vector<double> vNumerator, vector<double> vDenominator);

double distanceLvector(int power, vector<double> x,vector<double> y );
double distanceLvector(int power, vector<double> x,vector<double> y, int firstElements );
double distanceLvector(int power, vector<double> x,vector<double> y, vector<int> weights );

double powerLvector(int power, vector<double> x,vector<double> y );

double	normLvector(int power, vector<double> x);

double	maxElementVector(vector<double> x);
double	minElementVector(vector<double> x);
int		argminElementVector(vector<double> x);

vector<double> max_vector(vector<double> x, vector<double> y);


double standardDeviation(vector<double> x);
double standardDeviation(vector<unsigned long int> x);

vector<double> smoothVector(vector<unsigned long int> x, int lagSmooth);

vector<double> numericalDerivative(vector<double> x, vector<double> t);
vector<double> numericalSecondDerivative(vector<double> x, vector<double> t, int lagDeriv);

double findAvgMax(vector<unsigned long int> x, vector<double> t, int lagSmooth, double minSlopeAboutMax=0);
double findAvgMin(vector<unsigned long int> x, vector<double> t, int lagSmooth, double minSlopeAboutMin=0);



template <class T> double averageElements(vector<T> x) 
{
	/// Returns the mean (average) of all elements
	T a = 0;
	for (int i=0; i<x.size(); i++) a += x[i];
	return ((double)a)/((double)x.size());
}

template <class T> double meanElements(vector<T> x)
{
	/// Returns the mean (average) of all elements
	return averageElements(x);
}


template <class T> double varianceElements(vector<T> x)
{
	stopif(x.size()<2,"Cannot calculate variance on vector of size <2");
	
	double m = meanElements(x);
	
	T a = 0;
	for (int i=0; i<x.size(); i++)
	{
		a += (x[i]-m)*(x[i]-m);
	}
	return a/(x.size()-1);
}


template <class T> T sumElements(vector<T> x) 
{
	/// Returns the sum of all elements from a vector
	
	T s = 0;
	for (int i=0; i<x.size(); i++) 
	{
		s += x[i];
	}
	return s;
}

template <class T> T sumElements(vector< vector<T> > x)
{
	/// Returns the sum of all elements from a vector of vector
	/// (not necessarily a rectangular matrix, may have different column lengths)
	
	T s = 0;
	for(int j=0;j<x.size(); j++)
		for (int i=0; i<x[j].size(); i++)
		{
			s += x[j][i];
		}
	return s;
}


template <class T> T extractElementRandom(vector<T> x)
{
	// Extracts an element randomly
	int rndPos = uniformInt(0, (int)(x.size()-1));
	
	return x[rndPos];
}


template <class T> vector<T> deleteElement(vector<T> x, int positionElementToDelete) 
{
	// Rather use "erase()" method of std vector library
	
	vector<T> y;
	for (int i=0; i<x.size(); i++) 
	{
		if (i!=positionElementToDelete)
			y.push_back(x[i]);
	}
	return y;
}




template <class T> bool isElementPresent(vector<T> x, T elemValue)
{
	// CHECK IF A GIVEN VALUE IS AN ELEMENT OF A VECTOR
	// (if vector empty, returns "false" anyway)
	
	bool isPresent = false;
	
	if (x.size()>0) 
	{
		for (unsigned long i=0; i<x.size(); i++)
		{
			if (x[i]==elemValue)
			{
				isPresent=true;
				break;
			}
		}
	}
	return isPresent;
}


template <class T> unsigned long findIndexElement(vector<T> x, T elemValue) 
{
	/// Find the position of a value in a vector
	
	unsigned long s = 0;
	bool isPresent = false;
	
	for (unsigned long i=0; i<x.size(); i++) 
	{
		if (x[i]==elemValue)
		{
			s=i;
			isPresent=true;
			break;
		}
	}
	
	if (!isPresent) 
	{
		cout << endl << "ERROR [findIndexElement]: element <"<< elemValue;
		cout << "> not found in this vector:";
		displayVector(x);
		exit(1);
	}
	
	return s;
}

template <class T> vector<T> popElementValue(vector<T> x, T valueElementToDelete) 
{
	/// Removes the element of a vector that has a given value.
	
	vector<T> y=x;

	unsigned long ii = findIndexElement(x, valueElementToDelete);
	
	y.erase(y.begin()+ii);

	return y;
}


vector<string> trim(vector<string> x);



template <class T> vector<unsigned long> distribution(vector<T> data, vector<T> data_breaks)
{
	/// *** HISTOGRAM WITH LEFT CLOSED CONVENTION ***
	
	/// CALCULATES NUMBER OF ELEMENTS OF VECTOR 'data' THAT
	/// ARE BETWEEN THE BINS DEFINED BY 'data_breaks'
	
	/// 'data_breaks' ** MUST ** COVER THE WHOLE RANGE OF 'data' VALUES
	
	/// DISTRIBUTION RESULT HAS SIZE 'N-1' (N:size of 'data_breaks')
	
	
	unsigned long n = data_breaks.size();
	
	// Initialize counter for distribution
	vector<unsigned long> cnt(n-1, 0);
	
	// Convention: LEFT CLOSED
	// cnt[k]: # elements b/w
	// data_breaks[k] (INCLUDED) and data_breaks[k+1] (EXCLUDED)
	
	for (int i=0; i<data.size(); i++)
	{
		int k = 1;
		T currValue = data[i];
		
		// find the interval where the value falls in vector data_breaks
		while (data_breaks[k] <= currValue  &&  k<=n-1)
		{ k++;}
		
		string errmsg = "Data_breaks not broad enough. Value: "+to_string(currValue)+"  is not between [" + to_string(data_breaks[0]) + ";" + to_string(data_breaks[n-1]) + "]";
		stopif(k>n-1,errmsg);
		
		// once the interval is found, add this element to the count.
		cnt[k-1]++;
	}
	
	// consistency check
	unsigned long total = data.size();
	unsigned long total2 = (unsigned long)(sumElements(cnt));
	stopif(total != total2, "Consistency check FAILED!");
	
	return cnt;
}


template <class T> vector<double> distributionNormalized(vector<T> data, vector<T> data_breaks)
{
	/// CALCULATES PROPORTION OF ELEMENTS OF VECTOR 'data' THAT
	/// ARE BETWEEN THE BINS DEFINED BY 'data_breaks'
	
	
	vector<unsigned long> dis = distribution(data, data_breaks);
	unsigned long n = dis.size();
	
	double sum = (double)(sumElements(dis));
	
	vector<double> res(n);
	
	for (int i=0; i<n; i++)
	{
		res[i] = (double)(dis[i])/sum;
	}
	
	return res;
}




template <class T> vector<unsigned long> distribution2(vector<T> data, vector<T> data_breaks)
{
	// TO DO: REPLACE "distribution" by this function???
	
	// CALCULATES NUMBER OF ELEMENTS OF VECTOR 'data' THAT
	// ARE BETWEEN THE BINS DEFINED BY 'data_breaks'
	
	// 'data_breaks' DOESN't NEED TO COVER THE WHOLE RANGE OF 'data' VALUES
	
	
	// Make sure 'breaks' cover min-max range of data
	T themin = minElementVector(data);
	T themax = maxElementVector(data);
	
	if (themin<data_breaks[0])
	{
		//insert value that will cover smallest range
		typename std::vector<T>::iterator it;
		it = data_breaks.begin();
		it = data_breaks.insert ( it , themin-1 );
	}
	
	if (themax>data_breaks[data_breaks.size()-1])
	{
		//insert value that will cover largest range
		data_breaks.push_back(themax+1);
	}
	
	// DISTRIBUTION RESULT HAS SIZE 'N-1' (N:size of 'data_breaks')
	
	int n = data_breaks.size();
	
	// Initialize counter for distribution
	vector<unsigned long> cnt(n-1, 0);
	
	// Convention:
	// cnt[k]: # elements b/w data_breaks[k] and data_breaks[k+1]
	
	for (int i=0; i<data.size(); i++)
	{
		int k = 1;
		T currValue = data[i];
		
		// find the interval where the value falls in vector data_breaks
		while (data_breaks[k]<currValue && k<=n-1)
		{ k++;}
		
		if (k>n-1)
		{
			cout << endl << "ERROR [distribution2]: data_breaks not broad enough"<<endl;
			cout << "Value: "<<currValue<<" is not between ["<<data_breaks[0]<<";"<<data_breaks[n-1]<<"]"<<endl;
			exit(1);
		}
		// once the interval is found, add this element to the count.
		cnt[k-1]++;
	}
	
	// consistency check
	unsigned long total = data.size();
	unsigned long total2 = (unsigned long)(sumElements(cnt));
	
	if (total != total2)
	{
		cout << endl << "ERROR [distribution2]: consistency check FAILED!"<<endl;
		exit(1);
	}
	return cnt;
}



template <class T> vector<double> distributionNormalized2(vector<T> data, vector<T> data_breaks)
{
	// CALCULATES PROPORTION OF ELEMENTS OF VECTOR 'data' THAT
	// ARE BETWEEN THE BINS DEFINED BY 'data_breaks'
	
	
	vector<unsigned long> dis = distribution2(data, data_breaks);
	unsigned long n = dis.size();
	
	double sum = (double)(sumElements(dis));
	
	vector<double> res(n);
	
	for (int i=0; i<n; i++)
	{
		res[i] = (double)(dis[i])/sum;
	}
	
	return res;
}



vector<double> vector_seq(double start, double end, unsigned int size);
vector<double> vector_seq_by(double start, double end, double by);


void my_pause(double waitTime);


// ===================================================
// ================    SAMPLINGS      ================
// ===================================================

vector<double> partitionLinear(double a_min, double a_max, int nbPartitions);

vector<double> partitionLog(double a_min, double a_max, int nbPartitions);

dcMatrix	LatinHypercubeSampling(vector<double> Vmin, vector<double> Vmax,
							   int samplingNb);

vector<int> uniformIntVector(int seed, int size, int min, int max); // Vector of size "size", elements randomly valued with uniform distribution
vector<long> uniformIntVectorUnique(long size, long min, long max); // Vector of size "size", elements randomly valued with uniform distribution. Cannot have duplicate elements




// ===================================================
// ==================== CONVERSION ===================
// ===================================================

string int2string(int i);


// ===================================================
// ==================== Compensation =================
// ===================================================

vector<double> comp_adj(double a, double min_a, double max_a,
						double b, double min_b, double max_b,
						double c_new);



// ===================================================
// ================  MATHS FUNCTIONS   ===============
// ===================================================


double pseudo_beta(double x, double a, double b);
double pseudo_gamma(double x, double shape, double scale);



#endif
