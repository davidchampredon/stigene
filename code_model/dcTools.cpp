//
//  dcTools.cpp
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "dcTools.h"



void stopif(bool condition, string error_msg,
			int error_code, const char ff[])
{
	if (condition)
	{
		cerr << endl << " *=*=*=*=*=*=* ERROR *=*=*=*=*=*=* " << endl<<endl;
		cerr<<"In function: "<<string(ff)<<endl<<endl;
		cerr << error_msg <<endl;
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* " << endl;
		exit(error_code);
	}
}

int factorial(int i) 
{
    if (i <= 1)
        return i;
    return (i * factorial(i - 1));
}


long combination(int n, int k)
{
	assert(n>=k);
	return factorial(n)/factorial(k)/factorial(n-k);
}


double max(double a, double b)
{
	double c_new=b;
	if (a>b) c_new=a;
	
	return c_new;
}


double uniforme ()
{ //simulation d'une VA uniforme sur [0;1]
	
	return 1.0*random()/RAND_MAX;
}

void coutline(unsigned int n)
 {cout<<endl; for (int i=0;i<n;i++) cout<<"_"; cout<<endl;}


// ===================================================
// ================ FILES MANIPULATION ===============
// ===================================================

unsigned int nbLinesFile(string pathFile)
{
	ifstream f(pathFile.c_str());
	unsigned int n=0;
	string tmp;
	
	while (f.good()) 
	{
		n++;
		getline(f, tmp);
	}
	
	return n;
}


void delete_out_files(string path_out)
{
	string cmdline = "rm -f " + path_out + "/*.out";
	system(cmdline.c_str());
}


void vectorFromFile(vector<long>& res, const char * theFileName)
{
	
	ifstream thefile (theFileName); // declare file stream
	
	//assert(thefile.is_open());
	
	if (!thefile)
	{
		cout<<endl<<" ERROR [vectorFromFile]: This file is not found: "<<theFileName<<endl;
		exit(1);
	}
	
	vector<long> x;
	
	long line;
	
	x.clear();
	
	while (!thefile.eof())   
	{
        
		thefile >> line;
		if (thefile.eof()) break;
		
		x.push_back (line);	
		
	}
	
	thefile.close();
	res = x;
}

void vectorFromFile(vector<double>& res, const char * theFileName)
{
	ifstream thefile (theFileName); // declare file stream
	
	if (!thefile)
	{
		cout<<endl<<" ERROR [vectorFromFile]: This file is not found: "<<theFileName<<endl;
		exit(1);
	}
	
	
	vector<double> x;	
	double line;
	
	x.clear();
	
	while (!thefile.eof())   
	{
		thefile >> line;
		if (thefile.eof()) break;
		x.push_back (line);	
		//cout<<"DEBUG adding : " << line<<endl;
	}
	
	thefile.close();
	res = x;
}

void vectorFromFile(vector<int>& res, string theFileName)
{
	ifstream thefile (theFileName.c_str()); // declare file stream
	
	if (!thefile)
	{
		cout<<endl<<" ERROR [vectorFromFile]: This file is not found: "<<theFileName<<endl;
		exit(1);
	}
	
	
	vector<int> x;	
	int line;
	
	x.clear();
	
	while (!thefile.eof())   
	{
		thefile >> line;
		if (thefile.eof()) break;
		x.push_back (line);	
		//cout<<"DEBUG adding : " << line<<endl;
	}
	
	thefile.close();
	res = x;
}






typedef vector <double> record_t;
typedef vector <record_t> data_t;

typedef vector <string> record_t_string;
typedef vector <record_t_string> data_t_string;


//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.

istream& operator >> ( istream& ins, record_t& record )
{
	
	record.clear();	// make sure that the returned record contains only the stuff we read now
	
	string line;
	getline( ins, line );	// read the entire line into a string (a CSV record is terminated by a newline)
    
	stringstream ss( line );		// now we'll use a stringstream to separate the fields out of the line
	string field;
    
    //double tmp; int n;
    
	while (getline( ss, field, ',' ))
    {
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs( field );
		double f = 0.0;  // (default value is 0.0)
		fs >> f;
		
		record.push_back( f );		// add the newly-converted field to the end of the record
        
		//n=record.size();       // FIXME: to delete
		//tmp = record[n-1];
    }
	
	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	
	return ins;
}

istream& operator >> ( istream& ins, record_t_string& record )
{
	
	record.clear();	// make sure that the returned record contains only the stuff we read now
	
	string line;
	getline( ins, line );	// read the entire line into a string (a CSV record is terminated by a newline)
    
	stringstream ss( line );		// now we'll use a stringstream to separate the fields out of the line
	string field;
    
	while (getline( ss, field, ',' ))
    {
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs( field );
		string f = "";  // (default value is 0.0)
		fs >> f;
		
		record.push_back( f );		// add the newly-converted field to the end of the record
        
		//n=record.size();       // FIXME: to delete
		//tmp = record[n-1];
    }
	
	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	
	return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.

istream& operator >> ( istream& ins, data_t& data )
{
	data.clear();
	
	// For every record we can read from the file, append it to our resulting data
	record_t record;
	while (ins >> record)
    {
		data.push_back( record );
    }
	
	return ins;  
}

istream& operator >> ( istream& ins, data_t_string& data )
{
	data.clear();
	
	// For every record we can read from the file, append it to our resulting data
	record_t_string record;
	while (ins >> record)
    {
		data.push_back( record );
    }
	
	return ins;
}




void vectorFromCSVfile(vector<double>& res, const char * theFileName, int column)
{
	
	/// Read the nth column of a CSV file
	/// Expects 'double' data
	
	data_t data;
	
	// Here is the file containing the data. Read it into data.
	ifstream infile( theFileName );
	
	string errmsg = "This file is not found: " + string(theFileName);
	stopif(!infile, errmsg);

	infile >> data;
	
	infile.close();
	
	vector<double> x(data.size());
	
	for (unsigned i = 0; i < data.size(); i++)
		x[i] = data[i][column-1];
	
	res = x;
}



void vectorFromCSVfile_string(vector<string>& res, const char * theFileName, int column)
{
	/// Read the nth column of a CSV file
	/// Expects 'string' data
	
	data_t_string data;
	
	// Here is the file containing the data. Read it into data.
	ifstream infile( theFileName );
	
	if (!infile)
	{
		cout<<endl<<" ERROR [vectorFromCSVfile]: This file is not found: "<<theFileName<<endl;
		exit(1);
	}
	
	infile >> data;
	infile.close();
	
	vector<string> x(data.size());
	
	for (unsigned i = 0; i < data.size(); i++)
		x[i] = data[i][column-1];
	
	res = x;
}

void MatrixFromCSVfile(dcMatrix& M,string filename, unsigned long ncol)
{
	/// CONSTRUCT A MATRIX BY READING
	/// COLUMN VECTORS IN A CSV FILE
	/// MUST KNOW IN ADVANCE THE NuMBER OF COLUMNS TO BE READ

	dcMatrix A;
	
	for (int j=0;j<ncol;j++)
	{
		vector<double> vectcol;
		vectorFromCSVfile(vectcol, filename.c_str(), j+1);

		if (j==0) A.resize(vectcol.size(), ncol);
		
		A.setColumnValues(j, vectcol);
	}
	
	M = A;
}


dcMatrix distributionFromFile(string filename)
{
	/// RETRIEVE DISTRIBUTION
	/// FORMAT EXPECTED: CSV
	/// 1ST COLUMN: BREAKS
	/// 2ND COLUMN: PROPORTIONS
	/// PROPORTION[i] for BREAK[i] <= x < BREAK[i+1]
	
	vector<double> breaks;
	vector<double> proportion;
	
	vectorFromCSVfile(breaks, filename.c_str(), 1);
	vectorFromCSVfile(proportion, filename.c_str(), 2);
		
	dcMatrix res(breaks);
	res.addColVector(proportion);
	
	return res;
}


double getParameterFromFile(string paramName, string fileName)
{
	// THE FIRST COLUMN MUST HAVE THE NAME OF THE PARAMETER
	// THE 2ND COLUMN ITS VALUE
	// COLUMNS ARE SEPARATED BY A COMMA (",") 
	// WARNING: NO SPACE ALLOWED IN THE NAME (especially between name and coma)
	
	ifstream theFile(fileName.c_str());
	
	string errmsg = "Can't open file: "+fileName;
	stopif(!theFile,errmsg);
	
	
	bool found = false;
	
	vector<string> fileContent;
	string sValue;
	
	while ( theFile.good() ) 
	{
		string line;
		getline(theFile,line);
		
		if (line.find(paramName, 0) != string::npos)
		{
			// Enter here it at least the name of 
			// the parameter is in the string
			// (but not done yet! there may be variations on the name
			// for example, "rate_1", ""rate_2",etc)
			
			// Find the position of the comma
			// (delimits the name from its value)
			int pos = line.find(",");
			
			// Actual name read on this current line
			string currentName = line.substr(0, pos);
			
			// Check if the current name is the exact match
			if (currentName.compare(paramName)==0)
			{
				found = true;			
				sValue = line.substr(pos+1, line.length()-pos);
			}
		}
	}
	
	// Value has type "double"
	double res = NAN;
	
	if (found) res=atof(sValue.c_str());
	
	errmsg.clear();
	errmsg = "Parameter name '" + paramName + "' not found in file '" + fileName+"'";
	stopif(!found,errmsg);
	

	return res;
}


string getParameterFromFile_string(string paramName, string fileName)
{
	/// READ PARAMETER FROM FILE.
    
    /// THE FIRST COLUMN MUST HAVE THE NAME OF THE PARAMETER
	/// THE 2ND COLUMN ITS VALUE
	/// COLUMNS ARE SEPARATED BY A COMMA (",")
	/// WARNING: NO SPACE ALLOWED IN THE NAME (especially between name and coma)
	
	ifstream theFile(fileName.c_str());
	
	
	string errmsg = "Can't open file: "+fileName;
	stopif(!theFile,errmsg);
	bool found = false;
	
	vector<string> fileContent;
	string sValue;
	
	while ( theFile.good() )
	{
		string line;
		getline(theFile,line);
		
		if (line.find(paramName, 0) != string::npos)
		{
			// Enter here it at least the name of
			// the parameter is in the string
			// (but not done yet! there may be variations on the name
			// for example, "rate_1", ""rate_2",etc)
			
			// Find the position of the comma
			// (delimits the name from its value)
			int pos = line.find(",");
			
			// Actual name read on this current line
			string currentName = line.substr(0, pos);
			
			// Check if the current name is the exact match
			if (currentName.compare(paramName)==0)
			{
				found = true;
				sValue = line.substr(pos+1, line.length()-pos);
			}
		}
	}
	
	// Clean up if carriage return present
	if (!sValue.empty() && sValue[sValue.size() - 1] == '\r')
		sValue.erase(sValue.size() - 1);
	
	// Value has type "string"
	string res;
	
	if (found) res=sValue;
	
	errmsg.clear();
	errmsg = "Parameter name '" + paramName + "' not found in file '" + fileName+"'";
	stopif(!found,errmsg);

	return res;
}





vector<string> getFirstLineHeaders( istream& ins )
{
	/// READ THE FIRST LINE (HEADERS)
	/// OF A CSV FILE
	
	vector<string> res;
	
	string line;
	getline( ins, line );	// read the entire line into a string (a CSV record is terminated by a newline)
	
	stringstream ss( line );		// now we'll use a stringstream to separate the fields out of the line
	string field;
	
	while (getline( ss, field, ',' ))
	{
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs( field );
		string f = "";
		fs >> f;
		
		res.push_back( f );		// add the newly-converted field to the end of the record
		
		
	}
	
	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	
	return res;
}


unsigned int getNumberColumns(string filename)
{
	/// Returns the number of columns of a CSV file
	/// (based on the first line of this file)
	
	ifstream f(filename);
	
	vector<string> x = getFirstLineHeaders(f);
	return x.size();
}


void write_headers_if_emptyFile(string filename, string headers)
{
	/// Write headers if file is empty
	
	ifstream fcheck(filename.c_str());
	string tmp;
	getline(fcheck,tmp);
	fcheck.close();
	bool isempty =(tmp.length()==0);
	ofstream ff(filename.c_str(),ios::app);
	if (isempty)ff << headers <<endl;
}


// ===================================================
// ================ VECTOR OPERATIONS ================
// ===================================================


vector<string> trim(vector<string> x)
{
	/// Remove elements that are empty strings
	
	vector<string> res;
	
	for(int i=0;i<x.size();i++)
	{
		if(x[i]!="") res.push_back(x[i]);
	}
	return res;
}


vector<double> multVector(double a, vector<double> x)
{	
    vector<double> res(x.size());
    
	for (int i=0; i<x.size(); i++) 
    {
		res[i] = a * x[i];
	}
    
    return res;
}


vector<double> addVector(vector<double> a, vector<double> b)
{
	assert(a.size()==b.size());
	
	vector<double> res(a.size());
	
	for (int i=0; i<a.size(); i++) 
	{
		res[i] = a[i]+b[i];
	}
	
	return res;
}

vector<double> substractVector(vector<double> a, vector<double> b)
{
	assert(a.size()==b.size());
	
	vector<double> res(a.size());
	
	for (int i=0; i<a.size(); i++) 
	{
		res[i] = a[i]-b[i];
	}
	
	return res;
}

vector<double> divideVector(vector<double> vNumerator, vector<double> vDenominator)
{
	assert(vNumerator.size()==vDenominator.size());
    
    bool isZero = false;
    
    for (int c_new=0; c_new<vDenominator.size(); c_new++) {
        if (fabs(vDenominator[c_new])<1e-15) isZero=true;

    }
    if(isZero)   // Cannot divide by an element=0
    {
        cout << endl << "# # # ERROR: CANNOT DIVIDE WITH VECTOR WHOSE AT LEAST ONE ELEMENT IS 0"<<endl;
        displayVector(vDenominator);
        exit(1);
    }
	
	vector<double> res(vNumerator.size());
	
	for (int i=0; i<vNumerator.size(); i++) 
	{
		res[i] = vNumerator[i]/vDenominator[i];
	}
	
	return res;
}

double distanceLvector(int power, vector<double> x,vector<double> y )
{
	// Calculates the L-n distance between 2 vectors
	
	
	int n=x.size();
	string errmsg = "Vectors do not have the same size: "+to_string(n)+"!="+to_string(y.size());
	stopif(n != y.size(),errmsg);
	
	double res=0.0;
	
	for(int i=0; i<n; i++)
	{
		res += pow(x[i]-y[i],power);
	}
	return pow(res,1.0/power);
}

double distanceLvector(int power, vector<double> x,vector<double> y, int firstElements )
{
	// Calculates the L-n distance between 2 vectors
	// TAKES INTO ACCOUNT THE "N" FIRST ELEMENTS
	
	int n=x.size();
	assert(n == y.size());
	assert(firstElements<n+1);
	
	double res=0.0; 
	
	for(int i=0; i<firstElements; i++)
	{
		res += pow(x[i]-y[i],power);
	}
	
	return pow(res,1.0/power);
}


double distanceLvector(int power, vector<double> x,vector<double> y, vector<int> weights )
{
	// Calculates the L-n distance between 2 vectors
	// TAKES INTO ACCOUNT THE "weights" factor
	
	int n=x.size();
	assert(n == y.size());
	assert(n == weights.size());
	
	double res=0.0; 
	
	for(int i=0; i<n; i++)
	{
		res += weights[i] * pow(x[i]-y[i],power);
	}
	
	return pow(res,1.0/power);
}



double normLvector(int power, vector<double> x)
{
	// Calculates the L-n norm of a vector
	
	int n=x.size();
	
	double res=0.0; 
	
	for(int i=0; i<n; i++)
		res += pow(x[i],power);
	
	return pow(res,1.0/power);
}


double powerLvector(int power, vector<double> x,vector<double> y )
{
	// Calculates the L-n distance^n between 2 vectors
	
	
	int n=x.size();
	assert(n == y.size());
	
	double res=0.0;
	
	for(int i=0; i<n; i++)
		res += pow(x[i]-y[i],power);

	return res;
}

double maxElementVector(vector<double> x)
{
	/// return the maximum element of a vector
    double max=-999999999;
    
    for (int i=0; i<x.size(); i++) {
        if (max < x[i]) {
            max = x[i];
        }
    }
    return max;
}

double minElementVector(vector<double> x)
{
    double min = 999999999;
    
    for (int i=0; i<x.size(); i++) 
	{
        if (min > x[i]) 
		{
            min = x[i];
        }
    }
    return min;
}

int argminElementVector(vector<double> x)
{
    double min = 999999999;
	int argmin = 0;
    
    for (int i=0; i<x.size(); i++) 
	{
        if (min > x[i]) 
		{
            min = x[i];
			argmin = i;
        }
    }
    return argmin;
}


vector<double> max_vector(vector<double> x, vector<double> y)
{
	/// Return a vector composed of the max element by element of x and y
	stopif((x.size()!=y.size()), "vectors must be same size!");
	stopif(x.size()==0,"vector is empty!");
	
	vector<double> res;
	for(int i=0; i<x.size(); i++)
		res.push_back(max(x[i],y[i]));
	return res;
}



double standardDeviation(vector<double> x)
{
	// Returns the standard deviation of all elements
	
	double avg = averageElements(x);
	int N = x.size();
	
	double sd = 0 ;
	
	for (int i=0; i<N; i++)
	{
		sd += (x[i]-avg)*(x[i]-avg);
		
	}
	
	sd = sd/(N-1); // (N-1) not N, because unbiaised estimator
	return sqrt(sd);
}
double standardDeviation(vector<unsigned long int> x)
{
	// Returns the standard deviation of all elements
	
	double avg = averageElements(x);
	int N = x.size();
	
	double sd = 0 ;
	
	for (int i=0; i<N; i++)
	{
		sd += (x[i]-avg)*(x[i]-avg);
		
	}
	
	sd = sd/(N-1); // (N-1) not N, because unbiaised estimator
	return sqrt(sd);
}

/* **** **** */

vector<double> smoothVector(vector<unsigned long int> x, int lagSmooth)
{
	int N = x.size()-2*lagSmooth;
	if (N<=1)
	{
		cout<< endl << "CANNOT SMOOTH: VECTOR SIZE TOO SMALL AND LAG TOO LARGE"<<endl;
		exit(1);
	}
	
	vector<double> s(N);
	
	bool debug = false;
	if(debug) 
	{
		cout<<"[in smooth] x size="<<x.size()
		<<" ; s size="<<N<<" ; lagSmooth="<<lagSmooth <<endl;
	}
	
	for (int i=lagSmooth; i < x.size()-lagSmooth ; i++) 
	{
		unsigned long int sum=0;
		
		for (int k=-lagSmooth; k<=lagSmooth; k++) 
		{
			sum += (double)(x[i+k]);
			//cout << "sum_"<<i<<"_"<<k<<" = "<<sum<<endl; // DEBUG
		}
		
		s[i-lagSmooth] = sum/(2*lagSmooth+1);
		if(debug) {cout<<"s["<<i-lagSmooth<<"]="<<s[i-lagSmooth] <<endl;}	// DEBUG
	}
	
	return s;
}


vector<double> numericalDerivative(vector<double> x, vector<double> t)
{
	// dX/dt
	
	int n = x.size();
	if (t.size()!=n) 
	{
		cout << " numericalDerivative ERROR: vectors not the same size ***"<<endl;
		exit(1);
	}
	
	vector<double> dx(n-1);
	
	for (int i=0; i<n-1; i++) 
	{
		dx[i] = (x[i+1]-x[i])/(t[i+1]-t[i]);
	}
	
	return dx;
}


vector<double> numericalSecondDerivative(vector<double> x, vector<double> t, int lagDeriv = 1)
{
	// dX/dt
	
	int n = x.size();
	if (t.size()!=n) exit(1);
	
	if (n<=3) 
	{
		cout << endl<< "ERROR NUM 2nd DERIV: Vector not long enough!";
		exit(1);
	}
	
	int sizeDeriv2 = n-2*lagDeriv;
	vector<double> ddx(sizeDeriv2);
	
	for (int i=1; i<=sizeDeriv2; i++) 
	{
		ddx[i-1] = (x[i+lagDeriv] - 2*x[i] + x[i-lagDeriv]) / (t[i+lagDeriv]-t[i])/(t[i+lagDeriv]-t[i]);
	}
	
	return ddx;
}



double findAvgMax(vector<unsigned long int> x, vector<double> t, int lagSmooth, double minSlopeAboutMax)
{
	
	/// Find the average maximum among "noisy" values
	/// (smooth data first)
	bool debug = false;
	
	if (debug){
		cout<<" --- FindAvgMax --- "<<endl; //DEBUG
		cout << " Lag = "<< lagSmooth<<endl;
		cout<<"-inside x:";displayVector(x);
		cout<<"-inside t:";displayVector(t);
	}
	
	// Extract the right elements of the "time" vector 
	int n = t.size();
	
	vector<double> tt(n-2*lagSmooth);
	
	for (int i=lagSmooth; i<n-lagSmooth; i++) 
	{
		tt[i-lagSmooth] = t[i];
	}
	
	vector<double> s	= smoothVector(x, lagSmooth);
	if(debug) { cout << "smoothed vect:"; displayVector(s);}
	
	vector<double> ds	= numericalDerivative(s, tt);
	if(debug) { cout << "deriv vect:"; displayVector(ds);}
	
	int lagDeriv = 1;
	
	vector<double> dds	= numericalSecondDerivative(s, tt,lagDeriv);
	if(debug) { cout << "deriv2 vect:"; displayVector(dds);}
	
	if (1)
	{
		vectorToFile(s, "tmpSmooth.txt");
		vectorToFile(tt, "tmpSmooth_t.txt");
	}
	if (debug)
	{
		cout<<"-inside s:";displayVector(s);
		cout<<"-inside ds:";displayVector(ds);
		cout<<"-inside dds:";displayVector(dds);
	}
	
	vector<unsigned long int> maximum;
	vector<double> arg_maximum;
	
	// Looking for extrema
	
	int cnt = 0;
	
	for (int i=0; i<dds.size(); i++) 
	{		
		if(debug){cout<<i<<"_deriv= "<<(ds[i])<<endl;
			cout<<i<<"_deriv2= "<<(dds[i])<<endl;}
		
		if (ds[i]*ds[i+1]<=0 
			&&	dds[i]<0 
			&&	abs(ds[i])>= minSlopeAboutMax
			&&	abs(ds[i+1])>= minSlopeAboutMax) 
		{
			// If first derivative changes sign 
			// and 2nd deriv neg => it's a max
			// To avoid false positive, option to check 
			// the slope about the max (with test on minSlopeAboutMax)
			
			maximum.push_back(x[i+1+lagSmooth]);
			arg_maximum.push_back(t[i+1+lagSmooth]);
			cnt++;
			
			if (debug)	cout << "TMP_MAX_"<<i<<" = "<<x[i]<<endl; // DEBUG
		}
	}
	
	if(1) 
	{
		vectorToFile(maximum, "findAvgMax.out");
		vectorToFile(arg_maximum, "findAvgMax_arg.out");
	}
	
	double res=-9e17;
	if (cnt>0) res = averageElements(maximum);
	
	// If no maximum found:
	if (cnt==0 && ds[0]>0) res = x[x.size()-1];
	if (cnt==0 && ds[0]<0) res = x[0];
	
	return res;
}


double findAvgMin(vector<unsigned long int> x, vector<double> t, int lagSmooth, double minSlopeAboutMin)
{
	// NOTE: Mirror the findAvgMax code
	
	// Find the average maximum among "noisy" values
	// (smooth data first)
	bool debug = false;
	
	if (debug){
		cout<<" --- FindAvgMin --- "<<endl; //DEBUG
		cout << " Lag = "<< lagSmooth<<endl;
		cout<<"-FindAvgMin x:";displayVector(x);
		cout<<"-FindAvgMin t:";displayVector(t);
	}
	
	// Extract the right elements of the "time" vector 
	int n = t.size();
	
	vector<double> tt(n-2*lagSmooth);
	
	for (int i=lagSmooth; i<n-lagSmooth; i++) 
	{
		tt[i-lagSmooth] = t[i];
	}
	
	vector<double> s	= smoothVector(x, lagSmooth);
	if(debug) { cout << "smoothed vect:"; displayVector(s);}
	
	vector<double> ds	= numericalDerivative(s, tt);
	if(debug) { cout << "deriv vect:"; displayVector(ds);}
	
	int lagDeriv = 1;
	
	vector<double> dds	= numericalSecondDerivative(s, tt,lagDeriv);
	if(debug) { cout << "deriv2 vect:"; displayVector(dds);}
	
	if (1)
	{
		vectorToFile(s, "tmpSmooth.txt");
		vectorToFile(tt, "tmpSmooth_t.txt");
	}
	if (debug)
	{
		cout<<"-FindAvgMin s:";displayVector(s);
		cout<<"-FindAvgMin ds:";displayVector(ds);
		cout<<"-FindAvgMin dds:";displayVector(dds);
	}
	
	vector<unsigned long int> minimum;
	vector<double> arg_minimum;
	
	// Looking for extrema
	
	int cnt = 0;
	
	for (int i=0; i<dds.size(); i++) 
	{		
		if(debug){cout<<i<<"_deriv= "<<(ds[i])<<endl;
			cout<<i<<"_deriv2= "<<(dds[i])<<endl;}
		
		if (ds[i]*ds[i+1]<=0 
			&&	dds[i]>0 
			&&	abs(ds[i])>= minSlopeAboutMin
			&&	abs(ds[i+1])>= minSlopeAboutMin) 
		{
			// If first derivative changes sign 
			// and 2nd deriv neg => it's a max
			// To avoid false positive, option to check 
			// the slope about the max (with test on minSlopeAboutMax)
			
			minimum.push_back(x[i+1+lagSmooth]);
			arg_minimum.push_back(t[i+1+lagSmooth]);
			cnt++;
			
			if (debug)	cout << "TMP_MIN_"<<i<<" = "<<x[i]<<endl; // DEBUG
		}
	}
	
	if(1) 
	{
		vectorToFile(minimum, "findAvgMin.out");
		vectorToFile(arg_minimum, "findAvgMin_arg.out");
	}
	
	double res=9e17;
	if (cnt>0) res = averageElements(minimum);
	
	// If no min found:
	if (cnt==0 && ds[0]>0) res = x[x.size()-1];
	if (cnt==0 && ds[0]<0) res = x[0];
	
	return res;
}




vector<double> vector_seq(double start, double end, unsigned int size)
{
	/// Create a vector of size 'size' with values
	/// uniformly partitioned between 'start' and 'end'
	
	vector<double> x;
	
	for (unsigned int i=0; i<size; i++)
		x.push_back(start + (end-start)/(size-1)*(double)(i));
	
	return x;
}



vector<double> vector_seq_by(double start, double end, double by)
{
	/// Create a vector uniformly partitioned between 'start' and 'end'
	vector<double> x;
	
	x.push_back(start);
	
	for(unsigned int i=0; (x[i]+by)<end; ) {
		x.push_back(x[i]+by);
		i++;
	}
	return x;
}




// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////

void my_pause(double waitTime)
{
	time_t t1,t2;
	
	t1=time(NULL);
	t2=t1;
	
	while (difftime(t2,t1)<waitTime) {
		t2=time(NULL);
	}
	
	cout << "time1= "<<t1<<endl;
	cout << "time2= "<<t2<<endl;	
}



// ===================================================
// ====================  SAMPLING  ===================
// ===================================================

dcMatrix LatinHypercubeSampling(vector<double> Vmin, vector<double> Vmax,
                              int samplingNb)
{
	// For a given Min and Max values for each parameters,
	// returns a shuffled dcMatrix
	
	// Initialize the seed
	timeval tim;
    gettimeofday(&tim, NULL);
	srand(tim.tv_usec);
	
    assert(Vmin.size()==Vmax.size());

    int N = Vmin.size();
    
    dcMatrix V(samplingNb,N);
    
    // Populate values (in order)
	// _UNIFORMLY_ distributed
	
    for (int i=0; i<samplingNb; i++) 
    {
        for (int j=0; j<N; j++) 
        {
            V(i,j) = Vmin[j]+i*(Vmax[j]-Vmin[j])/(samplingNb-1);
        }
    }
    
    dcMatrix W(samplingNb,N); 
    
	unsigned int theseed = tim.tv_usec; // Fix the seed if reproductibility desired
    srand ( theseed );
    
    for (int j=0; j<N; j++) 
    {
		// Shuffle every elements of each column of V
		
        vector<double> x = V.extractColumn(j);
        random_shuffle(x.begin(), x.end());
        W.setColumnValues(j, x);
    }    
    return W;
}






vector<double> partitionLinear(double a_min, double a_max, int nbPartitions)
{
	vector<double> a(nbPartitions);
	
	assert(nbPartitions>1);
	
	for (int i=0; i<nbPartitions; i++)
	{
		a[i] = a_min + i*(a_max-a_min)/(nbPartitions-1);
	}
	
	return a;
}


vector<double> partitionLog(double a_min, double a_max, int nbPartitions)
{
	vector<double> a(nbPartitions);
	
	assert(nbPartitions>1);
	
	for (int i=0; i<nbPartitions; i++)
	{
		a[i] = exp(log(a_min) + i*(log(a_max)-log(a_min))/(nbPartitions-1));
	}
	
	return a;
}


vector<int>	uniformIntVector(int seed, int size, int min, int max)
{
	vector<int> res;
	
	for (int i=0; i<size; i++) 
	{
		res.push_back(uniformInt(min,max));
	}
	
	return res;
}




vector< long> uniformIntVectorUnique(long size, long min, long max)
{
	// Returns a vector of size "size",
	// with elements integers randomly
	// drawn from [min;max]
	// == No redundant elements ==
	
	// Check if it is possible to have uniqueness
	if (size > max-min+1)
	{
		cout << endl <<  "*** ERROR [uniformIntVectorUnique]: Cannot guarantee uniqueness";
		cout << "(size="<<size<<"> range:"<< max-min+1<<endl;
		exit(1);
	}
	
	// Vector x: min,min+1,...,max
	vector<long> x(max-min+1);
	x[0] = min;
	for (int i=1; i<x.size(); i++)
		x[i] = x[i-1]+1;
	
	vector<long> res;
	
	for (int i=0; i<size; i++)
	{
		long v_i = extractElementRandom(x);
		res.push_back(v_i);
		// pop element to ensure uniqueness
		x = popElementValue(x, v_i);
	}
	return res;
}

// ===================================================
// ==================== CONVERSION ===================
// ===================================================



string int2string(int i)
{
	string res;					// string which will contain the result
	ostringstream convert;		// stream used for the conversion
	convert << i;				// insert the textual representation of 'Number' in the characters in the stream
	res = convert.str();		// set 'Result' to the contents of the stream
	
	return res;	
}


// ===================================================
// ==================== MATH FUNCTIONS ===================
// ===================================================


// OBSOLETE?

//double pseudo_beta_density(double x, double a, double b)
//{
//	double cst = pow(((a-1)/(a+b-2)),(a-1))*pow((b-1)/(a+b-2),(b-1));
//	return(pow(x,(a-1))*pow((1-x),(b-1))/cst);
//}



// ===================================================
// ==== Compensation ====
// ===================================================



double comp_odds(double q, double nq, double xq)
{
	double r = (log(q)-log(nq))/(log(xq)-log(nq));
	
	return log(r/(1-r));
}

double comp_orig(double L, double nq, double xq)
{
	double r = 1/(1+exp(-L));
	
	return exp(log(nq) + r*(log(xq)-log(nq)));
}


double comp_analytic_root(double a,double na,double xa,
						  double b,double nb,double xb,
						  double c_new)
{
	double alpha_a = log(xa/na);
	double alpha_b = log(xb/nb);
	
	double r_a = log(a/na)/log(xa/na);
	double r_b = log(b/nb)/log(xb/nb);
	
	double beta_a = 1/r_a -1;
	double beta_b = 1/r_b -1;
	
	double K_ab = log(c_new/(na*nb));
	
	double A = K_ab*beta_a*beta_b;
	double B = beta_a*(K_ab-alpha_b) + beta_b*(K_ab-alpha_a);
	double C = K_ab - alpha_a-alpha_b;
	
	double Y_plus = (-B+sqrt(B*B-4*A*C))/2/A;
	double Y_minus = (-B-sqrt(B*B-4*A*C))/2/A;
	
	double res=0;
	
	if(Y_minus>0) res = -log(Y_minus);
		
	if(Y_plus>0) res = -log(Y_plus);
	
	return res;
}



vector<double> comp_adj(double a, double min_a, double max_a,
						double b, double min_b, double max_b,
						double c_new)
{
	double D = comp_analytic_root(a, min_a, max_a, b, min_b, max_b, c_new);
	
	double a_new = comp_orig(D + comp_odds(a, min_a, max_a), min_a, max_a);
	double b_new = comp_orig(D + comp_odds(b, min_b, max_b), min_b, max_b);
	
	vector<double> res(2);
	res[0] = a_new;
	res[1] = b_new;
	
	return res;
}





// ====== MATHS FUNCTIONS ======



double pseudo_beta(double x, double a, double b)
{
	/// Function inspired from the density
	/// of the beta distribution
	/// (maximum value [not integral] of fct is 1.0)
	
	double res = 0;
	
	if (x>=0 && x<=1){
	double cst =  pow((a-1)/(a+b-2),(a-1)) * pow((b-1)/(a+b-2),(b-1));
	res = pow(x,a-1) * pow(1-x,b-1) / cst;
	}
	
	return res;
}


double pseudo_gamma(double x, double shape, double scale)
{
	/// Function inspired from the density
	/// of the gamma distribution
	/// (maximum value [not integral] of fct is 1.0)
	
	stopif(x<=0, "Negative value not allowed");
	stopif(shape<=0, "Negative value for shape not allowed");
	stopif(scale<=0, "Negative value for scale not allowed");
	
	double tmp1 = pow(x,shape-1) * exp(-x/scale);
	double tmp2 = factorial((int)(shape)) * pow(scale,shape);
	
	return tmp1/tmp2;
}


