//
//  dcMatrix.cpp
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "dcMatrix.h"
#include "dcTools.h"
 

// ////////////////////////////////////////////////////////////////////
//              CONSTRUCTORS
// ////////////////////////////////////////////////////////////////////


dcMatrix::dcMatrix(string fileName)
{
    ifstream thefile (fileName.c_str());
	
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			thefile >> val[nbCols*i+j];
		}
	
}



dcMatrix::dcMatrix(vector<double> v)
{
	unsigned int nCol = v.size();
	nbRows = nCol;
	nbCols = 1;
	
	val.resize(nCol);
	
	for (unsigned int i=0; i<nCol; i++) 
	{
		val[i] = v[i];
	}
}


// ////////////////////////////////////////////////////////////////////
//              OPERATORS
// ////////////////////////////////////////////////////////////////////

double & dcMatrix::operator () (unsigned int i,unsigned int j)
{
	/// Retrieve matrix element
	/// top left element is M(0,0)
	/// bottom right element is M(n-1,m-1)
	
    assert(i>=0);assert(i<nbRows);
    assert(j>=0);assert(j<nbCols);
    
    return val[nbCols*i+j];
}

double & dcMatrix::operator [] (unsigned int i)
{
    assert(i>=0);assert(i<nbRows);assert(nbCols==1);
    return val[i];
}


// ////////////////////////////////////////////////////////////////////
//              Set Functions
// ////////////////////////////////////////////////////////////////////



void dcMatrix::RandomInit()
{
	srand ( time(NULL) );
	
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			val[nbCols*i+j]=rand() % 11 - 10;
		}
    
}

void dcMatrix::setAllValues(double value)
{
    for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = value;
		}
    
}

void dcMatrix::setValueFromMatrix(dcMatrix M)
{
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = M(i,j);
		}
}


// ////////////////////////////////////////////////////////////////////
//              File Functions
// ////////////////////////////////////////////////////////////////////


void dcMatrix::FromFile(const char* fileName)
{
	ifstream thefile (fileName);
	
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			thefile >> val[nbCols*i+j];
		}
	
}

void dcMatrix::FromFile(string fileName)
{
	/// READ A FILE AND FILL THE
	/// ELEMENT VALUES TO THE MATRIX
	/// *** WARNING *** MUST BE SAME SIZE!!!
	
	ifstream thefile (fileName.c_str());
	
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
			thefile >> val[nbCols*i+j];
}

void dcMatrix::FromFile_Rows(string fileName, unsigned int nrow)
{
	/// READ A FILE AND FILL THE
	/// ELEMENT VALUES TO THE MATRIX

	/// NUMBER OF ROWS PRE-SPECIFIED
	
	ifstream thefile (fileName.c_str());
	
	if (!thefile)
	{
		cout<<endl<<" ERROR [FromFile_Rows]: This file is not found: "<<fileName<<endl;
		exit(1);
	}
	
	vector<double> x;	
	double file_value;
	
	x.clear();
	
	while (!thefile.eof())   
	{
		thefile >> file_value;
		if (thefile.eof()) break;
		x.push_back (file_value);	
	}
	
	thefile.close();
	
	unsigned int n = x.size();
	
	if (n%nrow != 0 )
	{
		cout << endl << "ERROR [FromFile_Rows]: number of rows ("<< nrow<<") does not divide number of data("<<n<<")!"<<endl;
		exit(1);
	}

	this->resize(nrow, n/nrow);
	
	unsigned int cnt=0;
	
	for(unsigned int i=0;i<nbRows;i++)		
		for(unsigned int j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = x[cnt];
			cnt++;
		}	
	
}



void dcMatrix::WriteToFile(string fileName)
{
	ofstream thefile (fileName.c_str());
	
	for(unsigned int i=0;i<nbRows;i++)
    {
		for(unsigned int j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] << "\t";
		}
        thefile << endl;
    }
}

void dcMatrix::WriteToFileCSV(string fileName)
{
	ofstream thefile (fileName.c_str());
	
	for(unsigned int i=0;i<nbRows;i++)
    {
		for(unsigned int j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] ;
			if (j<nbCols-1) thefile<< ",";
		}
        thefile << endl;
    }
}

void dcMatrix::WriteToFileCSV(string fileName, vector<string> headers)
{
	ofstream thefile (fileName.c_str());
	
	string err_msg = "Headers size (" + to_string(headers.size()) + ") does not match dcMatrix columns size (" + to_string(nbCols) + ")!\n";
	err_msg = err_msg + "Cannot write matrix to this file: " + fileName;
	stopif(nbCols!=headers.size(), err_msg);
	
//	if (nbCols!=headers.size())
//	{
//		cout << "ERROR -- WriteToFileCSV: Headers size ("<<headers.size()<<") does not match dcMatrix columns size ("<<nbCols<<")!"<<endl;
//		cout << "Cannot write matrix to this file: "<<fileName<<endl;
//		exit(1);
//	}
	
	// Write headers on the first line
	for (unsigned int j=0; j<nbCols; j++) 
	{
		thefile << headers[j];
		if (j<nbCols-1) thefile<< ",";
	}
	thefile << endl;
	
	// Then write data
	for(unsigned int i=0;i<nbRows;i++)
    {
		for(unsigned int j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] ;
			if (j<nbCols-1) thefile<< ",";
		}
        thefile << endl;
    }
}



// ////////////////////////////////////////////////////////////////////
//              Vector Functions
// ////////////////////////////////////////////////////////////////////



void dcMatrix::addRowVector(vector<double> v)
{    
	dcMatrix res;
	
	if (nbRows==0 && nbCols==0)
	{
		// If dcMatrix is empty, then this vector initialize the dcMatrix
		// as a 1 row matrix
		
		dcMatrix tmp(1,v.size());
		
		for (unsigned int j=0; j<v.size(); j++) tmp(0,j) = v[j];
		res = tmp;
	}
	
	else 
	{
		// If dcMatrix is NOT empty
		
		if(nbCols != v.size())
		{
			cout<<"CAN'T ADD ROW VECTOR TO MATRIX, SIZES DO NOT MATCH: dcMatrix cols = "<< nbCols
			<<" vs Vector size = "<< v.size() << endl;
			exit(1);
		}
		
		dcMatrix tmp(nbRows+1,nbCols);
		
		for(unsigned int i=0; i<nbRows; i++)
		{
			for(unsigned int j=0; j<nbCols;j++)
			{
				tmp(i,j) = val[i*nbCols+j];//this->val[i*ncol+j];
			}
		}
		
		for(unsigned int j=0; j<nbCols;j++)
		{
			tmp(nbRows,j) = v[j];
		}
		
		res=tmp;
	}


    
	*this = res;
}


void dcMatrix::addRowVector(vector<unsigned long> v)
{
	vector<double> tmp;
	
	for (unsigned int i=0; i<v.size(); i++)
		tmp.push_back((double)(v[i]));
	
	addRowVector(tmp);
}


void dcMatrix::addColVector(vector<double> v)
{        
    unsigned int nrow = this->nbRows;
    unsigned int ncol = this->nbCols;

	string errmsg = "CAN'T ADD Col VECTOR TO MATRIX, SIZES DO NOT MATCH: dcMatrix rows = ";
	errmsg = errmsg + to_string(nrow) + " vs Vector size = " + to_string(v.size());

	
	stopif(nrow != v.size(),errmsg);
	
    dcMatrix tmp(nrow,ncol+1);
    
    for(unsigned int i=0; i<nrow; i++){
        for(unsigned int j=0; j<ncol;j++)
            tmp(i,j) = this->val[i*ncol+j];
    }
    
    for(unsigned int i=0; i<nrow; i++) tmp(i,ncol) = v[i];
    
    *this = tmp;
}

void dcMatrix::removeRow(unsigned int i_row)
{
	unsigned int c=0;
	for (unsigned int j=nbCols-1; c<nbCols; j--){
		val.erase(val.begin() + i_row*nbCols+j);
		c++;
	}
	nbRows--;
}


void dcMatrix::removeCol(unsigned int j_col)
{
	unsigned int c=0;
	for (unsigned int i=nbRows-1; c<nbRows; i--){
		val.erase(val.begin() + i*nbCols+j_col);
		c++;
	}
	nbCols--;
}


vector<double> dcMatrix::extractColumn(unsigned int j_col)
{
	string errmsg ="Cannot extract dcMatrix column("+ to_string(j_col)+") greater than size (0--"+ to_string(nbCols-1)+")!";
	stopif (j_col >= nbCols, errmsg);
	
	vector<double> v(nbRows);
	
    for(unsigned int i=0; i<nbRows; i++) v[i]= val[i*nbCols+ j_col];
    
	return v;
}

vector<double> dcMatrix::extractRow(unsigned int i_row)
{
	if (i_row >= nbRows) {
		cout << endl << " ERROR [dcMatrix::extractRow]:cannot extract row("<< i_row;
		cout << ") greater than size (0--"<< nbRows-1<< ")!"<<endl;
		exit(1);
	}
	
	vector<double> v(nbCols);
	
    for(unsigned int j=0; j<nbCols; j++) 
		v[j]= val[i_row*nbCols+ j];
    
	return v;
}


vector<double>	dcMatrix::extractRow_cond_minElement(unsigned int j_col)
{
	if (j_col >= nbCols) {
		cout << endl << " ERROR [dcMatrix::extractRow_cond_minElement]:cannot extract row("<< j_col;
		cout << ") greater than size (0--"<< nbCols-1<< ")!"<<endl;
		exit(1);
	}
	
	vector<double> col_select = extractColumn(j_col);
	
	//argminElementVector(col_select); DOESN'T WORK WHEN #include "dcTools.h" ???
	
	double min = 999999999;
	unsigned int argmin = 0;
    
    for (unsigned int i=0; i<col_select.size(); i++) 
	{
        if (min > col_select[i]) 
		{
            min = col_select[i];
			argmin = i;
        }
    }
	
	return extractRow(argmin);
}




void dcMatrix::setColumnValues(unsigned int colNb_start0, vector<double> v)
{
	string errmsg = "dcMatrix column size ("+to_string(nbRows)+") is different from vector size ("+to_string(v.size())+")" ;
    stopif(v.size() != nbRows,errmsg);
    
    for(unsigned int i=0; i<nbRows; i++) 
        val[i*nbCols+ colNb_start0] = v[i] ;
}


void dcMatrix::setRowValues(unsigned int rowNb_start0, vector<double> v)
{
	/// Set the values of a given row
	
	string errmsg = "dcMatrix row size ("+to_string(nbCols)+") is different from vector size ("+to_string(v.size())+")" ;
	stopif(v.size() != nbCols,errmsg);
    
    for(unsigned int j=0; j<nbCols; j++) 
        val[rowNb_start0*nbCols+ j] = v[j] ;
}


void dcMatrix::setRowValues(unsigned int rowNb_start0, vector<unsigned long> v)
{
	/// Set the values of a given row
	
	string errmsg = "dcMatrix row size ("+to_string(nbCols)+") is different from vector size ("+to_string(v.size())+")" ;
	stopif(v.size() != nbCols,errmsg);
	
	for(unsigned int j=0; j<nbCols; j++)
		val[rowNb_start0*nbCols+ j] = (double)(v[j]) ;
}




// ////////////////////////////////////////////////////////////////////
//              Usual dcMatrix Functions
// ////////////////////////////////////////////////////////////////////



void dcMatrix::display()
{
    cout<<endl;
	cout << "dcMatrix dimension: "<<nbRows<<"x"<<nbCols<<endl;
    for(unsigned int i=0;i<nbRows;i++)
    {
        cout<<"[ ";
        for(unsigned int j=0;j<nbCols;j++)
        {
            cout<<val[i*nbCols+j];
			if (j<nbCols-1) cout<<"\t";
        }
        cout<<"]"<<endl;
    }
}

bool dcMatrix::isSymetric()
{
	if (nbCols!=nbRows) return false;
	unsigned int nCol = nbRows;
    unsigned int i=1,j=0;
	
	if (nCol==1) return 1;
	
    
	for (i=1;i<nCol;i++)
        for (j=0;j<i;j++)
            if ( val[i*nbCols+j]!=val[j*nbCols+i] ) {break;}
	if ( (i==nCol) && (j==nCol-1) ) return true;
	else return false;
}


dcMatrix dcMatrix::transpose()
{
    dcMatrix B(nbCols,nbRows);
    
    for (unsigned int i=0;i<nbCols;i++)
        for (unsigned int j=0;j<nbRows;j++)
        {
            B(i,j)=val[j*nbCols+i];
        }
    return B;
}


double dcMatrix::determinant()
{
    assert (nbCols==nbRows);

    double s = 0;
    unsigned int nCol = nbCols;

    if (nCol==2) return val[0*nbCols+0]*val[1*nbCols+1]
                        - val[0*nbCols+1]*val[1*nbCols+0]; //A(0,0)*A(1,1)-A(0,1)*A(1,0);
    
    else
    {
        dcMatrix B(nCol-1);
        for(unsigned int k=0;k<nCol;k++)
        {
            for(unsigned int i=0;i<nCol;i++)
                for(unsigned int j=1;j<nCol;j++) 
                { 
                    if (i<k) B(i,j-1)=val[i*nbCols+j];
                    if (i>k) B(i-1,j-1)=val[i*nbCols+j];
                }
            s += val[k*nbCols+0]*B.determinant()*pow (-1,k);
        } 
        return s;
    }
}

dcMatrix dcMatrix::getMinor(unsigned int row, unsigned int col)
{
	if (nbCols!=nbRows)
	{
		cout << endl<< "ERROR [dcMatrix::getMinor] : non square matrix" <<endl;
		exit(1);
	}
	
	unsigned int n = nbCols;
	
	dcMatrix B(n-1);
	// calculate the cofactor of element (row,col)
	
	// indicate which col and row is being copied to dest
	unsigned int colCount=0,rowCount=0;
	
	for(unsigned int i = 0; i < n; i++ )
	{
		if( i != row )
		{
			colCount = 0;
			for(unsigned int j = 0; j < n; j++ )
			{
				// when j is not the element
				if( j != col )
				{
					B(rowCount,colCount) = val[i*n+j];
					colCount++;
				}
			}
			rowCount++;
		}
	}

	return B;
}


dcMatrix dcMatrix::inverse()
{	
	if (nbRows!=nbCols)
	{
		cout << "ERROR [dcMatrix::inverse] : cannot inverse a non-square matrix!"<<endl;
		exit(1);
	}
	
	unsigned int n = nbCols;
	
	// Copy this dcMatrix unsigned into A
	dcMatrix A(n);
	for (unsigned int i=0; i<n; i++) {
		for (unsigned int j=0; j<n; j++) {
			A(i,j) = val[n*i+j];
		}
	}
	
	// get the determinant of A
	double det = A.determinant();
	
	if (det==0)
	{
		cout << "ERROR [dcMatrix::inverse] : cannot inverse matrix (determinant	=0)!"<<endl;
		exit(1);
	}
	
    double oneOverDet = 1.0/det;
	
	dcMatrix AInverse(n);
	
	if (n==2)
	{
		AInverse(0,0) = A(1,1)*oneOverDet;
		AInverse(0,1) = -A(0,1)*oneOverDet;
		AInverse(1,0) = -A(1,0)*oneOverDet;
		AInverse(1,1) = A(0,0)*oneOverDet;
	}
	
	if (n>2)
	{
		dcMatrix Aminor(n-1);
		
		for(unsigned int j=0;j<n;j++)
		{
			for(unsigned int i=0;i<n;i++)
			{
				// get the co-factor (matrix) of A(j,i)
				Aminor = A.getMinor(j,i);
				
				AInverse(i,j) = oneOverDet*Aminor.determinant(); 
				
				if( (i+j)%2 == 1)
					AInverse(i,j) = -AInverse(i,j);
			}
		}
	}
	
	return AInverse;
}



dcMatrix dcMatrix::Cholesky()
{				
    if(nbCols!=nbRows) 		
    { cout<<"Cholesky: non square dcMatrix !"<<endl;
        exit(1);
    }
    
    if( ! this->isSymetric() ) 
    {  cout<<"Cholesky: non symetric dcMatrix !"<<endl;
        exit(1);
    }
	
    unsigned int nCol = nbCols;
    dcMatrix L(nCol);			
    float s;
    
    for(unsigned int i=0;i<nCol;i++)
    {
        s=0;
        for(unsigned int k=0;k<i;k++) s+=L(i,k)*L(i,k);
        L(i,i) = sqrt(val[i*nbCols+i]-s);
        
        for(unsigned int j=i+1;j<nCol;j++) 
        {
            s=0;
            for(unsigned int k=0;k<i;k++) s+=L(j,k)*L(i,k);
            L(j,i)=(val[i*nbCols+j]-s)/L(i,i); 
        }
    }
    return L;
} 


double dcMatrix::getMinimumValue()
{
    double min = 1e36;
    
    for (unsigned int i=0;i<nbCols;i++)
        for (unsigned int j=0;j<nbRows;j++)
        {
            double tmp = val[i*nbCols+j];
            if (tmp<min) min=tmp;
            
        }
    return min;
}


double dcMatrix::getMaximumValue()
{
    double maxi = -1e36;
    
    for (unsigned int i=0;i<nbCols;i++)
        for (unsigned int j=0;j<nbRows;j++)
        {
            double tmp = val[i*nbCols+j];
            if (maxi<tmp) maxi=tmp;
            
        }
    return maxi;
}

double dcMatrix::sumAllElements()
{
	double s=0;
	
	for(unsigned int i=0;i<nbRows;i++)
	{
		for(unsigned int j=0;j<nbCols;j++)
		{
			s=s+val[i*nbCols+j];
		}
        
	}
	
	return s;
}

double dcMatrix::sumLine(unsigned int i)
{
	double s=0.0;
	
	for(unsigned int k=0;k<nbCols;k++)
	{
        s += val[i*nbCols+k];
	}
	
	return s;
}

double dcMatrix::sumColumn(unsigned int j)
{
	double s=0.0;
	
	for(unsigned int k=0;k<nbRows ;k++)
	{
        s += val[k*nbCols+j];
	}
	
	return s;
}


double dcMatrix::sumColumnIf(unsigned int colToSum, unsigned int colToTest,
						   double lowerBound, double upperBound)
{
	double s=0.0;
	
	for(unsigned int k=0;k<nbRows ;k++)
	{
		if (val[k*nbCols+colToTest]>lowerBound &&
			val[k*nbCols+colToTest]<upperBound)
			
			s += val[k*nbCols+colToSum];
	}
	
	return s;
}

unsigned int dcMatrix::countColumnIf(unsigned int colToTest, double lowerBound, double upperBound)
{
	unsigned int c = 0;
	
	for(unsigned int k=0;k<nbRows ;k++)
	{
		if (val[k*nbCols+colToTest]>lowerBound &&
			val[k*nbCols+colToTest]<upperBound)
			
			c ++;
	}
	
	return c;
}
	

/* ************************************************************ */


dcMatrix operator + (dcMatrix &A,dcMatrix &B)
{
    assert(A.nbCols==B.nbCols && A.nbRows==B.nbRows);
    dcMatrix C(A.nbRows,A.nbCols);
    for (unsigned int i=0;i<A.nbRows;i++)
        for (unsigned int j=0;j<A.nbCols;j++)
        {
            C(i,j)=A(i,j)+B(i,j);
        }
    return C;   
}

dcMatrix operator - (dcMatrix &A,dcMatrix &B)
{
    assert(A.nbCols==B.nbCols && A.nbRows==B.nbRows);
    dcMatrix C(A.nbRows,A.nbCols);
    for (unsigned int i=0;i<A.nbRows;i++)
        for (unsigned int j=0;j<A.nbCols;j++)
        {
            C(i,j)=A(i,j)-B(i,j);
        }
    return C;   
}


dcMatrix operator * (dcMatrix &A,dcMatrix &B)
{
    assert(A.nbCols==B.nbRows);
    dcMatrix C(A.nbRows,B.nbCols);
    double s=0;
    for (unsigned int i=0;i<A.nbRows;i++)
        for (unsigned int j=0;j<B.nbCols;j++)
        {
            s=0;
            for(unsigned int k=0;k<A.nbCols;k++) s+=A(i,k)*B(k,j);
            C(i,j)=s;
        }
    return C;   
}

dcMatrix operator * (double a,dcMatrix &A)
{
    dcMatrix C(A.nbRows,A.nbCols);
    for (unsigned int i=0;i<A.nbRows;i++)
        for (unsigned int j=0;j<A.nbCols;j++)
        {
            C(i,j)=a*A(i,j);
        }
    return C;   
}





dcMatrix transpo(dcMatrix A)
{
    dcMatrix B(A.nbCols,A.nbRows);
    for (unsigned int i=0;i<A.nbCols;i++)
        for (unsigned int j=0;j<A.nbRows;j++)
        {
            B(i,j)=A(j,i);
        }
    return B;
}

unsigned int test_sym(dcMatrix A)
{
	if (A.nbCols!=A.nbRows) return 0;
	unsigned int nCol=A.nbRows,i=0,j=0;
	
	if (nCol==1) return 1;
	
    
	for (i=1;i<nCol;i++)
        for (j=0;j<i;j++)
            if ( A(i,j)!=A(j,i) ) {break;}
	if ( (i==nCol) && (j==nCol-1) ) return 1;
	else return 0;
}

dcMatrix Id(unsigned int nCol)
{
    dcMatrix I(nCol);
    for (unsigned int i=0;i<nCol;i++) 
        for (unsigned int j=0;j<=i;j++)
        {
            if(i==j) I(i,j)=1;
            else { I(i,j)=0; I(j,i)=0; }
        } 
    return I;
}



dcMatrix power(dcMatrix A,unsigned int nCol)
{
    if (A.nbCols!=A.nbRows)
    {cout<<"Power over non square matrix impossible!"<<endl;
        exit(1);}
    else
    {
        dcMatrix B(A.nbCols);
        B=Id(A.nbCols);
        for(unsigned int i=0;i<nCol;i++) {B=B*A;}
        return B;
    }
}

dcMatrix cholesky(dcMatrix A)	//renvoi la dcMatrix triangul L tq:
{				//si A symetrique,carre
    if(A.nbCols!=A.nbRows) 		//L*transpo(L)=A
    { cout<<"Cholesky(non memb): non square dcMatrix !"<<endl;
        exit(1);}
    if(!test_sym(A)) 
    {  cout<<"Cholesky(non memb): non symetric dcMatrix !"<<endl;
        exit(1);
    }
	
    unsigned int nCol=A.nbCols;
    dcMatrix L(nCol);			
    float s;
    for(unsigned int i=0;i<nCol;i++)
    {
        s=0;
        for(unsigned int k=0;k<i;k++) s+=L(i,k)*L(i,k);
        L(i,i)=sqrt(A(i,i)-s);
        
        for(unsigned int j=i+1;j<nCol;j++) 
        {
            s=0;
            for(unsigned int k=0;k<i;k++) s+=L(j,k)*L(i,k);
            L(j,i)=(A(i,j)-s)/L(i,i); 
        }
    }
    return L;
} 


bool same_size(dcMatrix A, dcMatrix B)
{
	bool res = true;
	unsigned int ra = A.getNbRows();
	unsigned int rb = B.getNbRows();
	unsigned int ca = A.getNbCols();
	unsigned int cb = B.getNbCols();
	
	if (ra!=rb || ca!=cb) res = false;
	
	return res;
}

double distance_Matrix(dcMatrix A, dcMatrix B, double power)
{
	if (!same_size(A,B))
	{
		cout << "ERROR [distance_Matrix]: matrix not the same size!"<<endl;
		exit(1);
	}
	
	double dist = 0.0;
	
	for (unsigned int i=0; i<A.getNbRows(); i++)
	{
		for (unsigned int j=0; j<A.getNbCols(); j++) {
			dist += pow(A(i,j)-B(i,j),power);
		}
	}
	
	return pow(dist,1.0/power);
}


dcMatrix rowBind(dcMatrix A, dcMatrix B)
{
	if (A.getNbCols()!=B.getNbCols())
	{
		cout << "ERROR [rowBind]: matrix not the same column size!"<<endl;
		exit(1);
	}
	
	dcMatrix M(A.getNbRows()+B.getNbRows(),A.getNbCols());
	unsigned int i=0;
	
	for (i=0; i<A.getNbRows(); i++) {
		for (unsigned int j=0; j<A.getNbCols(); j++) {
			M(i,j) = A(i,j);
		}
	}
	for (unsigned int ii=0; ii<B.getNbRows(); ii++) {
		for (unsigned int j=0; j<A.getNbCols(); j++) {
			M(i+ii,j) = B(ii,j);
		}
	}
	
	return M;	
}



double Det(dcMatrix A)
{
    assert (A.nbCols==A.nbRows);
    double s=0;
    unsigned int nCol=A.nbCols;
    if (nCol==2) return A(0,0)*A(1,1)-A(0,1)*A(1,0);
    else
    {
        dcMatrix B(nCol-1);
        for(unsigned int k=0;k<nCol;k++)
        {
            for(unsigned int i=0;i<nCol;i++)
                for(unsigned int j=1;j<nCol;j++) 
                { 
                    if (i<k) B(i,j-1)=A(i,j);
                    if (i>k) B(i-1,j-1)=A(i,j);
                }
            s+=A(k,0)*Det(B)*pow (-1,k);
        } 
        return s;
    }
}

