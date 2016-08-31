/*
 *  RV.cpp
 *
 *
 *  Created by David Champredon  on 13-04-19.
 *  Copyright 2013-14. All rights reserved.
 *
 */

#include "RV.h"
#include "dcTools.h"
#include "globalVar.h"

using namespace std;


void force_seed_reset()
{
	_RANDOM_GENERATOR.seed(_RANDOM_SEED);
}

void force_seed_reset(unsigned int manual_seed)
{
	// DEBUG
	// cout << "SEED IS SET TO: "<<manual_seed<<endl;
	// -----
	_RANDOM_GENERATOR.seed(manual_seed);
}



/* ****************************************************
 CONTINUOUS SUPPORT
 ****************************************************/


double uniform01()
{
	uniform_real_distribution<double> dist(0,1);
	return dist(_RANDOM_GENERATOR);
}


double expo(double lambda)
{
	exponential_distribution<double> dist(lambda);
	return dist(_RANDOM_GENERATOR);
	
}


double normal(double mean, double sd)
{
	normal_distribution<double> dist(mean, sd);
	return dist(_RANDOM_GENERATOR);
}

double gamma(double shape, double scale)
{
	gamma_distribution<double> dist(shape,scale);
	return dist(_RANDOM_GENERATOR);
}

double beta(double a, double b)
{
	gamma_distribution<double> X(a,1.0);
	gamma_distribution<double> Y(b,1.0);
	
	double xx = X(_RANDOM_GENERATOR);
	double yy = Y(_RANDOM_GENERATOR);
	
	return xx/(xx+yy);
}



/* ****************************************************
			DISCRETE SUPPORT
 ****************************************************/


int	uniformInt(int nmin, int nmax)
{
	uniform_int_distribution<> dist(nmin,nmax);
	return dist(_RANDOM_GENERATOR);
}

int geometric(double p)
{
	geometric_distribution<int> dist(p);
	return dist(_RANDOM_GENERATOR);
}

//unsigned long	binom(double p, unsigned long N)
//{
//	stopif( (p<0 || p>1),
//		   "probability ("+to_string(p)+") not in [0;1]");
//
//	binomial_distribution<unsigned long> dist(N,p);
//	unsigned long res = dist(_RANDOM_GENERATOR);
//    return res;
//}

// WARNING: PROBLEM WHEN UNSIGNED LONG IN OLDER VERSION OF COMPILER
int	binom(double p, int N)
{
    stopif( (p<0 || p>1),
           "probability ("+to_string(p)+") not in [0;1]");
    
    binomial_distribution<int> dist(N,p);
    int res = dist(_RANDOM_GENERATOR);
    return res;
}




int poisson(double expectedValue)
{
	poisson_distribution<int> dist(expectedValue);
	return dist(_RANDOM_GENERATOR);
}



vector<unsigned int> multinomial(unsigned int N, vector<double> proba)
{
	// Retrieve the dimension that defines the random variables vector
	unsigned long dim = proba.size();
	
	// DEBUG
	//cout << "binom_dim:" <<dim << endl;
	// cout << "binom_trials:" <<N << endl;
	
	// Healthy checks
	if (dim==0 || N==0){
		cout << endl << "ERROR [multinomial]: can't generate with dimension=0"<<endl;
		if (dim==0) cout << "size of proba vector = 0 !"<<endl;
		if (N==0) cout << "N = 0 !"<<endl;
		exit(1);
	}
	double s=0;
	for (int i=0; i<proba.size(); i++) {s+=proba[i];}
	if (fabs(s-1)>0.0001){
		cout << endl << "ERROR [multinomial]: probabilities do not add-up to 1"<<endl;
		displayVector(proba);
		exit(1);
	}
	
	// Partition interval [0;1] with respect to probabilities
	vector<double> a(dim+1);
	a[0]=0;
	for (int i=1; i<dim+1; i++) a[i] = a[i-1]+proba[i-1];
	
	// initiate the result vector
	vector<unsigned int> res(dim,0); // all values at 0
	
	// Draw the N trials
	for (int k=0; k<N; k++) {
		double u = uniform01();
		// Find the interval where "u" is.
		int m;
		for (m=1; a[m]<u; m++) {}
		
		// "u" is between a[m-1] and a[m]
		// so the (m-1)th value of the stochastic vector is increased
		
		res[m-1]++;
	}
	return res;
}



int	probaHistogramInt(vector<int> values, vector<double> probas)
{
	/// Returns (randomly) one of the element of 'values' based
	/// on its associated probability (defined by 'probas')
	
	// --- Integrity checks ---
	string errmsg ="size of 'values' (" + to_string(values.size()) + ") and 'probas' (" + to_string(values.size()) +") not the same!";
	stopif(values.size() != probas.size(),errmsg);
	double se = sumElements(probas);
	stopif(fabs(se-1.00)>1e-7,"Probabilities do not add up to 1!");
	
	// ----------------------
	
	double s	= probas[0];
	int k		= 0;
	int kk		= 0;
	bool found	= false;
	double u	= uniform01();
	
	while(k<probas.size() && !found){
		if (u<s) {
			kk=k;
			found = true;
		}
		k++;
		s += probas[k];
	}
	return values[kk];
}
	

