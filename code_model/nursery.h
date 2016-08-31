//
//  nursery.h
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-12-09.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#ifndef __LocalSTI__nursery__
#define __LocalSTI__nursery__

#include "dcTools.h"
#include "individuals.h"


class Nursery
{
	/// Record of all births and children following
	/// all pregnancies in a 'Simulation'
	
	vector<Individual>		_child;
	
	vector<unsigned long>	_uid_mother;	// mother's UID
	vector<unsigned long>	_age_mother;	// mother's age at child birth
	
	vector<double>			_dob;			// date of birth
	
	vector<STI>				_STI;			// STI pattern from Population
	
public:
	
	Nursery(){}
	
	Nursery(vector<STI> s) {_STI = s;}
	
	
	void	set_STI(vector<STI> s) {_STI = s;}
	
	void	add_child(Individual child,
					  unsigned long uid_mother,
					  double age_mother,
					  double dob,
				   vector<bool> STI_MTCT);
	
	unsigned long	census_infected(STIname stiname);
	
	
	void	saveToCSVfile(string filename);
};

#endif /* defined(__LocalSTI__nursery__) */
