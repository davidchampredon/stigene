//
//  nursery.cpp
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-12-09.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#include "nursery.h"


void Nursery::add_child(Individual child,
						unsigned long uid_mother,
						double age_mother,
						double dob,
						vector<bool> STI_MTCT)
{
	/// Add a child to the Nursery following birth from pregnancy.
	/// Includeds STI mother-to-child transmission
	
	// Set STI status _AT BIRTH_
	// (in Nursery, STI duration meaningless, just STI status at birth matters)
	child.set_STI(_STI);
	
	for (int i=0; i<STI_MTCT.size(); i++)
	{
		if(STI_MTCT[i])
			child.set_STIduration(_STI[i].get_name(), 1.0/365.0);
	}
	_child.push_back(child);
	
	// Add other informations:
	_uid_mother.push_back(uid_mother);
	_age_mother.push_back(age_mother);
	_dob.push_back(dob);
}



void Nursery::saveToCSVfile(string filename)
{
	ofstream f(filename.c_str());
	
	// headers
	
	f << "nursery_id,uid_mother,age_mother,dob,";
	
	for (int d=0; d<_STI.size(); d++){
		f<<STInameString(_STI[d].get_name());
		if(d<_STI.size()-1) f<<",";
	}
	f << endl;
	
	// values
	
	for (int i=0; i<_child.size(); i++){
		f << i+1 << ",";
		f << _uid_mother[i] << ",";
		f << _age_mother[i] << ",";
		f << _dob[i] << ",";
		
		// STI information
		for (int d=0; d<_STI.size(); d++){
			f << _child[i].get_STIduration()[d];
			if(d<_STI.size()-1) f<<",";
		}
		f << endl;
	}
	
}



unsigned long Nursery::census_infected(STIname stiname)
{
	/// Counts the number of children infected (by their mother)
	/// with a given STI
	
	int pos_sti = positionSTIinVector(stiname, _STI);
	unsigned long cnt = 0;
	
	for (int i=0; i<_child.size(); i++)
		if (_child[i].get_STIduration()[pos_sti]>0) cnt++;
	
	return cnt;
}














