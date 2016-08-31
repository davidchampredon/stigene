//
//  intervention.h
//  TpHIV
//
//  Created by David CHAMPREDON on 2015-01-05.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __TpHIV__intervention__
#define __TpHIV__intervention__

#include <stdio.h>
#include "dcTools.h"
#include "STI.h"

#endif /* defined(__TpHIV__intervention__) */



class Intervention{
	
	
	// Type of intervention (e.g., syndromic treatment, mass vaccination, ...)
	string	_type;
	
	// Which STI the intervention targets
	STIname _stiname;
	
	
	// Schedule of the intervention.
	// MUST be pair of dates representing start and end date
	// of each intervention period.
	//
	// For example:
	// _schedule[0]: start of the 1st intervention
	// _schedule[1]: end of the 1st intervention
	// _schedule[2]: start of the 2nd intervention
	// _schedule[3]: end of the 2nd intervention
	
	vector<double> _schedule;
	
	
	// Annual coverage rate of targeted population for intervention.
	//
	// ** Warning **
	// Although defined as an annual rate, this will be applied at each time step.
	// For example, if annual coverage is 50% and time step is 0.1, then
	// for a given time step, only 0.5*0.1=0.05=5% of the targeted population
	// will receive intervention.
	
	double			_annCvgRate;
	vector<double>	_annCvgRate_multi;
	
	
public:
	
	// ===== (PSEUDO) CONSTRUCTOR =====
	
	Intervention(){}
	
	Intervention(string file);
	
	
	// ===== SET FUNCTIONS =====
	
	void set_type(string s){_type = s;}
	void set_schedule(vector<double> x){_schedule = x;}
	
	void set_stiname(STIname s) {_stiname = s;}
	
	void set_annCvgRate(double x){_annCvgRate = x;}
	
	void set_annCvgRate_multi(vector<double> x){_annCvgRate_multi = x;}
	
	
	
	// ===== GET FUNCTIONS =====
	
	string				get_type(){return _type;}
	
	vector<double>		get_schedule(){return _schedule;}
	
	STIname				get_stiname(){return _stiname;}
	
	double				get_annCvgRate(){return _annCvgRate;}

	
	vector<double> get_annCvgRate_multi(){return _annCvgRate_multi;}

	
	// ===== UTILS =====
	
	void	displayInfo();
	
};

