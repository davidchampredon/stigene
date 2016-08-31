//
//  intervention.cpp
//  TpHIV
//
//  Created by David CHAMPREDON on 2015-01-05.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "intervention.h"


Intervention::Intervention(string file)
{
	/// CONSTRUCT AN Intervention OBJECT FROM A FILE
	
	_type = getParameterFromFile_string("_type", file);
	
	string stiname = getParameterFromFile_string("_stiname", file);
	_stiname = StringToSTIname(stiname);
	
	_schedule.push_back(getParameterFromFile("_schedule_0", file));
	_schedule.push_back(getParameterFromFile("_schedule_1", file));
	// TO DO: implement to manage more than 2 dates
	
	_annCvgRate = getParameterFromFile("_annCvgRate", file);
	
	
	// -- Integrity checks --
	
	stopif((_schedule.size()%2) != 0,
		   "Intervention schedule must be defined with pairs of start and end dates!");
	
	for (int i=0; i<_schedule.size()-1; i++)
	{
		stopif(_schedule[i]>=_schedule[i+1],
			   "Intervention schedule dates not ordered!");
	}
	
}



void Intervention::displayInfo()
{
	coutline(80);
	cout <<endl << " ===== INTERVENTION FEATURES ======"<<endl;
	cout <<endl;
	cout << "_stiname:"<< STInameString(_stiname) <<endl;
	cout << "_type: "<<_type<<endl;
	cout << "_schedule_0: "<<_schedule[0]<<endl;
	cout << "_schedule_1: "<<_schedule[1]<<endl;
	cout << "_annCvgRate: "<<_annCvgRate<<endl;
	
	coutline(80);
}