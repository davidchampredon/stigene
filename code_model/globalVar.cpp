//
//  globalVar.cpp
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 2014-08-28.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#include "globalVar.h"


//  ==== Directories ====


std::string _DIR_OUT		= "../outputs/";
std::string _DIR_CALIB		= "../calibration/";


// ==== Random seeds ====

unsigned int	_RANDOM_SEED	= 12345;
std::mt19937_64	_RANDOM_GENERATOR(_RANDOM_SEED);