//
//  globalVar.h
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 2014-08-28.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#ifndef __STIagent_AIR__globalVar__
#define __STIagent_AIR__globalVar__

#include <iostream>
#include <random>


//  ==== Directories ====

extern std::string _DIR_IN;		// all input files are saved in this folder
extern std::string _DIR_OUT;	// all "*.out" files are saved in this folder
extern std::string _DIR_CALIB;	// all calibration files are saved in this folder




// ==== Random seed ====

extern unsigned int	_RANDOM_SEED;		// seed for random number generators
extern std::mt19937_64 _RANDOM_GENERATOR;


#endif /* defined(__STIagent_AIR__globalVar__) */
