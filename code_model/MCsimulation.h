//
//  MCsimulation.h
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-11-17.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#ifndef __LocalSTI__MCsimulation__
#define __LocalSTI__MCsimulation__

//#include <stdio.h>
#include "simulation.h"


class MCsimulation
{
	/// Object storing all Monte Carlo trials of a simulation
	
	
	unsigned long		_nMC;	// total number of MC iterations
	vector<Simulation>	_simulation;
	
	double			_horizon;
	double			_timeStep;
	vector<double>	_schedule;
	
	string			_scenarioName;
	
public:
	
	// ============================
	// ==== CONSTRUCTORS  =========
	// ============================
	
	MCsimulation(){}
	MCsimulation(vector<Simulation> S);
	

		
	// ============================
	// ==== GET FUNCTIONS =========
	// ============================
	
	vector<Simulation>		get_simulation(){return _simulation;}
	Simulation				get_simulation(int i){return _simulation[i];}
	unsigned long			get_nMC() {return _nMC;}
	string					get_scenarioName() {return _scenarioName;}
	
	
	// ============================
	// ==== SET FUNCTIONS =========
	// ============================
	
	void	set_simulation(vector<Simulation> x){_simulation = x;}
	void	set_nMC(unsigned long n) {_nMC=n;}
	void	set_scenarioName(string x) {_scenarioName = x;}
	void	set_horizon(double x){_horizon=x;}
	void	set_timeStep(double x){_timeStep=x;}
	void	set_schedule(vector<double> x) {_schedule=x;}
	
	// ============================
	// ==== MISC   =========
	// ============================
	
	vector<double> population_final();
	
	// ============================
	// ==== STI RELATED   =========
	// ============================
	
	vector<double>	prevalence_final(STIname);
	vector<double>	cumul_incidence_final(STIname);
	vector<double>	cumul_incidence_MTCT_final(STIname);
	
	vector<double>	mean_prevalence_ts(STIname);
	double			mean_prevalence_final(STIname);
	
	vector<double>	mean_incidence(STIname);
	
	double			mean_cumul_incidence_MTCT(STIname);
	
	
	// ============================
	// ==== CALIBRATION   =========
	// ============================
	
	vector<double>	distance_from_targets();
	
	
	// ============================
	// ====	   UTILS      =========
	// ============================
	
	void displayInfo();
	
	// Not yet implemented -- see and copy from  'Simulation.h'
	void save_outputs_demog(string pathFolder); // Save to file demographic outputs (e.g. age distribution) ; formatted for comparison w/ calibration targets
	void save_outputs_prtnr(string pathFolder);	// Save to file partnership outputs (e.g. age gap distribution, single ratio,...)
	void save_outputs_sex(string pathFolder);	// Save to file sexual behavior outputs (e.g. number of lifetime sex partners distribution)
	void save_outputs_epi(string pathFolder);	// Save to file epidemiological outputs (e.g. STI prevalences)
	
};


// One stochastic run of the whole simulation
// (used for R library wrap)
Simulation	runSimulation_one_obj(Population P_init,
								  string filename_init_STI_prev,
								  vector<string> filename_interventions,
								  double horizon_prtn,
								  double timestep_prtn,
								  double horizon,
								  double timestep,
								  bool TraceNetwork,
								  int displayProgress,
								  unsigned int mc_id,
								  string folder_inputs,
								  string folder_calib);



vector<Simulation>	runSimulationMC_obj(unsigned int nMC,
										Population P_init,
										string filename_init_STI_prev,
										vector<string> filename_interventions,
										double horizon_prtn,
										double timestep_prtn,
										double horizon,
										double timestep,
										bool TraceNetwork,
										int displayProgress,
										string folder_inputs,
										string folder_calib,
										int jobnum=1);




// Calculates _average_ distance from calibration targets
double	average_distance_target(vector<Simulation>);


#endif /* defined(__LocalSTI__MCsimulation__) */
