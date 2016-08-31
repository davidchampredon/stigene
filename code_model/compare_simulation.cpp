//
//  compare_simulation.cpp
//  TpHIV
//
//  Created by David CHAMPREDON on 2015-01-08.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "compare_simulation.h"
#include "MCsimulation.h"
#include "dcTools.h"




double comp_simul_mean_MTCT(STIname stiname,
							MCsimulation S1,
							MCsimulation S2)
{
	/// Compares b/w 2 simulations
	/// incidence of mother-to-child transmission
	/// of a given STI.
	/// Returns mean(incidence_1 - incidence_2)
	
	
	double inc1 = S1.mean_cumul_incidence_MTCT(stiname);
	double inc2 = S2.mean_cumul_incidence_MTCT(stiname);
	
	return inc1-inc2;
}



double comp_simul_mean_cumul_incidence(STIname stiname,
									   MCsimulation S1,
									   MCsimulation S2)
{
	/// Compares b/w 2 simulations
	/// cumulative incidence at simulation horizon
	/// of a given STI.
	/// Returns mean(incidence_1 - incidence_2)
	
	vector<double> inc1 = S1.mean_incidence(stiname);
	vector<double> inc2 = S2.mean_incidence(stiname);
	
	return inc1[inc1.size()-1] - inc2[inc2.size()-1];
}



double	comp_simul_mean_prevalence(STIname stiname,
								   MCsimulation S1,
								   MCsimulation S2)
{
	/// Compares b/w 2 simulations
	/// prevalence at simulation horizon
	/// of a given STI.
	/// Returns mean(prev_1 - prev_2)
	
	double prev1 = S1.mean_prevalence_final(stiname);
	double prev2 = S2.mean_prevalence_final(stiname);
	
	return prev1-prev2;
}








dcDataFrame results_scenario_outcome(STIname stiname,
								  vector<MCsimulation> S)
{
	/// data frame with results of the comparison of outcomes
	/// b/w different interventions scenario
	/// for a given STI
	
	stopif(S.size()==0,"Vector of simulations empty!");
	vector<string> colname;
	
	// ==== Prevalence at horizon
	
	vector<double> prev0 = S[0].prevalence_final(stiname);
	colname.push_back("prev_final0");
	
	dcMatrix M(prev0);
	
	for (int i=1; i<S.size(); i++){
		vector<double> prev_i = S[i].prevalence_final(stiname);
		M.addColVector(prev_i);
		prev_i.clear();
		colname.push_back("prev_final"+int2string(i));
	}
	
	// ==== Cumulative incidence at horizon
	
	for (int i=0; i<S.size(); i++){
		vector<double> cum_inc_i = S[i].cumul_incidence_final(stiname);
		M.addColVector(cum_inc_i);
		cum_inc_i.clear();
		colname.push_back("cum_inc_final"+int2string(i));
	}
	
	// ==== Cumulative MTCT at horizon
	
	for (int i=0; i<S.size(); i++)	{
		vector<double> cum_inc_mtct_i = S[i].cumul_incidence_MTCT_final(stiname);
		M.addColVector(cum_inc_mtct_i);
		cum_inc_mtct_i.clear();
		colname.push_back("cum_inc_mtct_final"+int2string(i));
	}
	
	// ==== Population at horizon
	
	for (int i=0; i<S.size(); i++){
		vector<double> pop_final_i = S[i].population_final();
		M.addColVector(pop_final_i);
		pop_final_i.clear();
		colname.push_back("pop_final"+int2string(i));
	}
	
	// === Final data frame
	dcDataFrame D(M);
	D.set_colname(colname);
	
	return D;
}





void writeToFile_scenario_outcome(string filename,
								  STIname stiname,
								  vector<MCsimulation> S)
{
	/// Writes to a file the comparison of outcomes
	/// b/w different interventions scenario
	/// for a given STI
	
	dcDataFrame D = results_scenario_outcome(stiname,S);
	D.saveToCSV(filename, true);
	
	// Save scenario names to a file:
	string fname = _DIR_OUT + "scenario_names.out";
	ofstream g(fname.c_str());
	for (int i=0; i<S.size(); i++)
		g << i << "," << S[i].get_scenarioName() <<endl;
}





MCsimulation run_wrap_obj(unsigned int nMC,
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
						  int jobnum){
	
	vector<Simulation> S = runSimulationMC_obj(nMC,
											   P_init,
											   filename_init_STI_prev,
											   filename_interventions,
											   horizon_prtn,
											   timestep_prtn,
											   horizon,
											   timestep,
											   TraceNetwork,
											   displayProgress,
											   folder_inputs,
											   folder_calib,
											   jobnum);
	
	MCsimulation MCS;
	MCS.set_nMC(S.size());
	stopif(S.size()==0,"Simulation vector is empty!");
	MCS.set_simulation(S);
	
	MCS.set_horizon(S[0].get_horizon());
	MCS.set_timeStep(S[0].get_timeStep());
	MCS.set_schedule(S[0].get_schedule());
	
	return(MCS);
}


vector<dcDataFrame>	run_comp_interventions_obj(unsigned int nMC,
											   Population P_init,
											   string filename_init_STI_prev,
											   vector<string> scenario_list,
											   double horizon_prtn,
											   double timestep_prtn,
											   double horizon,
											   double timestep,
											   bool TraceNetwork,
											   int displayProgress,
											   int jobnum,
											   string folder_inputs,
											   string folder_calib)
{
	/// Run several scenarios and compare their outcomes
	/// Returns a data frame of scenario results for each STI
	
	// deals with eol character issues
	scenario_list = trim(scenario_list);
	unsigned long n_scenario = scenario_list.size();
	
	// Display basic info:
	coutline(100);
	cout << endl << " Inputs folder: "<< folder_inputs << endl;
	cout << endl << " Calibration folder: "<< folder_calib << endl;
	cout << endl << " Scenario list: "<<endl;
	displayVector(scenario_list);
	cout << endl << " Initial STI prevalence set from file: "<<filename_init_STI_prev<<endl;
	coutline(100);
	// ------------------
	
	vector<MCsimulation> scenario;
	
	for (int i=0; i<n_scenario; i++)
	{
		// retrieve relevant intervention wrapper
		vector<string> scenario_i;
		string file_i = folder_inputs + scenario_list[i];
		vectorFromCSVfile_string(scenario_i, file_i.c_str(), 1);
		
		// run the MC simulations
		MCsimulation tmp;
		
		tmp = run_wrap_obj(nMC,
						   P_init,
						   filename_init_STI_prev,
						   scenario_i,
						   horizon_prtn,
						   timestep_prtn,
						   horizon,
						   timestep,
						   TraceNetwork,
						   displayProgress,
						   folder_inputs,
						   folder_calib,
						   jobnum);
		
		tmp.set_scenarioName(scenario_list[i].c_str());
		scenario.push_back(tmp);
	}
	
	
	// ============
	// OUTPUT
	// ============
	
	vector<dcDataFrame> D;
	
	Population pop0 = scenario[0].get_simulation(0).get_population();
	unsigned long n_sti = pop0.get_nSTImodelled();
	
	for(int s=0; s<n_sti; s++){
		STIname stiname = pop0.get_STI()[s].get_name();
		D.push_back(results_scenario_outcome(stiname, scenario));
	}
	
	return(D);
}

