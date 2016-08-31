/// =================================================================
///
///   WRAPPING FOR R
///
///   Created 2015-10-06 by David Champredon
///
/// =================================================================

#include <random>
#include <Rcpp.h>
using namespace Rcpp;

#include "globalVar.h"
#include "dcTools.h"
#include "dcMatrix.h"
#include "RV.h"
#include "population.h"
#include "simulation.h"
#include "MCsimulation.h"
#include "compare_simulation.h"


List dcDataFrameToRcppList(dcDataFrame df,
						   bool addRowName){
	
	/// Converts a dcDataFrame to a Rccp list
	/// (helper function)
	
	unsigned long ncol = df.get_colname().size();
	
	// Translate the 'dcDataFrame' into a R list
	// (convert to data frame in R, I don't know how to do it here in Rcpp):
	Rcpp::List rcpplist;
	// each column of the dcDataFrame is a list:
	for(int j=0; j<ncol; j++)
		rcpplist.push_back(df.get_value().extractColumn(j));
	// set the associated column names:
	if(!addRowName) rcpplist.attr("names") = df.get_colname();
	
	// insert row names (as an additional column, don't know how to do differently)
	if(addRowName) {
		rcpplist.push_back(df.get_rowname());
		vector<string> cn = df.get_colname();
		cn.push_back("rowname");
		rcpplist.attr("names") = cn;
	}
	
	return rcpplist;
}

// WARNING:
// commented line below (rcpp export) is necessary
// to export the function defined just after.

// [[Rcpp::export]]
List stiagent_runsim(List params) {
	
	/// Run ONE SINGLE Monte-Carlo iteration
	/// of the STIAGENT model,
	/// given inputs parameter values
	/// in both 'inputs' and 'calibration' folders
	
	
	// Unpack parameters:
	//
	string _DIR_IN		= params["folder_inputs"]; //../inputs/";
	string _DIR_CALIB	= params["folder_calib"];
	string scenario_file= params["scenario_file"];
	int displayProgress	= params["displayProgress"];
	
	unsigned long founder_size	= params["founder_size"];
	double founder_femprop		= params["founder_femprop"];
	double founder_cswprop		= params["founder_cswprop"];
	
	
	// ==============================================
	// === Important comment ===
	//
	// ID of the current MC.
	// Used when several MC specified _and_
	// several scenarios that must be compared
	// hence must use same seed for ith MC iteration.
	//
	// For example ('seed=mc_id'):
	// Scen 1, MC 1  ==> seed=1
	// Scen 1, MC 2  ==> seed=2
	// Scen 1, MC 3  ==> seed=3
	//
	// Scen 2, MC 1  ==> seed=1  (<--same seed as scen 1, MC 1)
	// Scen 2, MC 2  ==> seed=2
	// Scen 3, MC 3  ==> seed=3
	//
	// etc...
	//
	// Hence, it is *NOT* the nuber of Monte Carlo iterations!!!!
	// ==============================================
	unsigned int MC_id	= params["MC_id"];
	
	
	// ======================
	// === Initialization ===
	// ======================
	
	// Initialize empty Population object
	Population P(0);
	bool debugInfo=true;
	
	P.setup_for_simulation(founder_size,
						   founder_femprop,
						   founder_cswprop,
						   _DIR_IN,
						    "in_STI.csv",
						    "in_STI_SFincrease.csv",
						    "in_HIVrebound.csv",
						    "in_STItreatment.csv",
						    "in_STI_vaccine.csv",
						    "in_STI_genomes.csv",
						   debugInfo);
	
	// Detailed population:
	dcDataFrame pop_init	= P.export_to_dataframe();
	Rcpp::List pop_init_R	= dcDataFrameToRcppList(pop_init,false);
	
	
	// ======================
	// === Run simulation ===
	// ======================
	
	// retrieve parameters:
	double horizon			= getParameterFromFile("horizon_years", _DIR_IN + "in_simulation.csv");
	double timeStep			= getParameterFromFile("timestep_days", _DIR_IN + "in_simulation.csv")/365.0;
	double horizon_prtn		= getParameterFromFile("horizon_prtn_years", _DIR_IN + "in_simulation.csv");
	double timestep_prtn	= getParameterFromFile("timestep_prtn_days", _DIR_IN + "in_simulation.csv")/365.0;
	

	bool TraceNetwork		= false;
	
	string file_init_STI	= _DIR_IN + "in_STI_initial_prevalence.csv";
	
	// Scenario
	// (includes all interventions that will
	//  be done during this simulation)
	vector<string> file_intervention;
	string file_interv_base =_DIR_IN + scenario_file;
	vectorFromCSVfile_string(file_intervention,file_interv_base.c_str(), 1);
	file_intervention = trim(file_intervention);
	
	
	// Run the simulation:
	Simulation Sobj = runSimulation_one_obj(P,
											file_init_STI,
											file_intervention,
											horizon_prtn,
											timestep_prtn,
											horizon, timeStep,
											TraceNetwork,
											displayProgress,
											MC_id,
											_DIR_IN,
											_DIR_CALIB);
	
	// =============
	//    Outputs
	// =============

	Population POP = Sobj.get_population();
	
	unsigned long nSTI = POP.get_nSTImodelled();
	
	// data frame of time series:
	dcDataFrame df_sim		= Sobj.get_df_sim();
	Rcpp::List df_sim_R		= dcDataFrameToRcppList(df_sim,false);

	// Detailed population:
	dcDataFrame pop_last	= POP.export_to_dataframe();
	Rcpp::List pop_last_R	= dcDataFrameToRcppList(pop_last,false);

	// Detailed interventions:
	dcDataFrame df_interv	= Sobj.get_df_interv();
	Rcpp::List df_interv_R	= dcDataFrameToRcppList(df_interv,true);
	
	
	vector<double> cuminc_mtct_final;
	vector<string> stiname_str;
	for(unsigned long i=0; i<nSTI; i++){
		STIname stiname = POP.get_STI()[i].get_name();
		stiname_str.push_back(STInameString(stiname));
		cuminc_mtct_final.push_back(Sobj.get_nursery().census_infected(stiname));
	}
	

	// infectivity curves:
	Rcpp::List IC;
	vector<std::string> fm;
	fm.push_back("female"); fm.push_back("male");
	
	for(int i=0; i<nSTI;i++){
		STIname stiname = POP.get_STI()[i].get_name();

		Rcpp:List tmp_ic;
		tmp_ic.push_back(POP.get_infectivityCurve(stiname, female));
		tmp_ic.push_back(POP.get_infectivityCurve(stiname, male));
		tmp_ic.attr("names") = fm;
		IC.push_back(tmp_ic);
	}
	IC.attr("names") = stiname_str;
	
	// Genomes
	Rcpp::List genomes;
	
	// TO DO: make a loop and retrieveSTI name...
	vector<string> g_HIV = Sobj.get_population().get_population_genomes(0);
	vector<string> g_Tp  = Sobj.get_population().get_population_genomes(1);
	genomes.push_back(g_HIV);
	genomes.push_back(g_Tp);
	
	// =========================================================================
	// =========================================================================
	//
	// R List storing all outputs:
	//
	return List::create(Named("df_sim") = df_sim_R,
						Named("cuminc_mtct_final")= cuminc_mtct_final,
						Named("popsize_alive") = POP.census_alive(),
						Named("STInames") = stiname_str,
						Named("infCurves") = IC,
						Named("population_inital") = pop_init_R,
						Named("population") = pop_last_R,
						Named("df_interv") = df_interv_R,
						Named("seed") = MC_id,
						Named("genomes") = genomes
						);
}




// [[Rcpp::export]]
List stiagent_comp_interv(List params) {
	
	///  ---- WARNING: MAYBE BE OBSOLETE (2015-11-11) ----
	///  ---- TO DO: check if really needed now that new code does better...
	
	///
	/// Run several scenarios and compare their outcomes.
	/// Returns an R list with size the number of STI modeled.
	/// Each element of the list is itself a list storing the
	/// outcome values for all scenarios.
	///
	/// For example, if HIV and Tp are the two STIs modeled and
	/// 3 scenarios are defined,the output list will look like
	/// this in R:
	///
	/// x <- stiagent_comp_interv(params)  # <== call this function in R
	///
	/// then structure of 'x' is:
	///
	/// x
	///  |_ HIV
	///  |     |_ prev_final0 : 0.031
	///  |     |_ prev_final1 : 0.036
	///  |     |_ prev_final2 : 0.023
	///  |     |_ cuminc_final0 : 0.12
	///  |     |_ cuminc_final1 : 0.00
	///  |     |_ cuminc_final2 : 0.04
	///  |     |_ ...etc...
	///  |
	///  |_ Tp
	///        |_ prev_final0 : 0.22
	///        |_ prev_final1 : 0.07
	///        |_ prev_final2 : 0.0
	///        |_ cuminc_final0 : 0.87
	///        |_ cuminc_final1 : 0.23
	///        |_ cuminc_final2 : 0.11
	///        |_ ...etc...
	///
	
	
	// Folders where parameters are saved:
	string _DIR_IN		= params["folder_inputs"]; //../inputs/";
	string _DIR_CALIB	= params["folder_calib"];
	int jobnum			= params["jobnum"];
	
	unsigned long founder_size	= params["founder_size"];
	double founder_femprop		= params["founder_femprop"];
	double founder_cswprop		= params["founder_cswprop"];
	
	// ======================
	// === Initialization ===
	// ======================
	
	// Initialize empty Population object
	Population P(0);
	
	bool debugInfo=true;
	
	P.setup_for_simulation(founder_size,
						   founder_femprop,
						   founder_cswprop,
						   _DIR_IN,
						   "in_STI.csv",
						   "in_STI_SFincrease.csv",
						   "in_HIVrebound.csv",
						   "in_STItreatment.csv",
						   "in_STI_vaccine.csv",
						   "in_STI_genomes.csv",
						   debugInfo);
	
	
	// Scenario files files
	string path_scenario = _DIR_IN + "in_scenario.csv";
	vector<string> file_scenario;
	vectorFromCSVfile_string(file_scenario, path_scenario.c_str(), 1);
	
	//	cout<<endl<<"Scenario files:";
	//	displayVector(file_scenario);
	
	// ======================
	// === Run simulation ===
	// ======================
	
	// retrieve parameters:
	double horizon			= getParameterFromFile("horizon_years", _DIR_IN + "in_simulation.csv");
	double timeStep			= getParameterFromFile("timestep_days", _DIR_IN + "in_simulation.csv")/365.0;
	double horizon_prtn		= getParameterFromFile("horizon_prtn_years", _DIR_IN + "in_simulation.csv");
	double timestep_prtn	= getParameterFromFile("timestep_prtn_days", _DIR_IN + "in_simulation.csv")/365.0;
	unsigned int iter_mc	= getParameterFromFile("MCiter", _DIR_IN + "in_simulation.csv");
	bool TraceNetwork		= false;
	
	int displayProgress		= 0;
	
	string file_init_STI	= _DIR_IN + "in_STI_initial_prevalence.csv";
	
	// DELETE?
	//vector<string> file_intervention;
	//string file_interv_base =_DIR_IN + "in_interv_baseline_wrapper.csv";
	//vectorFromCSVfile_string(file_intervention,file_interv_base.c_str(), 1);
	//file_intervention = trim(file_intervention);
	
	vector<dcDataFrame> df_comp_interv;
	
	df_comp_interv = run_comp_interventions_obj(iter_mc,
												P,
												file_init_STI,
												file_scenario,
												horizon_prtn,
												timestep_prtn,
												horizon,
												timeStep,
												TraceNetwork,
												displayProgress,
												jobnum,
												_DIR_IN,
												_DIR_CALIB
												);
	
	
	std::vector<string> list_name;
	unsigned long nsti = df_comp_interv.size();
	
	for(int i=0; i<nsti; i++)
		list_name.push_back(STInameString(P.get_STI()[i].get_name()));
		
	// Returns a list of lists:
	Rcpp::List x;
	
	// The response variable for each scenario:
	for(int i=0; i<nsti; i++)
		x.push_back(dcDataFrameToRcppList(df_comp_interv[i],false));
	
	// The name of each scenario:
	x.push_back(trim(file_scenario));
	list_name.push_back("scenarioName");
	
	// Name each list element:
	x.attr("names") = list_name;
	
	return(x);
}




















