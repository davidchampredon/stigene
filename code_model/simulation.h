/*
 *  simulation.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-08-30.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef simulation_h
#define simulation_h


#include "population.h"
#include "dcDataFrame.h"
#include "nursery.h"
#include "intervention.h"

class Simulation 
{
	
	// Monte-Carlo parameters
	double			_horizon;
	double			_timeStep;
	double			_simulationTime;
	vector<double>	_schedule;
	
	unsigned int	_MC_trial_iter;		// number of this Monte-Carlo trial (out of the total number of trials)
	
	// Population
	Population		_population_initial;	// store the initial population (used in calibration)
	Population		_population;
	
	// New borns
	unsigned long	_newborn_timestep;		// number of newborns (during a Simulation timestep)
	Nursery			_nursery;

	// STI related
	dcMatrix		_STI_incidence;		// rows=time ; columns = STI
	dcMatrix		_STI_prevalence;	// rows=time ; columns = STI
	
	// Interventions (Treatments, Vaccinations)
	vector<Intervention>	_intervention;

	// Data frames that record outputs
	dcDataFrame		_df_sim;		// time series
	dcDataFrame		_df_interv;		// intervention details
	
	
	// ===============================
	// === TARGETS FOR CALIBRATION ===
	// ===============================
	
	int		_nMC_calibration;	// Number of MC runs for the simulation; used for calibration purposes
	
	// Vector of times during the simulation
	// when a calibration is expected
	vector<double>				_calibrationTime;
	
	
	// For each calibration time, there can be
	// an arbitrary number of output to calibrate
	
	// filename of the associated target
	vector< vector<string> >	_calibrationTargetFile;
	// type of target (e.g., prevalence, age distribution, etc)
	vector< vector<string> >	_calibrationTargetType;

	// Relative weight in the glabal calibration score of each calibration target
	vector< vector<double> >	_calibrationWeights;
	
	// Weighted Distances model from target (needs to be calculated)
	vector< vector<double> >	_calibrationDistances;
	
	
	bool			_save_trace_files;
	
public:
	
	// === CONSTRUCTOR ===
	
	Simulation (){}
	
	Simulation (double horizon, double timeStep, Population P, 
				int nMC_calibration);
	
	void		set_population(Population P) {_population=P;}
	
		
	// ==== SET FUNCTIONS ====
	
	void		set_horizon(double x) {_horizon = x;}
	void		set_timeStep(double x) {_timeStep = x;}
	void		set_MC_trial_iter(unsigned int i){_MC_trial_iter=i;}
	
	void		set_newborn_timestep(unsigned long n) {_newborn_timestep=n;}
	
	void		set_intervention(vector<Intervention> v) {_intervention = v;}
	
	void		set_calibrationTime(vector<double> x) {_calibrationTime=x;}
	void		set_calibrationOutput(int i, vector<string> s) {_calibrationTargetFile[i] = s;}
	
	void		set_save_trace_files(bool x) {_save_trace_files = x;}
	
	// =======================
	// =======================
	// ==== GET FUNCTIONS ====
	// =======================
	// =======================
	
	Population			get_population() {return _population;}
	Nursery				get_nursery() {return _nursery;}
	
	unsigned int		get_MC_trial_iter() {return _MC_trial_iter;}
	
	double				get_horizon(){return _horizon;}
	double				get_timeStep(){return _timeStep;}
	vector<double>		get_schedule(){return _schedule;}
	
	dcMatrix			get_STI_incidence(){return _STI_incidence;}
	dcMatrix			get_STI_prevalence(){return _STI_prevalence;}
	
	
	vector<double>			get_calibrationTime() {return _calibrationTime;}
	vector<vector<string> >	get_calibrationTargetFile() {return _calibrationTargetFile;}
	vector<vector<double> >	get_calibrationDistances() {return _calibrationDistances;}
	
	
	// FOR ALL AGE (AGE GAP) DISTRIBUTIONS:
	// 1st col = age breaks
	// 2nd col = proportion of population
	// in associated age range age[i]<prop[i]<age[i+1]
	
	
	dcDataFrame			get_df_sim(){return _df_sim;}
	dcDataFrame			get_df_interv(){return _df_interv;}
	
	// =======================
	// =======================
	// ==   UPDATE EVENTS   ==
	// =======================
	// =======================
	
	void			update_pregnancies(double timestep);
	void			set_pregnant(unsigned long uid);
	void			set_gestationDuration(unsigned long uid,double timestep);
	void			increment_gestationDuration(unsigned long uid,double timestep);
	void			increment_nChildBorn(unsigned long uid);

	
	// =======================
	// =======================
	// ====      STI      ====
	// =======================
	// =======================
	
	void			STI_set_initial_prevalence(string filename);
	vector<bool>	MTCT(unsigned long uid);
	double			STI_cumul_incidence(STIname s,int time_i);
	
	
	
	// =======================
	// =======================
	// ====  INTERVENTION ====
	// =======================
	// =======================
	
	void	activate_intervention(int i);
	
	
	// =======================
	// =======================
	// ====  TREATMENT    ====
	// =======================
	// =======================
	
	void	treat_indiv(unsigned long uid, STIname sti);
//	void	treat_all(STIname sti);
//	void	treat_any_proportion(STIname sti, double proportion);
//	void	treat_symptomatic(STIname sti);
	
	void	cure_indiv(unsigned long uid, STIname sti);
	void	update_cure(STIname sti);
	
	
	// =======================
	// =======================
	// ====  VACCINATION  ====
	// =======================
	// =======================
	
	void	vaccinate_indiv(unsigned long uid, STIname stiname);
	void	update_vacc(STIname stiname);

//	void	vaccinate_all(STIname stiname);
//
//	void	vaccinate_general(STIname stiname, double propPerYear);
	
	
	// =======================
	// =======================
	// ==== RUN SIMULATION ===
	// =======================
	// =======================
	
	// Run all events (births, partnerships, sex acts, deaths,...)
	// for a given time step
	
	void		runAllEvents_timeStep_obj(int numTimeStep,
										  bool doSex,
										  bool logIndivInfo);
	
	// Run all events until '_horizon'
	// returns output to be calibrated
	// One stochastic realization only
	
	void		runAllEvents_horizon_obj(bool doSex,
										 bool logIndivInfo,
										 bool traceNetwork,
										 int displayProgress);

	
	// =====================
	// ==== CALIBRATION ====
	// =====================

	void		set_calibration_schedule(string filename_calibrationTime,
										 string filename_all_calib_output,
										 string filename_weight_target);
	
	bool			isCalibrationTime(double t);
	unsigned int	whichCalibrationTime(double t);

	

	
	// ==== Other =====
	
	void		displayInfo();

	
};

string calibration_file_to_type(string filename);



#endif