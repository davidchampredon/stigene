/*
 *  simulation.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-08-30.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "simulation.h"
#include "globalVar.h"



Simulation::Simulation(double horizon,
					   double timeStep,
					   Population P,
					   int nMC_calibration)
{
	_horizon = horizon;
	_timeStep = timeStep;
	
	_schedule.clear();
	for (double t=0; t<_horizon; t+=_timeStep)
		_schedule.push_back(t);
	
	_population = P;
	_population_initial = P;
	
	_population_initial.STI_update_templateToAllIndividuals();
	_population.STI_update_templateToAllIndividuals();
	
	_nMC_calibration = nMC_calibration;
	
	_newborn_timestep = 0;
	
	_STI_incidence.resize(1,_population.get_STI().size());
	_STI_prevalence.resize(1,_population.get_STI().size());
	
	_save_trace_files = false;
}


// ===========================================================
// ===========================================================

void Simulation::STI_set_initial_prevalence(string filename){
	_population.STI_set_initial_prevalence(filename);
}

void Simulation::set_pregnant(unsigned long uid){
	_population.set_pregnant(uid);
}

void Simulation::increment_gestationDuration(unsigned long uid, double x){
	_population.increment_gestationDuration(uid, x);
}

void Simulation::set_gestationDuration(unsigned long uid, double x){
	_population.set_gestationDuration(uid, x);
}

void Simulation::increment_nChildBorn(unsigned long uid){
	_population.increment_nChildBorn(uid);
}


vector<bool> Simulation::MTCT(unsigned long uid)
{
	/// Returns a vector of STI transmission to child
	/// (vector size = number of STIs modelled)
	/// ('true' = transmission to child)
	
	Individual I = _population.get_individual(uid);
	
	// MTCT at previous time step
	vector<bool> prev_mtct = I.get_STI_MTCT();
	
	// Loop through all STIs
	unsigned long nsti = _population.get_STI().size();
	vector<bool> mtct = prev_mtct;
	
	// DEBUG
	//displayVector(mtct);
	
	for (int i=0; i<nsti; i++)
	{
		STIname stiname = _population.get_STI()[i].get_name();
		
		if(!prev_mtct[i]){
			double stiduration = I.get_STIduration()[i];
			if(stiduration>0){
				// retrieve proba MTCT
				double proba = _population.STI_proba_MTCT(stiname, stiduration);
				// Draw if transmission occurs
				mtct[i] = (uniform01()<proba? true:false);
			}
		}
	}
	return mtct;
}



void Simulation::update_pregnancies(double timestep)
{
	/// Calculate the number of new pregnancies
	/// during a simulation time step.
	/// Loop through all females who had at least
	/// sex act and draw the number of newly
	/// pregnant women from a binomial distribution
	/// with probability _probaPregnantPerSexAct
	///
	/// Also increase the gestation duration by 'timestep'
	
	// --- NEW PREGNANCIES
	
	// Retrieve all female with potential to get pregnant
	vector<unsigned long> uid_pot_preg = _population.get_UID_pot_preg(); //_population.pregnantPotentialFemales();
	
	// DEBUG
//	unsigned int aa = uid_pot_preg.size();
//	unsigned int bb = _population.get_UID_pot_preg().size();
//	cout<< _simulationTime << " -- DEBUG  count= "<<aa<<" ; bookkeeping= "<<bb<<" ; diff = "<< bb-aa<< endl;
	
	// =====
	
	if (uid_pot_preg.size()>0)
	{
		// loop through all mothers-to-be
		// calculates the proba to become pregnant
		// based on number of sex acts (types 1 or 2).
		// Draw random variable for actual pregnancy
		
		double ppsex = _population.get_probaPregnantPerSexAct();
		
		for (int i=0; i<uid_pot_preg.size(); i++)
		{
			int nsex = _population.get_individual(uid_pot_preg[i]).get_nSexActs_period();    //.nSexActType1or2(); <--- TO DO: WARNING condom not taken into account here!!!
			
			double proba_pregnant = 1 - pow(1 - ppsex, nsex);
			
			double u = uniform01();
			if (u<proba_pregnant)
				set_pregnant(uid_pot_preg[i]);
		}
	}
	
	// --- EXISTING PREGNANCIES
	
	vector<unsigned long> uid_preg = _population.census_pregnant_UID();
	
	_newborn_timestep = 0;  // reset at each time step
	
	for (int i=0; i<uid_preg.size(); i++)
	{
		double gest = _population.get_individual(uid_preg[i]).get_gestationDuration();
		double age_mother = _population.get_individual(uid_preg[i]).get_age();
		
		// Mother-to-child STI transmission.
		// Updated at each time step but
		// effective only at delivery
		vector<bool> sti_mtct = MTCT(uid_preg[i]);
		_population.set_STI_MTCT(uid_preg[i],sti_mtct);

		// increase pregnancy if birth not soon enough
		if (gest < 0.75-timestep)
			increment_gestationDuration(uid_preg[i],timestep);
		
		// trigger child birth if close to gestation period
		if (gest >= 0.75-timestep){
			Individual child;
			_nursery.add_child(child,
							   uid_preg[i],
							   age_mother,
							   _population.get_simulationTime(),
							   sti_mtct);
			
			// Birth book-keeping
			set_gestationDuration(uid_preg[i],0);
			increment_nChildBorn(uid_preg[i]);
			_population.update_STI_mtct_cumcount(sti_mtct);

			vector<bool> reset_mtct(_population.get_nSTImodelled(),false);
			_population.set_STI_MTCT(uid_preg[i],reset_mtct);
			
			_newborn_timestep++;
		}
		// DEBUG
		//cout<< "UID "<< uid_preg[i] << " : " << sti_mtct[0]<<" ; " << sti_mtct[1]<<endl;
		// ------
	}
}


double Simulation::STI_cumul_incidence(STIname s,int time_i)
{
	/// Calculates cumulative incidence up to i^th time point
	
	int sti_pos = positionSTIinVector(s, _population.get_STI());
	stopif(time_i>=_schedule.size(), "Cumulative incidence cannot be calculated beyond simulation horizon!");
	
	double sum = 0.0;
	for (int i=0; i<=time_i; i++)
		sum += _STI_incidence(i,sti_pos);
	
	return sum;
}



void Simulation::runAllEvents_timeStep_obj(int numTimeStep,
										   bool doSex,
										   bool logIndivInfo)
{
	/// Execute all simulation events during a single time step
	
	bool debugflag = false;
	double t = _population.get_simulationTime();
	
	if (debugflag) cout << "SIZE = "<<_population.get_size()<<endl;
	
	// Youth arrivals
	if (debugflag) cout << "youth arrivals"<<endl;
	_population.youthArrivals(_timeStep,_save_trace_files);
	
	// Commercial sex workers (CSW)
	if (debugflag) cout << "CSW recruitment"<<endl;
	_population.CSWrecruitment(_timeStep);
	if (debugflag) cout << "CSW cessation"<<endl;
	_population.CSWcessation(_timeStep);
	
	// Resets variables at the start of this period
	if (doSex) _population.reset_sexActs();
	
	// Partnerships formation and dissolution
	if (debugflag) cout << "dissolve partnerships"<<endl;
	if (_population.totalNumberPartnerships()>0)
		_population.dissolvePartnerships(_timeStep,_save_trace_files);
	
	// Partnerships formations and dissolutions
	if (debugflag) cout << "spousalScanAllCasual"<<endl;
	_population.formPartnerships(_timeStep,_save_trace_files);
	_population.spousalScanAllCasual();
	
	// DEBUG
	_population.debug_check_partnerUID();
	
	// Sex acts
	if (debugflag) cout << "update_sexActs"<<endl;
	if (doSex) _population.update_sexActs(_timeStep, _save_trace_files);
	
	// STI transmissions
	if (debugflag) cout << "STI_transmissions"<<endl;
	vector<unsigned long> incidence;
	if (t>0 && doSex){
		incidence = _population.STI_transmissions(_timeStep,_save_trace_files);
	}
	
	// Pregnancies and MTCT
	if (debugflag) cout << "Pregnancies"<<endl;
	if (t>0 && doSex) update_pregnancies(_timeStep);
	
	// Duration, age updates
	if (debugflag) cout << "updateAllDurations_updateAllAges"<<endl;
	_population.updateAllDurations(_timeStep);
	_population.updateAllAges(_timeStep);
    
    // Genetic mutations
    if (doSex) _population.update_genomes_mutations(_timeStep);
	
	// Natural clearance of all STIs
	if (debugflag) cout << "STI_update_naturalClearance"<<endl;
	_population.STI_update_naturalClearance();

    if (doSex) _population.update_genomes_clear();
    
	// STI Treatment and vaccine waning:
	for (int sti=0; sti<_population.get_STI().size(); sti++){
		STIname stiname = _population.get_STI()[sti].get_name();
		update_cure(stiname);
		update_vacc(stiname);
	}

	// Deaths
	if (debugflag) cout << "deathEvents"<<endl;
	_population.deathEvents(_timeStep,false /*_save_trace_files*/); //force not write trace file because huge!
	
	// update incidence
	if (doSex){
		if (numTimeStep==0){
			unsigned long n_sti = _population.get_STI().size();
			for (int j=0;j<n_sti; j++)
				_STI_incidence(0,j) = 0.0;
		}
		if (numTimeStep>0){
			_STI_incidence.addRowVector(incidence);
		}
	}
	
	// update reproductive numbers for all STIs
	if (doSex) {
		for (int sti=0; sti<_population.get_nSTImodelled(); sti++) {
			STIname stiname = _population.get_STI()[sti].get_name();
			_population.update_secondary_cases(stiname);
		}
	}
	
	vector<double> v;
	
	if(doSex)	// TO DO: do not impose no log file when no sex
	{
		//
		// * * * * WARNING * * * *
		//
		// ORDER MUST BE SAME AS DEFINED
		// BY HEADERS IN "runAllEvents_horizon_obj"
		
		v.push_back(t);
		v.push_back(_population.census_alive());
		v.push_back(_population.census_dead());
		v.push_back(_population.totalNumberPartnerships());
		v.push_back(_population.totalNumberSpousalPartneships());
		v.push_back(_population.census_Females());
		
		int nSTI = _population.get_nSTImodelled();
		for (int i=0; i<nSTI; i++)
			v.push_back(_population.census_STIinfected()[i]);
		
		v.push_back(_population.census_STIcoinfected(HIV, Tp));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,0));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,1));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,2));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,9));

		v.push_back(_population.census_CSW().size());
		
		for (int r=0; r<= _population.get_maxRiskGroup(); r++)
			v.push_back(_population.census_riskGroup()[r]);
		
		v.push_back(_population.census_circum());
		v.push_back(_newborn_timestep);
		
		v.push_back(_population.census_pregnant(0));
		v.push_back(_population.census_pregnant(1));
		v.push_back(_population.census_pregnant(2));
		v.push_back(_population.census_pregnant(9));
		
		v.push_back(_population.get_STI_mtct_cumcount()[0]);
		v.push_back(_population.get_STI_mtct_cumcount()[1]);
		
		// HIV & Tp prevalence:
		v.push_back(_population.STI_prevalence(HIV));
		v.push_back(_population.STI_prevalence(Tp));
		
		// HIV prevalence by risk group
		v.push_back(_population.STI_prevalence(HIV, 0));
		v.push_back(_population.STI_prevalence(HIV, 1));
		v.push_back(_population.STI_prevalence(HIV, 2));
		v.push_back(_population.STI_prevalence(HIV, 9));
		
		// Tp prevalence by risk group
		v.push_back(_population.STI_prevalence(Tp, 0));
		v.push_back(_population.STI_prevalence(Tp, 1));
		v.push_back(_population.STI_prevalence(Tp, 2));
		v.push_back(_population.STI_prevalence(Tp, 9));
		
		// Reproductive numbers
		v.push_back(_population.Reff_cum_mean(HIV));
		v.push_back(_population.Reff_cum_mean(Tp));
		
		// Count of sex acts
		v.push_back(_population.count_sexActs(0));
		v.push_back(_population.count_sexActs(1));
		v.push_back(_population.count_sexActs(2));
		v.push_back(_population.count_sexActs(9));
		
		// * * * * WARNING * * * *
		//
		// ORDER MUST BE SAME AS DEFINED
		// BY HEADERS IN "runAllEvents_horizon_obj"
		// * * * * * * * * * * * *
		
		// Update data frame:
		_df_sim.addrow(to_string(numTimeStep), v);
	}
}


void Simulation::runAllEvents_horizon_obj(bool doSex,
										  bool logIndivInfo,
										  bool traceNetwork,
										  int displayProgress)
{
	/// Run a simulation with all events until specified horizon
	
	unsigned long nSTI = _population.get_nSTImodelled();
	_population.set_timeStep(_timeStep);
	_nursery.set_STI(_population.get_STI());
	
	// Number of interventions
	unsigned long n_intervention = _intervention.size();
	
//	// DEBUG
//	if ( _MC_trial_iter==1){
//		cout << " -- Number of interventions:"<<n_intervention<<endl;
//		for (int i=0; i<n_intervention; i++) {
//			_intervention[i].displayInfo();
//		}
//	}
	
	
	// Set-up data frame that will hold simulation outputs
	vector<string> colnames;
	
	if (doSex)
	{
		// Headers of the data frame
		//
		//  * * * WARNING * * * *
		// headers must be consistent with
		// values writtten in "runAllEvents_timeStep"
		//
		
		colnames.push_back("time");
		colnames.push_back("nAlive");
		colnames.push_back("nDead");
		colnames.push_back("nPartn");
		colnames.push_back("nSp");
		colnames.push_back("nFemale");

		for (int i=0; i<nSTI; i++)
			colnames.push_back(STInameString(_population.get_STI()[i].get_name()));
		
		colnames.push_back("nHIVTp");
		colnames.push_back("nHIVTp0");
		colnames.push_back("nHIVTp1");
		colnames.push_back("nHIVTp2");
		colnames.push_back("nHIVTp9");
		
		colnames.push_back("nCSW");
		
		for (int r=0; r<= _population.get_maxRiskGroup(); r++)
			colnames.push_back("nRskGrp"+int2string(r));
		
		colnames.push_back("nCircum");
		colnames.push_back("nNewBorn");
		
		colnames.push_back("nPregnantRisk0");
		colnames.push_back("nPregnantRisk1");
		colnames.push_back("nPregnantRisk2");
		colnames.push_back("nPregnantRisk9");
		
		colnames.push_back("mtctHIV");
		colnames.push_back("mtctTp");

		colnames.push_back("HIVprev");
		colnames.push_back("Tpprev");
		colnames.push_back("HIVprevRisk0");
		colnames.push_back("HIVprevRisk1");
		colnames.push_back("HIVprevRisk2");
		colnames.push_back("HIVprevRisk9");
		colnames.push_back("TpprevRisk0");
		colnames.push_back("TpprevRisk1");
		colnames.push_back("TpprevRisk2");
		colnames.push_back("TpprevRisk9");
		
		colnames.push_back("Reff_HIV");
		colnames.push_back("Reff_Tp");
		
		colnames.push_back("nSexActRisk0");
		colnames.push_back("nSexActRisk1");
		colnames.push_back("nSexActRisk2");
		colnames.push_back("nSexActRisk9");
		
		_df_sim.set_colname(colnames);
	}

	
	
//	if (logIndivInfo)
//		ff << "time,UID,alive,UIDpartner0, nSexP0, HIVdur,HIVinf"<<endl;
//	// ----
//	
	
	// ===================================================
	// ===    Loop through time until horizon
	// ===================================================
	
//	int cnt = 0;
	int t_i = 0;
	int n_timesteps = (int)(_horizon/_timeStep);
	

	for (double t=0; t<_horizon; t+=_timeStep)
	{
		// write network info
		if (traceNetwork){
//			string file_cnt = filepop+ int2string(cnt) + ".out";
//			_population.FileOutput(file_cnt);
//			cnt++;
		}

		// ---------------------------
		// Display simulation progress
		
		double prevSize = _population.get_size();

		if (displayProgress==1){
			cout << " time: "<< t << " size: "<< prevSize;
			cout << " alive: " << _population.census_alive() <<endl;
		}
		if (displayProgress==11){
			if (t_i%(n_timesteps/10)==0)
			{
				cout << " time: "<< t << " size: "<< prevSize;
				cout << " alive: " << _population.census_alive();// <<endl;
				
				cout << " ; prev_HIV: " << _population.STI_prevalence(HIV);
				cout << " ; prev_Tp: " << _population.STI_prevalence(Tp) <<endl;
				//				cout << " ; Reff_HIV: " << _population.Reff_cum_mean(HIV);
				//				cout << " Reff_Tp: " << _population.Reff_cum_mean(Tp) << endl;
			}
		}
		if(displayProgress==2){
			cout << "alive: " <<	_population.census_alive() << endl;
			cout << "partner0: " << _population.census_Partnered(0) << endl;
			cout << "partner1: " << _population.census_Partnered(1) << endl;
			cout << "partner2: " << _population.census_Partnered(2) << endl;
			cout << "singleRatio="<<_population.census_ratioSingles()<<endl;
		}
		// ---------------------------
		
		// Records the time of simulation
		_simulationTime = t;
		_population.set_simulationTime(t);
		
		
		// ================================================
		// === Execute all events during the time step  ===
		// ================================================
		
		runAllEvents_timeStep_obj(t_i,
								  doSex,
								  logIndivInfo);
		
		// Record STI prevalences
		if (doSex)
		{
			vector<double> tmp;
			for (int s=0; s<_population.get_STI().size(); s++)
			{
				STIname stiname = _population.get_STI()[s].get_name();
				
				if (t==0) _STI_prevalence(t_i,s) = _population.STI_prevalence(stiname);
				if (t>0) tmp.push_back(_population.STI_prevalence(stiname));
			}
			
			if (t>0) _STI_prevalence.addRowVector(tmp);
		}
		
		
		// ========================
		// ==== INTERVENTIONS =====
		// ========================
		
		// Intervention relevant when sex act occurs
		if (doSex && n_intervention>0)
		{
			// loop through all interventions
			for (int i=0; i<n_intervention; i++)
			{
				vector<double> sched = _intervention[i].get_schedule();
				
				// loop through the intervention schedule start/end dates
				for (int k=0; k<sched.size()/2; k++)
				{
					if (sched[2*k]<t && t<sched[2*k+1])
					{
						activate_intervention(i);
						//						//DEBUG
						//						cout << "time: "<<t<<" intervention #"<<i<<" activated (MC iter:"<<_MC_trial_iter<<")"<<endl;
					}
				}
			}
			
		}
		
		// Increase time step counter
		t_i++;
		
	} // end for loop on time
	
	if(_save_trace_files) {}
//		_nursery.saveToCSVfile(_DIR_OUT + "nursery_"+int2string(_MC_trial_iter)+ ".out");
}



/* ***************************************************/
/* ************** C A L I B R A T I O N **************/
/* ***************************************************/


void Simulation::set_calibration_schedule(string filename_calibrationTime,
										  string filename_all_calib_target,
										  string filename_weight_target)
{
	/// Set the calibration schedule
	/// which consists of:
	/// - a vector of times
	/// - for each of those times, a vector of
	/// strings identifying the filename of output
	/// to be calibrated (e.g. age distribution)
	
	/// Example to illustrate:
	///
	/// calibration times:
	/// ____
	/// 1.2
	/// 3.4
	/// 7.8
	///
	/// outputs to be calibrated at those times:
	/// _________________________________________________
	/// ageDist1.csv  |  ageDist2.csv   |  ageDist3.csv
	/// ageGap1.csv   |  ageGap2.csv    |  ageGap3.csv
	///               |  prevTp2.csv    |  prevTp3.csv
	///               |                 |  PrevHIV3.csv
	///
	/// Each file name must have a pre-specified name such that
	/// the _type_ of target is recognized and we can construct
	/// the associated table of calibration types:
	///
	/// _________________________________________________
	/// AGE           |  AGE            |  AGE
	/// AGE           |  GAP            |  GAP
	///               |  TP             |  TP
	///               |                 |  HIV
	///
	/// Finally, each calibration target has its weight
	/// in the global calibration score. For example:
	///
	/// _________________________________________________
	/// 1             |  1              |  1
	/// 1             |  1              |  1
	///               |  1              |  2
	///               |                 |  2
	
	
	// clean up
	_calibrationTime.clear();
	_calibrationTargetFile.clear();
	_calibrationTargetType.clear();
	_calibrationWeights.clear();
	_calibrationDistances.clear();
	
	// read calbration times from the file
	vectorFromCSVfile(_calibrationTime, filename_calibrationTime.c_str(), 1);
	
	// Integrity check
	unsigned int ncaltimes=getNumberColumns(filename_all_calib_target);
	stopif(_calibrationTime.size()!=ncaltimes,
		   "Calibration schedule files not consistent!");
	
	
	// Read all calibration target files
	
	vector< vector<string> > s_filename;
	s_filename.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		vectorFromCSVfile_string(s_filename[j],
								 filename_all_calib_target.c_str(),
								 j+1);
		s_filename[j] = trim(s_filename[j]);
		//displayVector(s_filename[j]);
	}
	_calibrationTargetFile = s_filename;
	
	
	// Convert filenames into calibration types
	
	vector< vector<string> > s_type;
	s_type.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		for(int i=0;i<s_filename[j].size(); i++)
			s_type[j].push_back(calibration_file_to_type(s_filename[j][i]));
		
		//displayVector(s_type[j]);
	}
	_calibrationTargetType = s_type;
	
	
	// Retrieve calibration weights
	
	_calibrationWeights.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		vectorFromCSVfile(_calibrationWeights[j],
						  filename_weight_target.c_str(),
						  j+1);
		
		//displayVector(_calibrationWeights[j]);
	}
	
	// Set correct sizes for calibration distances vectors.
	// Initialized values are set to obviously dummy ones (i.e. negative)
	// to facilitate bug detection
	
	_calibrationDistances.resize(ncaltimes);
	
	for (int i=0; i<ncaltimes; i++) {
		_calibrationDistances[i].resize(_calibrationTargetType[i].size(),-9E9);
	}
	
}

string calibration_file_to_type(string filename)
{
	/// converts a file name into a calibration type
	/// (the file name must obey the rules in this function)
	
	string res="";
	
	if(filename.substr(0,7)	=="ageDist")		res = "AGEDIST";
	if(filename.substr(0,10)=="ageGapDist")		res = "GAPDIST";
	if(filename.substr(0,6)	=="prevTp")			res = "TP";
	
	// USED?:
	if(filename.substr(0,7)	=="prevHIV")		res = "HIV";
	// -----
	
	if(filename.substr(0,12)=="HIVpos_age_f")	res = "HIVDAGEF";
	if(filename.substr(0,12)=="HIVpos_age_m")	res = "HIVDAGEM";
	if(filename.substr(0,14)=="HIV_prev_age_f")	res = "HIVPREVAGEF";
	if(filename.substr(0,14)=="HIV_prev_age_m")	res = "HIVPREVAGEM";
	
	if(filename.substr(0,7)	=="stiPrev")		res = "STIPREV";
	if(filename.substr(0,11)=="stiPrevRisk")	res = "STIPREVRISK";
	
	if(filename.substr(0,11)=="singleRatio")	res = "SINGLERATIO";
	
	if(filename.substr(0,9)=="age1sex_f")		res = "AGE1SEXF";
	if(filename.substr(0,9)=="age1sex_m")		res = "AGE1SEXM";
	
	if(filename.substr(0,7)=="lftNP_f")			res = "LFTNPF";
	if(filename.substr(0,7)=="lftNP_m")			res = "LFTNPM";
	
	if(filename.substr(0,8)=="visitCSW")		res = "VISITCSW";
	if(filename.substr(0,11)=="visitCSWage")	res = "VISITCSWAGE";
	if(filename.substr(0,9)=="condomCSW")		res = "CONDOMCSW";
	
	stopif(res=="","Cannot convert calibration file name ("+ filename+ ") to a calibration type!");
	return res;
}


bool Simulation::isCalibrationTime(double t)
{
	/// Check if this is a calibration time
	
	bool res = false;
	
	for(int i=0;i<_calibrationTime.size();i++)
	{
		if (fabs(t-_calibrationTime[i])<_timeStep/2) res=true;
	}
	return res;
}

unsigned int Simulation::whichCalibrationTime(double t)
{
	/// Given this is a calibration time, which set is it?
	
	unsigned int res = 0;
	
	for(int i=0;i<_calibrationTime.size();i++)
	{
		if (fabs(t-_calibrationTime[i])<_timeStep/2) res=i;
	}
	return res;
}


// ===================================================================
// ===================================================================
// =====  INTERVENTION  =====
// ===================================================================
// ===================================================================

void Simulation::activate_intervention(int i)
{
	/// ACTIVATE TREATMENT OR VACCINATION BASED ON INTERVENTION FEATURES
	/// (assumes simulation time is between start and end dates of intervention)
	/// God-like treatment: assumes we know everyone infected and
	/// everyone has access to intervention
	
	// Integrity checks
	stopif(_intervention[i].get_annCvgRate()<=0, "Proportion of intervention must be positive");
	
	// retrieve intervention features
	double target_dt	= _intervention[i].get_annCvgRate()*_timeStep;
	STIname sti			= _intervention[i].get_stiname();
	string interv_type	= _intervention[i].get_type();
	
	
	// * WARNING *
	// the proportion is understood as a rate and applied every time step
	// For example, propPerYear = 0.5 means 50% of infected will be treated over one year.
	// If the time step is 0.1 year, then the target proportion for that timestep
	// will be 0.1*0.5 = 0.05
	
	
	// check if sampling is actually needed
	// (when prop>1 it means everyone gets treatment during that time step)
	bool do_sample = (target_dt < 1.0);
	
	// count individuals (for log files)
	unsigned long cnt = 0;  // <-- number treated
	unsigned long cnt2 = 0; // <-- number targeted
	
	// scan the intervention type
	bool doTreat_mass		= (interv_type=="treatment_mass");
	bool doTreat_symptom	= (interv_type=="treatment_symptom");
	bool doVacc_mass		= (interv_type=="vaccination_mass");
	bool doVacc_symptom		= (interv_type=="vaccination_symptom");
	bool doVacc_female		= (interv_type=="vaccination_female");
	bool doVacc_femaleYoung	= (interv_type=="vaccination_femaleYoung");
	bool doVacc_hiRisk		= (interv_type=="vaccination_highRisk");
	
	bool doTreatment	= (interv_type.substr(0,9)=="treatment");
	bool doVaccination	= (interv_type.substr(0,11)=="vaccination");
	
	// Integrity checks
	bool type_known =	doTreat_mass || doTreat_symptom ||
						doVacc_mass || doVacc_symptom ||
						doVacc_female || doVacc_femaleYoung ||
						doVacc_hiRisk;
	stopif(!type_known, "Unknown intervention type!");
	
	unsigned int sti_i = positionSTIinVector(sti, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		
		// == Filter who is targeted by intervention ==
		
		Individual indiv = _population.get_individual(uid);
		
		if(indiv.isAlive()){
			
			bool indivIsTargeted	= false;
			
			bool isSymptomatic		= indiv.get_STIsymptom()[sti_i];
			bool alreadyVacc		= indiv.get_STI_vacc()[sti_i];  
			
			if(doTreat_mass)	indivIsTargeted = indiv.get_STIduration()[sti_i]>0; // slow code: indiv.STI_infected(sti);
			if(doTreat_symptom) indivIsTargeted = indiv.get_STIsymptom()[sti_i]; // slow code: is_symptomatic(sti);
			if(doVacc_mass){
				// everyone is targeted
				// but exclude the ones already vaccinated
				indivIsTargeted =  !alreadyVacc;
			}
			if(doVacc_symptom){
				// only symptomatic (for THAT STI) individuals
				indivIsTargeted = isSymptomatic && !alreadyVacc;
			}
			if(doVacc_female){
				// only (all) females
				indivIsTargeted = (indiv.get_gender()==female) && !alreadyVacc;
			}
			if(doVacc_femaleYoung){
				// only _young_ females
				double young_age = 15.0;
				
				bool tmp1 = (indiv.get_gender()==female);
				bool tmp2 = (indiv.get_age()<young_age);
				bool tmp3 = !alreadyVacc;
				indivIsTargeted = tmp1 && tmp2 && tmp3;
			}
			if(doVacc_hiRisk){
				// only the highest risk group
				int mxrg = _population.get_maxRiskGroup();
				indivIsTargeted = (indiv.get_riskGroup() >= mxrg) && !alreadyVacc;
			}
			
			// == Apply intervention on filtered individuals ==
			
			if(indivIsTargeted ){
				cnt2++;
				// no sampling because rate of intervention>1
				// (speed up code execution)
				if (!do_sample){
					if(doTreatment)		treat_indiv(uid, sti);
					if(doVaccination)	vaccinate_indiv(uid, sti);
					cnt++;
				}
				
				// random sample
				if (do_sample && (uniform01()<= target_dt) ){
					if(doTreatment)		treat_indiv(uid, sti);
					if(doVaccination)	vaccinate_indiv(uid, sti);
					cnt++;
				}
			}
		}
	} // end loop on individuals
	
	
	// Intervention information
	vector<double> v;
	vector<string> cnames;

	cnames.push_back("iMC");
	v.push_back(_MC_trial_iter);
	cnames.push_back("time");
	v.push_back(_population.get_simulationTime());
	cnames.push_back("n_reached");
	v.push_back(cnt);
	cnames.push_back("n_targeted");
	v.push_back(cnt2);
	
	string rowname = interv_type + "_" + STInameString(sti);
	_df_interv.addrow(rowname, v);
	if(_df_interv.get_nrows()==1) _df_interv.set_colname(cnames);
}




// ===================================================================
// ===================================================================
// =====  TREATMENT  =====
// ===================================================================
// ===================================================================


void Simulation::treat_indiv(unsigned long uid, STIname sti){
	_population.treat_indiv(uid, sti);
}


void Simulation::cure_indiv(unsigned long uid, STIname stiname)
{
	/// Cure and individual from a given STI
	/// (must be treated beforehand and Adherence>0.8)
	
	stopif(!_population.get_individual(uid).STI_treated(stiname),
		   "Cannot cure without treatment!");
	stopif(!(_population.get_individual(uid).get_STIduration(stiname)>0),
		   "Cannot cure inexistent STI!");
	
	// Cure only if:
	// - treated
	// - treatment is biologically successfull
	// - treatment duration larger than optimal duration
	// - STI is actually curable(!)
	
	
	int i_sti = positionSTIinVector(stiname, _population.get_STI());
	double optim_duration = _population.get_STI()[i_sti].get_optimalTreatmentDuration();
	
	if (_population.get_individual(uid).get_STItreatDuration(stiname) > optim_duration &&
		_population.get_individual(uid).get_STItreatTMS(stiname)==1)
	{
		double A = _population.get_individual(uid).get_STItreatAdherence(stiname);
		if(A>0.8) {
			_population.cure_indiv(uid, stiname);
			//DEBUG
			//cout << "UID "<<uid<< " cured from "<<STInameString(stiname)<<endl;
		}
	}
}


void Simulation::update_cure(STIname sti)
{
	/// Update the cure for all individuals treated against that STI
	/// (mostly when treatment duration has reached optimal duration
	
	int sti_i = positionSTIinVector(sti, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		Individual tmp = _population.get_individual(uid);
		if (tmp.isAlive() &&
			(tmp.get_STItreatDuration()[sti_i]>0) &&
			tmp.get_STIduration()[sti_i]>0) cure_indiv(uid,sti);
		
	}
}



// ===================================================================
// ===================================================================
// =====  VACCINATION  =====
// ===================================================================
// ===================================================================


void Simulation::vaccinate_indiv(unsigned long uid, STIname stiname)
{
	/// Vaccinate an individual against a given STI
	
	_population.vaccinate_indiv(uid, stiname);
}


void Simulation::update_vacc(STIname stiname){
	/// Update immunity provided by vaccination
	
	int sti_i = positionSTIinVector(stiname, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		
		Individual tmp = _population.get_individual(uid);
		
		if (tmp.isAlive() &&
			tmp.get_STI_vacc()[sti_i]){
			
			// immunity resulting from vaccination
			// (1 if success; 0 if failed)
			double imm_vax = tmp.get_STI_immunity(stiname);
			// time since vaccination:
			double tv = _simulationTime - tmp.get_STI_vacc_time()[sti_i];
			// waning rate:
			double w = _population.get_STI()[sti_i].get_vacc_waneRate();
			// immunity:
			double imm = imm_vax * exp(-tv * w);
			// update the value of waning immunity:
			_population.set_STI_immunity(uid, stiname, imm);
			
			// update susceptibility factor:
			double prev_sf = _population.get_STIsusceptFactor(stiname,uid);
			_population.set_STIsusceptFactor(uid, stiname, prev_sf*(1-imm));
		}
	}
}



// ===================================================================
// ===================================================================
// =====  MISCELLENAOUS  =====
// ===================================================================
// ===================================================================

void Simulation::displayInfo()
{
	coutline(80);
	cout<<"   * * * SIMULATION INFORMATION * * * "<<endl;
	coutline(80);
	
	cout << " Horizon:\t"<<_horizon<<endl;
	cout << " TimeStep:\t"<<_timeStep<<endl;
	cout << endl;
	
	cout << " Population information:";
	_population.displayInfo(false);
	
	coutline(80);
	
	cout << " Calibration Schedule"<<endl;
	
	cout << "Calibration Times:";
	displayVector(_calibrationTime);
	
	cout << "Output to be calibrated:"<<endl;
	
	//	for (int i=0; i<_calibrationTime.size(); i++)
	//	{
	//		cout <<endl<< "-- at time " << _calibrationTime[i]<<":";
	//		displayVector(_calibrationTargetFile[i]);
	//	}
	
	coutline(80);
}


