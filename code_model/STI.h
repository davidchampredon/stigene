//
//  STI.h
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 13-09-15.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#ifndef __STIagent_STI__
#define __STIagent_STI__

#include <iostream>
#include <fstream>
#include <string>
#include "RV.h"
#include "dcMatrix.h"
#include "dcTools.h"

#include "gender.h"

//#include "individuals.h"

using namespace std;


// Names of all STI of interest here
// Use abbreviations
enum STIname
{
	// *** W A R N I N G ***
	// changing these names
	// and their order
	// has consequences!!!
	
	HIV,HSV2,HPV,Ct,Ng,Tp,Hd,Bv,Tv,Cv
};

// * * * WARNING * * * 
// UPDATE THE FUNCTIONS BELOW
// IF STIname enum IS CHANGED

string				STInameString(STIname sti);
STIname				StringToSTIname(string stiname);
vector<STIname>		StringToSTIname_vector(vector<string> s);



// Types of STI
enum STItype {virus,bacteria,fungus,protozoa};


class STI
{
	STIname		_name;
	STItype		_type;
	
	string		_genome;
    double      _pathogen_generationTime; // in days
    double      _subst_perNucleotide_perGeneration; // substitutions per nucleotide per generation
	
	// --- Durations ---
	double		_latentDuration;			// latent duration
	double		_latentDuration_2;			// latent duration for 2nd phase of the disease (e.g. secondary syphilis)
	double		_latentDuration_3;			// latent duration for 3rd phase of the disease (e.g. early latent syphilis)
	
	double		_infectiousDuration;		// Duration of infectivity
	double		_infectiousDuration_2;		// Duration of infectivity (secondary stage for syphilis Tp)
	double		_infectiousDuration_2c;		// Duration of infectivity (secondary stage w/ condolymata for syphilis Tp)
	double		_infectiousDuration_3;		// Duration of infectivity (early latent stage for syphilis Tp)
	
	double		_infectiousDuration_f_a;	// Duration of infectivity for asymptomatic females (used for Tv)
	double		_infectiousDuration_f_s;	// Duration of infectivity for symptomatic females (used for Tv)
	double		_infectiousDuration_m_a;	// Duration of infectivity for asymptomatic males (used for Tv)
	double		_infectiousDuration_m_s;	// Duration of infectivity for symptomatic males (used for Tv)
	
	double		_naturalClearanceDuration;	// Mean Duration of natural clearance by immune system

	
	// --- Symptomatic features ---
	bool		_is_symptomatic;			// If this STI is actually symptomatic (when applied to a given Individual)
	double		_proba_symptomatic_female;	// Probability of symptomatic cases in females
	double		_proba_symptomatic_male;	// Probability of symptomatic cases in males
	
	
	// --- Recurrence or persistence ---
	bool		_is_recurrent;				// If this STI is actually recurrent (when applied to a given Individual)
	bool		_is_recurrent_2;			// Second recurrence (i.e. third event -- used for early latent syphilis)
	double		_proba_recurrence;			// Probability of recurrence of symptoms (For syphilis: proba progression to secondary)
	double		_proba_recurrence_2;		// Probability of recurrence of symptoms (For syphilis: proba progression to early latent)
	double		_recurrence_freq;			// Number of recurrences PER YEAR
	double		_recurrence_duration;		// Recurrence duration
	double		_recurrence_gap;			// Duration before recurrence occurs
	
	
	// --- Transmissions ---
	double          _probaMaxSexTransm;		// Maximum probability of transmission per sex act
	double          _minInfectiousness;		// Minimum values for infectivity curve
	double			_maxInfectiousness_Tp1;	// Maximum values for infectivity curve (Syphilis primary stage)
	double			_maxInfectiousness_Tp2;	// Maximum values for infectivity curve (Syphilis secondary stage)
	double			_maxInfectiousness_TpEL;// Maximum values for infectivity curve (Syphilis early latent stage)
	
	vector<double>	_probaSexTransm_SAT;	// Reduction factor from max proba of transmission per sex act
											// to take into account the type of sex act:
											// _probaSexTransm_SAT[0]: reduction compared to high risk sex when condom used (should be tiny number)
											// _probaSexTransm_SAT[1]: reduction compared to high risk sex when no condom used, low risk sex
											// _probaSexTransm_SAT[2]: =1 (no condom used, high risk sex)
	
	vector<double>	_HIVparam;				// Specific parameters for HIV
	
	double			_proba_MTCT;				// probability of mother to child transmission ("vertical")
	
	// --- Infectivity Curve shape parameters ---
	vector<double>	_shape_param;			// Shape parameters for infectivity curve (size depends on STI)
	double			_asymptom_IC_reduc;		// Reduction factor of the Infectivity Curve when infection is asymptomatic
	
	
	// --- Susceptibility ---
	double			_circum_SF_reduction;	// Susceptibility reduction factor when individual circumcised
	
	
	// --- Treatment ---
    bool            _isCurable;                 // Whether this STI is curable or not
    double          _optimalTreatmentDuration;  // Minimum treatment duration for successful cure/heavy suppression
	double          _proba_treatmentFailure;    // Probability of microbiological failure of treatment
	vector<double>  _adherence_param;			// distribution parameters for adherence
    vector<double>  _TREstar_param;             // Treatment Reduction effect parameters (full adherence, microbio non-failure)
	
	// --- Vaccine ---
	double			_proba_vaccineFailure;		// Probability the vaccine forthat STI does not trigger immune response
	double			_VRE;						// Vaccine Reduction Effect on infectivity curve
	double			_vacc_waneRate;				// Immunity exponential waning rate
	
public:
	
	// === CONSTRUCTOR ===
	
	STI (){}
	STI(STIname name, string filename);
	
	STI(STIname name,
		string STIfeatures_filename,
		string treatment_filename);
	
	
	void load_common_param(STIname sti_name, string filename);
	void load_treatment_param(STIname sti_name, string filename);
	void load_vaccine_param(STIname sti_name, string filename);
	
    
    // === GET FUNCTIONS ===
	
	STIname			get_name() {return _name;}
    
    string          get_genome() {return _genome;}

    double          get_pathogen_generationTime() {return _pathogen_generationTime;}
    double          get_subst_perNucleotide_perGeneration() {return _subst_perNucleotide_perGeneration;}
	
	double			get_proba_symptomatic_female() {return _proba_symptomatic_female;}
	double			get_proba_symptomatic_male()   {return _proba_symptomatic_male;}
	double			get_proba_recurrence()         {return _proba_recurrence;}
	
	bool			get_is_symptomatic() {return _is_symptomatic;}
	bool			get_is_recurrent()   {return _is_recurrent;}
	
	double			get_latentDuration()    {return _latentDuration;}
	double			get_infectiousDuration(){return _infectiousDuration;}

	double			get_probaMaxSexTransm()  {return _probaMaxSexTransm;}
	vector<double>	get_probaSexTransm_SAT() {return _probaSexTransm_SAT;}
	double			get_proba_MTCT()         {return _proba_MTCT;}
	
	double			get_naturalClearanceDuration() {return _naturalClearanceDuration;}
	
	// Treatment:
    double			get_proba_treatmentFailure()   {return _proba_treatmentFailure;}
    vector<double>	get_adherence_param()          {return _adherence_param;}
    vector<double>	get_TREstar_param()            {return _TREstar_param;}
	double			get_optimalTreatmentDuration() {return _optimalTreatmentDuration;}
	double			get_circum_SF_reduction()      {return _circum_SF_reduction;}

	// Vaccine:
	double			get_proba_vaccineFailure() {return _proba_vaccineFailure;}
	double			get_VRE() {return _VRE;}
	double			get_vacc_waneRate(){return _vacc_waneRate;}
	
	// === SET FUNCTIONS ===
	
    void            set_genome(string g) {_genome = g;}
    
    void            set_pathogen_generationTime(double x) {_pathogen_generationTime = x;}
    void            set_subst_perNucleotide_perGeneration(double x) {_subst_perNucleotide_perGeneration = x;}
    
	void			set_HIVparam(vector<double> x) {_HIVparam = x;}
	void			set_probaMaxSexTransm(double proba){_probaMaxSexTransm = proba;} // Sets the max proba of transmission
	void			set_probaSexTransm_SAT(vector<double> x) {_probaSexTransm_SAT=x;}
	void			set_minInfectiousness(double x) {_minInfectiousness=x;}
	void			set_naturalClearanceDuration(double x) {_naturalClearanceDuration=x;}
	void			set_proba_symptomatic_female(double proba) {_proba_symptomatic_female = proba;}
	void			set_proba_symptomatic_male(double proba) {_proba_symptomatic_male = proba;}
	void			set_latentDuration(double x) {_latentDuration = x;}
	void			set_infectiousDuration(double x) {_infectiousDuration = x;}
	
	void			set_proba_recurrence(double x) {_proba_recurrence=x;}
	void			set_recurrence_freq(double x) {_recurrence_freq=x;}
	void			set_recurrence_duration(double x) {_recurrence_duration=x;}
	
	void			set_is_symptomatic(bool x) {_is_symptomatic=x;}
	void			set_is_recurrent(bool x) {_is_recurrent = x;}
	void			set_is_recurrent_2(bool x) {_is_recurrent_2=x;}
	
	void			set_proba_treatmentFailure(double x) {_proba_treatmentFailure=x;}
	void			set_optimalTreatmentDuration(double x) {_optimalTreatmentDuration=x;}
	void			set_adherence_param(vector<double> x) {_adherence_param=x;}
	
	void			set_shape_param(vector<double> x) {_shape_param = x;}

	void			set_proba_vaccineFailure(double x) {_proba_vaccineFailure = x;}
	void			set_VRE(double x) {_VRE = x;}
	void			set_vacc_waneRate(double x) {_vacc_waneRate=x;}
	
	// === INFECTIVITY ===
	
	double			infectivityCurve(double stiDuration, Gender g);
	void			check_infectivityCurve(string nameFile);
	
	
    // === TREATMENT ===
	
    double          TREstar(double treatmentDuration);
    
	//  === OTHERS ===
	
	void			displayInfo();

};


// *** OUTSIDE CLASS ***

int positionSTIinVector(STIname s, vector<STI> v);

void check_ALL_infectivityCurves();

double IC(string fct_family, double t, vector<double> param);



#endif /* defined(__STIagent_STI__) */
