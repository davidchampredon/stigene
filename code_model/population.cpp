/*
 *  population.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-07-20.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "population.h"
#include "globalVar.h"
#include "dcDataFrame.h"
#include "genetics.h"

/* ****************************************** */
/* ************** CONSTRUCTORS ************** */
/* ****************************************** */


Population::Population(unsigned long size)
{
	_size = size;
	
	vector<Individual> indiv(_size);
	
	_individual = indiv;
}


void Population::STI_setAllParameters(string filename)
{
	/// SET-UP ALL STIs FOR A THAT POPULATION
	
	// ===== STIs =====
	
	vector<STIname> stiName;
	
	if (getParameterFromFile("do_Bv", filename))   stiName.push_back(Bv);
	if (getParameterFromFile("do_Ct", filename))   stiName.push_back(Ct);
	if (getParameterFromFile("do_Cv", filename))   stiName.push_back(Cv);
	if (getParameterFromFile("do_Hd", filename))   stiName.push_back(Hd);
	if (getParameterFromFile("do_HIV", filename))  stiName.push_back(HIV);
	if (getParameterFromFile("do_HPV", filename))  stiName.push_back(HPV);
	if (getParameterFromFile("do_HSV2", filename)) stiName.push_back(HSV2);
	if (getParameterFromFile("do_Ng", filename))   stiName.push_back(Ng);
	if (getParameterFromFile("do_Tp", filename))   stiName.push_back(Tp);
	if (getParameterFromFile("do_Tv", filename))   stiName.push_back(Tv);
	
	_nSTImodelled = stiName.size();
	
	// Construct the population level template '_STI'
	
	_STI.resize(_nSTImodelled);
	
	for (int s=0; s<_nSTImodelled; s++){
		// Construct each STI (infectious period, proba symptomatic, etc)
		// from 'filename'
		STI tmp(stiName[s], filename);
        cout << "DEBUG: " << filename <<endl;
		_STI[s] = tmp;
	}
	_secondary_cases.resize(_nSTImodelled);
}


void Population::STI_set_genomes(string filename){
    /// Set the genomes of the STIs (for the index case!)
    
    for(int i=0; i<_STI.size(); i++){
        string stiname_i = STInameString(_STI[i].get_name());
        
        // Set genome:
        string genome = getParameterFromFile_string(stiname_i, filename);
        _STI[i].set_genome(genome);
        
        // Set generation time of pathogen:
        double gt = getParameterFromFile(stiname_i + "_generationTime", filename);
        _STI[i].set_pathogen_generationTime(gt);
        
        // Set rate of substitution per nucleotide per generation time:
        double sng = getParameterFromFile(stiname_i + "_sng", filename);
        _STI[i].set_subst_perNucleotide_perGeneration(sng);
    }
}

void Population::set_RebHIV_fromFile(string filename)
{
	_RebHIV.clear();
	
	for(int i=0; i<_nSTImodelled; i++)
		_RebHIV.push_back(getParameterFromFile(STInameString(_STI[i].get_name()), filename));
}


void Population::initFromFile(string pathFile,
							  string STI_filename,
							  string _STI_SFincrease_filename,
							  string RebHIV_filename)
{
	
	/// WARNING ===> MAYBE BE OBSOLETE WITH R WRAPPING
	
	
	/// INITIALIZE POPULATION FROM A CSV FILE
	/// WHICH WAS GENERATED OUTSIDE C++ (SEE 'generateIndividuals.R')
	
	/// FILE MUST CONTAIN:
	/// UID, GENDER, AGE, MAX_NUMBER_PARTNERS,RISKGROUP,
	/// N_LIFETIME_PARTNERS, N_LIFTIME_SPOUSES, IS_DIVORCED, IS_WIDOWED
	
	/// **** ORDER OF COLUMNS IS IMPORTANT ****
	
	/// (BUT NO INFO ON STI NEEDED)
	
	
	// == Read population features as column vectors in file ==
	
	bool removeHeader = true; // The input file is expected to have a header
	
	vector<double> uidtmp;
	vector<double> gendertmp;
	vector<double> agetmp;
	vector<double> maxPartnertmp;
	vector<double> riskGrouptmp;
	vector<double> nLifetimePartnertmp;
	vector<double> isWidowtmp;
	vector<double> isDivorcedtmp;
	vector<double> isCircumtmp;
	
	
	// converting values to right format
	
	vectorFromCSVfile(uidtmp, pathFile.c_str(), 1);
	vectorFromCSVfile(gendertmp, pathFile.c_str(), 2);
	vectorFromCSVfile(agetmp, pathFile.c_str(), 3);
	
	vectorFromCSVfile(maxPartnertmp, pathFile.c_str(), 4);
	vectorFromCSVfile(riskGrouptmp, pathFile.c_str(), 5);
	vectorFromCSVfile(nLifetimePartnertmp, pathFile.c_str(), 6);
	vectorFromCSVfile(isDivorcedtmp, pathFile.c_str(), 8);
	vectorFromCSVfile(isWidowtmp, pathFile.c_str(), 9);
	vectorFromCSVfile(isCircumtmp, pathFile.c_str(), 10);
	
	
	vector<unsigned long> uid		= convertTo_unsignedLong(uidtmp,removeHeader);
	vector<Gender> gender			= convertTo_Gender(gendertmp,removeHeader);
	vector<int> maxPartner			= convertTo_int(maxPartnertmp,removeHeader);
	vector<double> age				= agetmp;
	if (removeHeader) age = deleteElement(agetmp,0); // 0 because remove header which is located as first item
	
	vector<int> riskGroup			= convertTo_int(riskGrouptmp,removeHeader);
	vector<int> nLifetimePartner	= convertTo_int(nLifetimePartnertmp,removeHeader);
	vector<int> isDivorced			= convertTo_int(isDivorcedtmp,removeHeader);
	vector<int> isWidow			= convertTo_int(isWidowtmp,removeHeader);
	vector<int> isCircum			= convertTo_int(isCircumtmp,removeHeader);
	
	_size = uid.size();
	
	
	// == Set all STIs parameters ==
	STI_setAllParameters(STI_filename);
	set_STI_SFincrease(_STI_SFincrease_filename);
	set_RebHIV_fromFile(RebHIV_filename);
	
	
	// ==== Individuals of the population ====
	
	// This construction assumes no STIs
	// (STI status will be determined with epidemic process)
	vector<double> stiDuration(_nSTImodelled,0.0);
	vector<bool> stiSymptom(_nSTImodelled,false);
	
	
	// Initialize all defined values for members of class Individual
	vector<Individual> indiv;
	
	for (int i=0; i<_size; i++){
		// Individual is created
		// no STIs or partners yet
		
		Individual tmp(uid[i],gender[i],age[i],
					   maxPartner[i],riskGroup[i],nLifetimePartner[i],
					   isWidow[i], isDivorced[i], isCircum[i],
					   stiDuration, stiSymptom,
					   _STI,_RebHIV);
		
		stopif(maxPartner[i]==0,"Creating a individual with maxPartner=0 but must be >=1");
		indiv.push_back(tmp);
	}
	_individual = indiv;
	
	
	// === Partneship matrix: no partnerships have been formed at this stage ===
	
	// retrieves the UIDs of all males and females
	vector<unsigned long> females	= getUID(female);
	vector<unsigned long> males		= getUID(male);
	
	_all_UID_femaleUID	= females;
	_all_UID_maleUID	= males;
	
	unsigned long nF = females.size();
	unsigned long nM = males.size();
	
	// Integrity check:
	string errmsg = "Number of females(" + to_string(nF) +")+males("+ to_string(nM) +") != total population("+ to_string(_size) +")!!!";
	stopif (nF+nM != _size,errmsg);
	
	_partnershipsMatrix.resize(0, 0);
	_totalNumberPartnerships = 0;
	_totalNumberSpousalPartnerships = 0;
	
	// === CSW ===
	
	// risk group value assigned to CSW (should be distinctively high)
	_CSWriskGroup = 9;
	
	// Circumcision for the general population
	// set at the average of this starting population
	_proportion_circum = averageElements(isCircumtmp);
}


void Population::initFromFile2(string pathFile,
							   string STI_filename,
							   string _STI_SFincrease_filename,
							   string RebHIV_filename,
							   string STI_treatment_filename,
							   string STI_vaccine_filename){
	
	/// WARNING ===> MAYBE BE OBSOLETE WITH R WRAPPING
	
	/// INITIALIZE POPULATION FROM FILE
	/// AND LOAD STI TREATMENT & VACCINE PARAMETERS
	
	initFromFile( pathFile,  STI_filename,
				 _STI_SFincrease_filename,RebHIV_filename);
	
	for (int s=0; s<_nSTImodelled; s++){
		_STI[s].load_treatment_param(_STI[s].get_name(), STI_treatment_filename);
		_STI[s].load_vaccine_param(_STI[s].get_name(), STI_vaccine_filename);
	}
}


void Population::set_and_check_UID(){
	
	// retrieves the UIDs of all males and females
	_all_UID_femaleUID	= getUID(female);
	_all_UID_maleUID	= getUID(male);
	
	unsigned long nF = _all_UID_femaleUID.size();
	unsigned long nM = _all_UID_maleUID.size();
	
	// Integrity check:
	string errmsg = "Number of females(" + to_string(nF) +")+males("+ to_string(nM) +") != total population("+ to_string(_size) +")!!!";
	stopif (nF+nM != _size,errmsg);
}


void Population::setup_for_simulation(unsigned long founder_size,
									  double founder_female_ratio,
									  double founder_prop_csw,
									  string folder_inputs,
									  string file_STI_features,
									  string file_STI_SFincrease,
									  string file_STI_HIVrebound,
									  string file_STI_treatment,
									  string file_STI_vaccine,
                                      string file_STI_genome,
									  bool debugInfo)
{
	/// SET UP THE INITIAL POPULATION
	/// SO THAT A SIMULATION CAN BE LAUNCHED
	/// WITH THIS POPULATION AS THE STARTING POINT
	
	
	// Set all parameters of this population
	setAllParameters(folder_inputs);

	// Set all STIs parameters
	STI_setAllParameters(folder_inputs + file_STI_features);
    STI_set_genomes     (folder_inputs + file_STI_genome);
    set_STI_SFincrease  (folder_inputs + file_STI_SFincrease);
	set_RebHIV_fromFile (folder_inputs + file_STI_HIVrebound);

	
	// Set treatment and vaccine parameters for each STI:
	for (int s=0; s<_nSTImodelled; s++){
		_STI[s].load_treatment_param(_STI[s].get_name(), folder_inputs+file_STI_treatment);
		_STI[s].load_vaccine_param(_STI[s].get_name(), folder_inputs+file_STI_vaccine);
		_STI_mtct_cumcount.push_back(0);
	}
	
	// Founder ppopulation: no partnerships, no STIs.
	// Reminder:
	// - The simulation will first run to form parnterships
	// according to the partnership rules
	// - Then STI will be introduced and the epidemic will start
	
	create_founder_population(founder_size,
							  founder_female_ratio,
							  founder_prop_csw);
	
	set_and_check_UID();
	
	// update the list of UIDs of females
	// that could potentially become pregnant
	set_UID_pot_preg(pregnantPotentialFemales());
	
	_rec_sexact.clear();
}



void Population::create_founder_population(unsigned long size,
										   double female_ratio,
										   double prop_csw){
	
	/// CREATE THE INITIAL POPULATION
	/// BASED ON POPULATION PARAMETERS (must be already set)
	
	// force same seed such that
	// always start w/ same population for all MC iterations
	force_seed_reset(123456);		// for c++ random functions defined in RV.cpp
	
	_size = size;
	
	// Create UIDs:
	vector<unsigned long> uid;
	for(unsigned long i=0; i<size; i++) uid.push_back(i);

	// Genders:
	unsigned long n_fem = (unsigned long) (size * female_ratio);
	//unsigned long n_mal = size - n_fem;
	vector<Gender> g;
	for(unsigned long i=0; i<n_fem; i++) g.push_back(female);
	for(unsigned long i=n_fem; i<size; i++) g.push_back(female);
	
	// Ages:
	vector<double> age;
	for(unsigned long i=0; i<size; i++)
		age.push_back(_ageSexMin + uniform01()*(_ageSexMax-_ageSexMin));
	
	// -- Risk groups --
	
	vector<unsigned int> rskgrp;
	
	// csw first (because we know firsts are females and csw must be females):
	unsigned long n_csw = (unsigned long) (n_fem * prop_csw);
	_CSWriskGroup = 9;
	for(unsigned long i=0; i<n_csw; i++)
		rskgrp.push_back(_CSWriskGroup);
	
	// Draw risk group values according
	// to the pre-specified risk group proportions:
	vector<unsigned int> rsk_n =  multinomial(size - n_csw, _propRiskGroup); //multinomial_gsl(r, size-n_csw, _propRiskGroup);
	vector<unsigned int> rsk_tmp;
	
	for(unsigned int i=0; i<rsk_n.size(); i++)
		for(unsigned int j=0; j<rsk_n[i]; j++)
			rsk_tmp.push_back(i);
	// Shuffle all risk group values
	// bc they are ordered and thus want to avoid
	// having all same values for females:
	auto engine = std::default_random_engine{};
	shuffle(rsk_tmp.begin(), rsk_tmp.end(),engine);
	
	// Finally, add the (shuffled) risk groups to the
	// first ones created (for CSW):
	for(unsigned int i=0; i<rsk_tmp.size(); i++) rskgrp.push_back(rsk_tmp[i]);
	
	// Maximum number of concurrent sex partners:
	vector<unsigned int> maxCurrSexPartners;
	
	for(unsigned long i=0; i<size; i++){
		double p_maxPrtn = proba_nMaxCurrSexPartner(g[i],rskgrp[i]);
		stopif((p_maxPrtn<=0 || p_maxPrtn>=1),
			   "'_nMaxCurrSexPartner_param' are not properly set because p_maxPrtn is <0 or >1");
		double mxp = (rskgrp[i]==_CSWriskGroup)?999:(1+geometric(p_maxPrtn));
		maxCurrSexPartners.push_back(mxp);
	}
	
	// Proportion circumcised:
	double pc = _proportion_circum;
	vector<double> isCircum;
	// females are not circumcised:
	for(unsigned int i=0; i<n_fem;i++) isCircum.push_back(0);
	// follow the pre-specified circumcision proportion for males:
	for(unsigned int i=n_fem; i<size;i++)
		isCircum.push_back(uniform01()<pc?1:0);
	
	
	// Other variables are here for legacy. TO DO: clean this up!
	vector<unsigned int> nLifetimePartner(size,0);
	vector<unsigned int> nLifetimeSpouse(size,0);
	vector<unsigned int> isDivorced(size,0);
	vector<unsigned int> isWidow(size,0);

	
	// == Create population ==
	
	// No STIs in founding population
	// (will be explictly introduced once
	// the partnership dynamics have reached some kind of equilibrium):
	vector<double>		stiDuration(_nSTImodelled,0.0);
	vector<bool>		stiSymptom(_nSTImodelled,false);
	
	_individual.clear();
	for (int i=0; i<_size; i++){
		// Individual is created
		// no STIs or partners yet
		
		Individual tmp(uid[i],
					   g[i],
					   age[i],
					   maxCurrSexPartners[i],
					   rskgrp[i],
					   nLifetimePartner[i],
					   isWidow[i], isDivorced[i], isCircum[i],
					   stiDuration,
					   stiSymptom,
					   _STI,
                       _RebHIV);
		
		_individual.push_back(tmp);
	}
	
	// No partnerships (because founder population)
	_partnershipsMatrix.resize(0, 0);
	_totalNumberPartnerships = 0;
	_totalNumberSpousalPartnerships = 0;
}




void Population::initSingleDuration(){
	
	for (int i=0; i<_size; i++){
		int nConc = _individual[i].get_nCurrSexPartner();
		
		double duration = 0;
		
		if(nConc==0){
			double age	= _individual[i].get_age();
			double unif = uniform01();
			duration = (age-_ageSexMin)*unif;
		}
		_individual[i].set_singleDuration(duration);
	}
}






/* ****************************************** */
/* ************** GET FUNCTIONS ************* */
/* ****************************************** */


vector<double> Population::get_age()
{
	vector<double> a(_size);
	
	for (int i=0; i<_size; i++)
		a[i] = _individual[i].get_age();
	return a;
}


double Population::getAgeGap(unsigned long uid1, unsigned long uid2)
{
	
	/// CALCULATES THE AGE GAP BETWEEN TWO SEX PARTNERS
	/// BY DEFINITION:
	/// AGE GAP = AGE MALE - AGE FEMALE
	
	if (_individual[uid1].get_gender()==_individual[uid2].get_gender())
	{
		cout <<endl<< "ERROR: cannot calculate Age Gap because ";
		cout << uid1<<" and "<<uid2<< "have same gender!"<<endl;
		exit(1);
	}
	
	unsigned long uid_female = uid1;
	unsigned long uid_male = uid2;
	
	if (_individual[uid1].get_gender()==male)
	{
		uid_female = uid2;
		uid_male = uid1;
	}
	
	double ageGap = _individual[uid_male].get_age()-_individual[uid_female].get_age();
	
	// DEBUG
	//cout << "["<<uid_male<<"--"<<uid_female<< "]Age gap=" <<
	//_individual[uid_male].get_age()<<"-"<<_individual[uid_female].get_age()<<"="<<ageGap<<endl;
	
	return ageGap;
}

vector<unsigned long> Population::getUID(Gender g)
{
	vector<unsigned long> V;
	
	for (unsigned long i=0; i<_size; i++)
	{
		if(_individual[i].get_gender() == g)
		{
			V.push_back(_individual[i].get_UID());
		}
	}
	
	return V;
}

vector<unsigned long> Population::getUID_available(Gender g)
{
	vector<unsigned long> V;
	
	for (unsigned long i=0; i<_size; i++)
	{
		if(_individual[i].isAlive()
		   && _individual[i].get_gender() == g
		   && (_individual[i].get_nCurrSexPartner() < _individual[i].get_nMaxCurrSexPartner() ))
		{
			V.push_back(_individual[i].get_UID());
		}
	}
	
	return V;
}

vector<Individual> Population::subsetIndividual(Gender g)
{
	// Retrieve all Individuals for a given Gender
	
	vector<Individual> V;
	
	for (int i=0; i<_size; i++)
	{
		if(_individual[i].get_gender() == g)
		{
			V.push_back(_individual[i]);
		}
	}
	
	return V;
}




/* ****************************************** */
/* ************** SET FUNCTIONS ************* */
/* ****************************************** */






void Population::setIndividual_UID(vector<unsigned long> uid)
{
	if (_size>0 && _size!=uid.size())
	{
		cout << "ERROR: setIndividual_UID: setting UID with unmatched size"<<endl;
		exit(1);
	}
	
	if (_size==0)
	{
		_size = uid.size();
		_individual.resize(_size);
		
		for (int i=0; i<_size; i++)
			_individual[i].set_UID(uid[i]);
	}
}

void Population::setIndividual_gender(vector<Gender> gender)
{
	if (_size>0 && _size!=gender.size())
	{
		cout << "ERROR: setIndividual_gender: setting Genders with unmatched size"<<endl;
		exit(1);
	}
	
	if (_size==0)
	{
		_size = gender.size();
		_individual.resize(_size);
		
		for (int i=0; i<_size; i++)
			_individual[i].set_gender(gender[i]);
	}
	
}

void Population::setIndividual_age(vector<double> age)
{
	if (_size>0 && _size!=age.size())
	{
		cout << "ERROR: setIndividual_age: setting Ages with unmatched size (ages size="
		<< age.size()<<" vs. Population size="<<_size<<")"<<endl;
		exit(1);
	}
	
	if (_size==0)
	{
		_size = age.size();
		_individual.resize(_size);
		
		for (int i=0; i<_size; i++)
			_individual[i].set_age(age[i]);
	}
}

void Population::setIndividual_nMaxCurrSexPartner(vector<int> nMaxCurrSexPartner)
{
	if (_size>0 && _size!=nMaxCurrSexPartner.size())
	{
		cout << "ERROR: setIndividual_nMaxCurrSexPartner: setting Max n partner with unmatched size"<<endl;
		exit(1);
	}
	
	if (_size==0)
	{
		_size = nMaxCurrSexPartner.size();
		_individual.resize(_size);
		
		for (int i=0; i<_size; i++)
		{
			_individual[i].set_nMaxCurrSexPartner(nMaxCurrSexPartner[i]);
			// DEBUG
			if (nMaxCurrSexPartner[i]==0)
			{
				;//DEBUG
			}
		}
	}
}


void Population::setIndividual_nCurrSexPartner()
{
	for (int i=0; i<_size; i++)
		_individual[i].set_nCurrSexPartner(0);
}



/******************* PARAMETERS *******************/



void Population::setAllParameters(string folder_inputs)
{
	/// SET ALL PARAMETERS AT ONCE
	/// FILE NAMES ARE HARD CODED (CHANGE THAT?)
	
	// File names
	
	string populationFeatures_file		= folder_inputs + "in_populationFeatures.csv";
	
	string demographics_file			= folder_inputs + "in_paramDMG.csv";
	string formation_file				= folder_inputs + "in_paramFORM.csv";
	string spousal_file					= folder_inputs + "in_paramSPOUSAL.csv";
	string dissolution_file				= folder_inputs + "in_paramDISSOL.csv";
	
	string param_sexActivity_file		= folder_inputs + "in_paramSexActivity.csv";
	string param_CSW_file				= folder_inputs + "in_CSW.csv";
	
	
	// Set values of parameters
	
	set_Population_features(populationFeatures_file);
	
	set_Demographic_parameters(demographics_file);
	set_Formation_parameters(formation_file);
	set_SpousalProgress_parameters(spousal_file);
	set_Dissolution_parameters(dissolution_file);
	
	set_SexualActivity_parameters(param_sexActivity_file);
	set_CSW_demographics_parameters(param_CSW_file);
}



void Population::set_Population_features(string filename)
{
	_maxAge = getParameterFromFile("maxAge", filename);
	
	// Define maximum risk group value
	_maxRiskGroup = getParameterFromFile("maxRiskGroup", filename);
	
	// Define age limits of sexual activity
	_ageSexMin = getParameterFromFile("minSexAge", filename);
	_ageSexMax = getParameterFromFile("maxSexAge", filename);
}


void Population::set_Demographic_parameters(string filename)
{
	_birthRate			= getParameterFromFile("birthRate",filename);
	_infantMortality	= getParameterFromFile("infantMortalityRate", filename);
	_childMortality		= getParameterFromFile("childMortalityRate", filename);
	
	_deathParam_Weibull_shape		= getParameterFromFile("death_Weibull_shape", filename);
	_deathParam_Weibull_scale		= getParameterFromFile("death_Weibull_scale", filename);
	_deathParam_Weibull_shape_hiv	= getParameterFromFile("death_Weibull_shape_hiv", filename);
	_deathParam_Weibull_scale_hiv	= getParameterFromFile("death_Weibull_scale_hiv", filename);
	
	_maxPregnantAge					= getParameterFromFile("_maxPregnantAge",filename);
	_probaPregnantPerSexAct			= getParameterFromFile("_probaPregnantPerSexAct",filename);
	
}


void Population::set_Formation_parameters(string filename)
{
	_formation_MaxRate = getParameterFromFile("formation_MaxRate", filename);
	
	_formation_age_fullstart	= getParameterFromFile("formation_age_fullstart", filename);		// Age when female is fully age-attractive
	_formation_age_pivot		= getParameterFromFile("formation_age_pivot", filename);			// Parameter representing age where female is half age-attractive (logistic fct)
	_formation_age_shape		= getParameterFromFile("formation_age_shape", filename);			// Shape parameter for female age attractiveness
	_formation_age_fmin			= getParameterFromFile("formation_age_fmin", filename);
	
	_formation_agegap_mean	= getParameterFromFile("formation_agegap_mean", filename);
	_formation_agegap_fmin	= getParameterFromFile("formation_agegap_fmin", filename);
	
	_formation_agegap_shape.push_back(getParameterFromFile("_formation_agegap_shape_1", filename)); //min age gap
	_formation_agegap_shape.push_back(getParameterFromFile("_formation_agegap_shape_2", filename)); // "a"
	_formation_agegap_shape.push_back(getParameterFromFile("_formation_agegap_shape_3", filename)); // "d"
	
	vector<double> tmp(0);
	tmp.push_back(getParameterFromFile("formation_RiskGroup_1", filename));
	tmp.push_back(getParameterFromFile("formation_RiskGroup_2", filename));
	
	_formation_RiskGroup = tmp;
	
	tmp.clear();
	tmp.push_back(getParameterFromFile("formation_PartnDeficit_1", filename));
	_formation_PartnDeficit = tmp;
	
	_formation_STIsymptom_f = getParameterFromFile("formation_STIsymptom_f", filename);
	_formation_STIsymptom_m = getParameterFromFile("formation_STIsymptom_m", filename);
}


void Population::set_SpousalProgress_parameters(string filename)
{
	_spousalProgress_maxRate		= getParameterFromFile("spousalProgress_maxRate", filename);
	_spousalProgress_meanAge_f		= getParameterFromFile("spousalProgress_meanAge_f", filename);
	_spousalProgress_varAge_f		= getParameterFromFile("spousalProgress_varAge_f", filename);
	_spousalProgress_meanGap		= getParameterFromFile("spousalProgress_meanGap", filename);
	_spousalProgress_varGap			= getParameterFromFile("spousalProgress_varGap", filename);
	_spousalProgress_durationK1		= getParameterFromFile("spousalProgress_durationK1", filename);
	_spousalProgress_durationK2		= getParameterFromFile("spousalProgress_durationK2", filename);
	_spousalProgress_meanDiffAgeGap	= getParameterFromFile("spousalProgress_meanDiffAgeGap", filename);
	_spousalProgress_varDiffAgeGap	= getParameterFromFile("spousalProgress_varDiffAgeGap", filename);
}



void Population::set_Dissolution_parameters(string filename)
{
	_dissolution_MaxRate		= getParameterFromFile("dissolution_MaxRate", filename);
	_dissolution_spouse			= getParameterFromFile("dissolution_spouse", filename);
	_dissolution_STI_symptom	= getParameterFromFile("dissolution_STI_symptom", filename);
	
	vector<double> tmp(0);
	tmp.push_back(getParameterFromFile("dissolution_RiskGroup_1", filename));
	_dissolution_RiskGroup	= tmp;
	
	tmp.clear();
	tmp.push_back(getParameterFromFile("dissolution_duration_1", filename));
	tmp.push_back(getParameterFromFile("dissolution_duration_2", filename));
	tmp.push_back(getParameterFromFile("dissolution_duration_3", filename));
	_dissolution_duration	= tmp;
	
	tmp.clear();
	tmp.push_back(getParameterFromFile("dissolution_age_1", filename));
	tmp.push_back(getParameterFromFile("dissolution_age_2", filename));
	tmp.push_back(getParameterFromFile("dissolution_age_probMin", filename));
	_dissolution_age = tmp;
	
	tmp.clear();
	tmp.push_back(getParameterFromFile("dissolution_PartnerDeficit_1", filename));
	tmp.push_back(getParameterFromFile("dissolution_PartnerDeficit_2", filename));
	_dissolution_PartnerDeficit = tmp;
}



void Population::set_SexualActivity_parameters(string param_sexActivity_file)
{
	/// SET ALL PARAMETERS RELATED TO SEXUAL BEHAVIOUR
	/// READ FROM A FILE
	
	
	/// Proportion of risk groups
	/// risk groups: 0,1,2,...,_maxRiskGroup
	/// (does not include CSW because managed seperately)
	
	_propRiskGroup.clear();
	for (int i=0; i<=_maxRiskGroup; i++){
		string s = "_propRiskGroup_" + int2string(i);
		_propRiskGroup.push_back(getParameterFromFile(s,param_sexActivity_file));
	}
	
	// Maximum rate of sexual intercourse
	set_sexActMaxRate(getParameterFromFile("_sexAct_maxRate_male",param_sexActivity_file),
					  getParameterFromFile("_sexAct_maxRate_female",param_sexActivity_file));
	
	// Maximum number of concurrent partners
	_nMaxCurrSexPartner_param_f.clear();
	_nMaxCurrSexPartner_param_f.push_back(getParameterFromFile("_nMaxCurrSexPartner_Female_1", param_sexActivity_file));
	_nMaxCurrSexPartner_param_f.push_back(getParameterFromFile("_nMaxCurrSexPartner_Female_2", param_sexActivity_file));
	
	_nMaxCurrSexPartner_param_m.clear();
	_nMaxCurrSexPartner_param_m.push_back(getParameterFromFile("_nMaxCurrSexPartner_Male_1", param_sexActivity_file));
	_nMaxCurrSexPartner_param_m.push_back(getParameterFromFile("_nMaxCurrSexPartner_Male_2", param_sexActivity_file));
	
	// Preference for sex act with spouses
	set_sexAct_proba_distribute_partnerTypes_prefSpouse(getParameterFromFile("_sexAct_proba_distribute_partnerTypes_prefSpouse",param_sexActivity_file));
	
	// reduction of sexual activity associated with age
	vector<double> sexred_age(3);
	sexred_age[0] = getParameterFromFile("_sexAct_reduce_age_param_1",param_sexActivity_file);
	sexred_age[1] = getParameterFromFile("_sexAct_reduce_age_param_2",param_sexActivity_file);
	sexred_age[2] = getParameterFromFile("_sexAct_reduce_age_param_3",param_sexActivity_file);
	set_sexAct_reduce_age_param(sexred_age);
	
	// reduction of sexual activity associated with risk group
	vector<double> sexred_risk(1);
	sexred_risk[0] = getParameterFromFile("_sexAct_reduce_risk_param_1",param_sexActivity_file);
	set_sexAct_reduce_risk_param(sexred_risk);
	
	// reduction of sexual activity associated with STI symptoms
	vector<double> sexred_sympt(2);
	sexred_sympt[0] = getParameterFromFile("_sexAct_reduce_STIsymptom_param_male",param_sexActivity_file);
	sexred_sympt[1] = getParameterFromFile("_sexAct_reduce_STIsymptom_param_female",param_sexActivity_file);
	set_sexAct_reduce_STIsymptom_param(sexred_sympt);
	
	// reduction of sexual activity associated with number of partners
	vector<double> sexred_np(1);
	sexred_np[0] = getParameterFromFile("_sexAct_reduce_nPartner_param",param_sexActivity_file);
	set_sexAct_reduce_nPartner_param(sexred_np);
	
	// Preference parameter for sex acts with sex workers
	vector<double> x(3);
	x[0] = getParameterFromFile("_sexAct_proba_distribute_partnerTypes_sexWorkParam_1",param_sexActivity_file);
	x[1] = getParameterFromFile("_sexAct_proba_distribute_partnerTypes_sexWorkParam_2",param_sexActivity_file);
	set_sexAct_proba_distribute_partnerTypes_sexWorkParam(x);
	
	// parameters defining the probability of intercourse with sex workers
	vector<double> xx(2);
	xx[0] = getParameterFromFile("_sexAct_proba_sexWorker_param_1",param_sexActivity_file);
	xx[1] = getParameterFromFile("_sexAct_proba_sexWorker_param_2",param_sexActivity_file);
	set_sexAct_proba_sexWorker_param(xx);
	
	// reduction factor of number of sexual acts due to commercial sex costs
	set_sexAct_CostSexWork_reduction(getParameterFromFile("_sexAct_CostSexWork_reduction",param_sexActivity_file));
	
	// Probability for condom use
	vector<double> condom;
	condom.push_back(getParameterFromFile("_sexAct_condom_param_1",param_sexActivity_file));
	condom.push_back(getParameterFromFile("_sexAct_condom_param_2",param_sexActivity_file));
	set_sexAct_condom_param(condom);
	
	// Ratio low risk vs high risk sex
	set_sexAct_TypeLowHigh(getParameterFromFile("_sexAct_TypeLowHigh",param_sexActivity_file));
	
}


void Population::set_CSW_demographics_parameters(string param_CSW_file)
{
	set_CSW_maxRecruitmentRate(getParameterFromFile("_CSW_maxRecruitmentRate",param_CSW_file));
	set_CSW_recruitment_saturPrm_1(getParameterFromFile("_CSW_recruitment_saturPrm_1",param_CSW_file));
	set_CSW_recruitment_saturPrm_2(getParameterFromFile("_CSW_recruitment_saturPrm_2",param_CSW_file));
	set_CSW_recruitment_ageMin(getParameterFromFile("_CSW_recruitment_ageMin",param_CSW_file));
	set_CSW_recruitment_ageMax(getParameterFromFile("_CSW_recruitment_ageMax",param_CSW_file));
	set_CSW_cessationRate(getParameterFromFile("_CSW_cessationRate",param_CSW_file)) ;
	set_CSW_ageMax(getParameterFromFile("_CSW_ageMax",param_CSW_file));
}





void Population::set_spousalProgress_Age_param(double m_af, double v_af, double m_g, double v_g)
{
	_spousalProgress_meanAge_f = m_af;
	_spousalProgress_varAge_f = v_af;
	_spousalProgress_meanGap = m_g;
	_spousalProgress_varGap = v_g;
}

void Population::set_spousalProgress_Duration_param(double k1, double k2)
{
	_spousalProgress_durationK1 = k1;
	_spousalProgress_durationK2 = k2;
}

void Population::set_spousalProgress_DiffAgeGap_param(double meanDelta, double varDelta)
{
	_spousalProgress_meanDiffAgeGap = meanDelta;
	_spousalProgress_varDiffAgeGap = varDelta;
}

void Population::setSpousalParameters(string spousalFile)
{
	// WARNING: Number of parameters read in file and order defined here
	
	vector<double> prm(8);
	vectorFromCSVfile(prm, spousalFile.c_str(), 2);
	
	set_spousalProgress_Age_param(prm[0],prm[1],prm[2],prm[3]);
	set_spousalProgress_Duration_param(prm[4],prm[5]);
	set_spousalProgress_DiffAgeGap_param(prm[6],prm[7]);
}



void Population::increment_nChildBorn(unsigned long uid)
{
	/// Event associated with child birth for mother
	
	stopif(_individual[uid].get_gender()==male, "A child cannot be born from a male!");
	
	_individual[uid].set_isPregnant(false);
	_individual[uid].set_nChildBorn(_individual[uid].get_nChildBorn()+1);
	
	// Update list of females potential pregnant
	if(_individual[uid].get_age()<_maxPregnantAge) _UID_pot_preg.push_back(uid);
}

void Population::set_STI_SFincrease(string filename)
{
	dcMatrix tmp(_nSTImodelled);
	tmp.FromFile(filename);
	
	_STI_SFincrease = tmp;
}


void Population::set_HIVparam(vector<double> param)
{
	// SETS SPECIFIC HIV PARAMETER INTO
	// THE POPULATION-LEVEL STI TEMPLATE
	
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==HIV)
		{
			_STI[i].set_HIVparam(param);
		}
	}
}

void Population::set_STI_probaMaxSexTransm(STIname name, double proba)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_probaMaxSexTransm(proba);
		}
	}
}

void Population::set_STI_probaSexTransm_SAT(STIname name, vector<double> param)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_probaSexTransm_SAT(param);
		}
	}
}

void Population::set_STI_naturalClearanceDuration(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_naturalClearanceDuration(x);
		}
	}
}

void Population::set_STI_proba_symptomatic_female(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_proba_symptomatic_female(x);
		}
	}
}

void Population::set_STI_proba_symptomatic_male(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_proba_symptomatic_male(x);
		}
	}
}

void Population::set_STI_latentDuration(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_latentDuration(x);
		}
	}
}

void Population::set_STI_infectiousDuration(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_infectiousDuration(x);
		}
	}
}




void Population::set_STI_proba_recurrence(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_proba_recurrence(x);
		}
	}
}
void Population::set_STI_recurrence_freq(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_recurrence_freq(x);
		}
	}
}
void Population::set_STI_recurrence_duration(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_recurrence_duration(x);
		}
	}
}


void Population::set_STI_minInfectiousness(STIname name, double x)
{
	for (int i=0; i<_nSTImodelled; i++)
	{
		if (_STI[i].get_name()==name)
		{
			_STI[i].set_minInfectiousness(x);
		}
	}
}



/* ****************************************** */
/* *********  PARAMETERS UPDATE  ************ */
/* ****************************************** */



void Population::UpdateSelectedParameter(string paramName,
										 double value)
{
	bool found = false;
	
	// =========================
	// DEMOGRAPHIC PARAMETERS
	// =========================
	
	if (paramName == "_birthRate") {set_birthRate(value); found=true;}
	if (paramName == "_deathParam_Weibull_scale") {set_deathParam_Weibull_scale(value); found=true;}
	if (paramName == "_deathParam_Weibull_scale_hiv") {set_deathParam_Weibull_scale_hiv(value); found=true;}
	if (paramName == "_deathParam_Weibull_shape") {set_deathParam_Weibull_shape(value); found=true;}
	if (paramName == "_deathParam_Weibull_shape_hiv") {set_deathParam_Weibull_shape_hiv(value); found=true;}
	
	
	// =========================
	// PARTNERSHIP PARAMETERS
	// =========================
	
	
	// --- FORMATION ---
	
	if (paramName == "_formation_MaxRate") {set_formation_MaxRate(value); found=true;}
	
	if (paramName == "_formation_age_fullstart") {set_formation_age_fullstart(value); found=true;}
	if (paramName == "_formation_age_pivot") {set_formation_age_pivot(value); found=true;}
	if (paramName == "_formation_age_shape") {set_formation_age_shape(value); found=true;}
	if (paramName == "_formation_age_fmin") {set_formation_age_fmin(value); found=true;}
	if (paramName == "_formation_agegap_mean") {set_formation_agegap_mean(value); found=true;}
	if (paramName == "_formation_agegap_var") {set_formation_agegap_var(value); found=true;}
	if (paramName == "_formation_agegap_fmin") {set_formation_agegap_fmin(value); found=true;}
	if (paramName == "_formation_agegap_shape_1") {set_formation_agegap_shape(value,0); found=true;}
	if (paramName == "_formation_agegap_shape_2") {set_formation_agegap_shape(value,1); found=true;}
	if (paramName == "_formation_agegap_shape_3") {set_formation_agegap_shape(value,2); found=true;}
	
	if (paramName == "_formation_RiskGroup_1") {set_formation_RiskGroup(value,0); found=true;}
	if (paramName == "_formation_RiskGroup_2") {set_formation_RiskGroup(value,1); found=true;}
	if (paramName == "_formation_PartnDeficit_1") {set_formation_PartnDeficit(value,0); found=true;}
	if (paramName == "_formation_PartnDeficit_2") {set_formation_PartnDeficit(value,1); found=true;}
	if (paramName == "_formation_STIsymptom_f") {set_formation_STIsymptom_f(value); found=true;}
	if (paramName == "_formation_STIsymptom_m") {set_formation_STIsymptom_m(value); found=true;}
	
	
	// --- SPOUSAL ---
	
	if (paramName == "_spousalProgress_maxRate") {set_spousalProgress_maxRate(value); found=true;}
	if (paramName == "_spousalProgress_meanAge_f") {set_spousalProgress_meanAge_f(value); found=true;}
	if (paramName == "_spousalProgress_varAge_f") {set_spousalProgress_varAge_f(value); found=true;}
	if (paramName == "_spousalProgress_meanGap") {set_spousalProgress_meanGap(value); found=true;}
	if (paramName == "_spousalProgress_varGap") {set_spousalProgress_varGap(value); found=true;}
	if (paramName == "_spousalProgress_durationK1") {set_spousalProgress_durationK1(value); found=true;}
	if (paramName == "_spousalProgress_durationK2") {set_spousalProgress_durationK2(value); found=true;}
	if (paramName == "_spousalProgress_meanDiffAgeGap") {set_spousalProgress_meanDiffAgeGap(value); found=true;}
	if (paramName == "_spousalProgress_varDiffAgeGap") {set_spousalProgress_varDiffAgeGap(value); found=true;}
	
	
	
	// --- DISSOLUTION ---
	
	if (paramName == "_dissolution_MaxRate") {set_dissolution_MaxRate(value); found=true;}
	if (paramName == "_dissolution_spouse") {set_dissolution_spouse(value); found=true;}
	if (paramName == "_dissolution_STI_symptom") {set_dissolution_STI_symptom(value); found=true;}
	if (paramName == "_dissolution_ageConcPartn") {set_dissolution_ageConcPartn(value); found=true;}
	
	if (paramName == "_dissolution_RiskGroup_1") {set_dissolution_RiskGroup(value,1); found=true;}
	
	if (paramName == "_dissolution_duration_1") {set_dissolution_duration(value,0); found=true;}
	if (paramName == "_dissolution_duration_2") {set_dissolution_duration(value,1); found=true;}
	if (paramName == "_dissolution_duration_3") {set_dissolution_duration(value,2); found=true;}
	
	if (paramName == "_dissolution_age_1") {set_dissolution_age(value,0); found=true;}
	if (paramName == "_dissolution_age_2") {set_dissolution_age(value,1); found=true;}
	if (paramName == "_dissolution_age_3") {set_dissolution_age(value,2); found=true;}
	
	if (paramName == "_dissolution_PartnerDeficit_1") {set_dissolution_PartnerDeficit(value,0); found=true;}
	if (paramName == "_dissolution_PartnerDeficit_2") {set_dissolution_PartnerDeficit(value,1); found=true;}
	
	
	
	// =========================
	// CSW PARAMETERS
	// =========================
	
	
	if (paramName == "_CSW_maxRecruitmentRate") {set_CSW_maxRecruitmentRate(value); found=true;}
	if (paramName == "_CSW_cessationRate") {set_CSW_cessationRate(value); found=true;}
	
	
	
	// =========================
	// SEXUAL BEHAVIOUR
	// =========================
	
	if (paramName == "_propRiskGroup_0") {set_propRiskGroup(value,0); found=true;}
	if (paramName == "_propRiskGroup_1") {set_propRiskGroup(value,1); found=true;}
	if (paramName == "_propRiskGroup_2") {set_propRiskGroup(value,2); found=true;}
	
	if (paramName == "_nMaxCurrSexPartner_param_f_1") {set_nMaxCurrSexPartner_param_f(value,0); found=true;}
	if (paramName == "_nMaxCurrSexPartner_param_f_2") {set_nMaxCurrSexPartner_param_f(value,1); found=true;}
	if (paramName == "_nMaxCurrSexPartner_param_m_1") {set_nMaxCurrSexPartner_param_m(value,0); found=true;}
	if (paramName == "_nMaxCurrSexPartner_param_m_2") {set_nMaxCurrSexPartner_param_m(value,1); found=true;}
	
	if (paramName == "_sexAct_maxRate_female") {set_sexAct_maxRate_female(value); found=true;}
	if (paramName == "_sexAct_maxRate_male") {set_sexAct_maxRate_male(value); found=true;}
	
	if (paramName == "_sexAct_reduce_age_param_1") {set_sexAct_reduce_age_param(value,0); found=true;}
	if (paramName == "_sexAct_reduce_age_param_2") {set_sexAct_reduce_age_param(value,1); found=true;}
	if (paramName == "_sexAct_reduce_age_param_3") {set_sexAct_reduce_age_param(value,2); found=true;}
	
	if (paramName == "_sexAct_reduce_risk_param_1") {set_sexAct_reduce_risk_param(value,0); found=true;}
	if (paramName == "_sexAct_reduce_risk_param_2") {set_sexAct_reduce_risk_param(value,1); found=true;}
	
	if (paramName == "_sexAct_reduce_STIsymptom_param_1") {set_sexAct_reduce_STIsymptom_param(value,0); found=true;}
	if (paramName == "_sexAct_reduce_STIsymptom_param_2") {set_sexAct_reduce_STIsymptom_param(value,1); found=true;}
	
	if (paramName == "_sexAct_reduce_nPartner_param_1") {set_sexAct_reduce_nPartner_param(value,0); found=true;}
	if (paramName == "_sexAct_reduce_nPartner_param_2") {set_sexAct_reduce_nPartner_param(value,1); found=true;}
	
	if (paramName == "_sexAct_proba_distribute_partnerTypes_prefSpouse")
	{set_sexAct_proba_distribute_partnerTypes_prefSpouse(value); found=true;}
	
	if (paramName == "_sexAct_proba_distribute_partnerTypes_sexWorkParam_1")
	{set_sexAct_proba_distribute_partnerTypes_sexWorkParam(value,0); found=true;}
	if (paramName == "_sexAct_proba_distribute_partnerTypes_sexWorkParam_2")
	{set_sexAct_proba_distribute_partnerTypes_sexWorkParam(value,1); found=true;}
	
	if (paramName == "_sexAct_proba_sexWorker_param_1") {set_sexAct_proba_sexWorker_param(value,0); found=true;}
	if (paramName == "_sexAct_proba_sexWorker_param_2") {set_sexAct_proba_sexWorker_param(value,1); found=true;}
	
	if (paramName == "_sexAct_CostSexWork_reduction") {set_sexAct_CostSexWork_reduction(value); found=true;}
	
	
	if (paramName == "_sexAct_condom_param_1") {set_sexAct_condom_param(value,0); found=true;}
	if (paramName == "_sexAct_condom_param_2") {set_sexAct_condom_param(value,0); found=true;}
	
	if (paramName == "_sexAct_TypeLowHigh") {set_sexAct_TypeLowHigh(value); found=true;}
	
	// If parameter was not found -> error message
	string errmsg = "Parameter '" + paramName +"' not found";
	stopif(!found,errmsg);
}


void Population::UpdateSelectedParameter(vector<string> paramName,
										 vector<double> value)
{
	stopif(paramName.size() != value.size(), "'paramName' and 'value' vector size must be equal");
	
	for (int i=0; i<paramName.size(); i++)
		UpdateSelectedParameter(paramName[i],value[i]);
}


void Population::UpdateSelectedParameter_file(string filename, unsigned int paramSet)
{
	/// UPDATE THE VALUE OF SELECTED PARAMETERS
	/// DEFINED IN A FILE
	///
	/// _COMPULSORY_ FILE FORMAT:
	///
	/// FIRST COLUMN: PARAMETER NAME
	/// SECOND COLUMN: VALUES FOR EACH PARAMETERS (SET #1)
	/// THIRD COLUMN: VALUES FOR EACH PARAMETERS (SET #2)
	/// ...
	/// Nth COLUMN: VALUES FOR EACH PARAMETERS (SET #N-1)
	
	bool headers = false;
	dcDataFrame df(filename,headers);
	
	vector<string>  paramName = df.get_rowname(); //df.get_colname();
	vector<double>  values = df.get_value().extractColumn(paramSet-1);
	
	UpdateSelectedParameter(paramName,values);
}





/* ****************************************** */
/* ************** DEMOGRAPHICS ************** */
/* ****************************************** */


void Population::addIndividual(Individual I)
{
	/// Add an individual to the population modelled
	/// Individual "I" must have all its attributes (age, gender,etc)
	/// already defined.
	/// BUT, the UID will be assigned here and any other value
	/// initialized before will be discarded
	
	// The new UID is set to maxUID+1:
	unsigned long newUID = getMaxUID()+1;
	I.set_UID(newUID);
	I.set_dateInPopulation(_simulationTime);
	
	// increase the size
	_size ++;
	
	// Add the individual to the list (vector) of individuals recorded
	_individual.push_back(I);
	
	// update the list of potential pregnant females:
	if(I.get_gender()==female && I.get_age()<_maxPregnantAge && !I.get_isPregnant()) _UID_pot_preg.push_back(newUID);
	
	// Update what's necessary  ....
	
}

unsigned long Population::getMaxUID(){
	/// Returns the largest UID (alive or not)
	return _size-1;
}



double	Population::probaSex_sexworker(int riskgroup)
{
	double c1 = _sexAct_proba_sexWorker_param[0];
	double c2 = _sexAct_proba_sexWorker_param[1];
	
	return c1*exp(-c2*(_maxRiskGroup-riskgroup));
}


double Population::proba_nMaxCurrSexPartner(Gender g, int riskgroup)
{
	/// Calculate probability used in geometric distribution
	/// of max number allowed of concurrent partnerships
	
	double p_maxPrtn = -9.9999; //initialize with obviously wrong value
	if(g==female) p_maxPrtn= _nMaxCurrSexPartner_param_f[0]*exp(-_nMaxCurrSexPartner_param_f[1]*riskgroup);
	if(g==male) p_maxPrtn= _nMaxCurrSexPartner_param_m[0]*exp(-_nMaxCurrSexPartner_param_m[1]*riskgroup);
	
	return p_maxPrtn;
}




void Population::youthArrivals(double prd, bool save_trace_file)
{
	/// ARRIVALS OF YOUTHS STARTING TO BE SEXUALLY ACTIVE
	/// (EQUIVALENT TO A BITRH PROCESS)
	
	/// Arrival of youths is linked to 'true' birth rate
	/// and survival rates up to teenage
	
	int y = (int)_ageSexMin;
	
	// Reduced mortality rate after 5 years old
	// TO DO: do not hard code!
	double reducMortalityAfter5yrs = 0.2;
	
	double rate = _birthRate
	*(1-_infantMortality)
	*pow(1-_childMortality,4)
	*pow(1-_childMortality*reducMortalityAfter5yrs,y-5);
	
	unsigned long N = census_alive();
	
	// Expected number of new arrivals during the period
	double lambda = rate*prd*(double)(N);
	int nPrd = poisson(lambda);
	
	// Trace file
	if (save_trace_file){
		ofstream f(_DIR_OUT + "births.out",ios::app);
		f <<_simulationTime<< ","<< N <<","<<nPrd<< ","<<lambda<< "," << _birthRate<< "," << rate<<endl;
	}
	
	// Set features of young new comers
	
	for (int i=0; i<nPrd; i++){
		// Modelled population enters
		// at minimum age of sexual debut
		double age = _ageSexMin;
		
		// 50% chance female
		double u_gender = uniform01();
		Gender g = female;
		if (u_gender<0.5) g = male;
		
		// ==== Sexual Activity ====
		
		// Risk Group
		// is selected according
		// to the population-wide proportions
		vector<int>		rskgrp;
		vector<double>	proba_rskgrp;
		
		for (int i=0; i<=_maxRiskGroup; i++){
			rskgrp.push_back(i);
			proba_rskgrp.push_back(_propRiskGroup[i]);
		}
		
		int riskgroup = probaHistogramInt(rskgrp, proba_rskgrp);
		
		// Maximum number of concurrent sex partners:
		double p_maxPrtn = proba_nMaxCurrSexPartner(g,riskgroup);
		stopif((p_maxPrtn<=0 || p_maxPrtn>=1),
			   "'_nMaxCurrSexPartner_param' are not properly set because p_maxPrtn is <0 or >1");
		int maxSexPartner = 1+geometric(p_maxPrtn);
		stopif(maxSexPartner==0, "maxSexPartner==0");
		
		// No history of sexual activity
		int nLifetimePartners = 0;
		
		// Assumed not infected by any STIs
		vector<double>	stiDuration(_nSTImodelled,0.0);
		vector<bool>	stiSymptom(_nSTImodelled,false);
		
		// Assumed not widow nor divorced
		bool isWidow	= false;
		bool isDivorced	= false;
		
		// Circumcision
		// Chance of male being circumcised same as
		// proportion of already circumcised men in the population
		bool isCircum = false;
		if (g==male){
			double u_circum = uniform01();
			if (u_circum < _proportion_circum) isCircum=true;
		}
		
		Individual I(0, g, age,
					 maxSexPartner,
					 riskgroup,
					 nLifetimePartners,
					 isWidow,isDivorced, isCircum,
					 stiDuration,stiSymptom,
					 _STI,_RebHIV);
		
		stopif(maxSexPartner==0,"Creating an individual with maxSexPartner=0 but must be >=1");
		
		// Individual added to the population
		addIndividual(I);
	}
}


void Population::CSWrecruitment(double prd)
{
	/// Recruits commercial sex workers in the population
	
	// Retrieve the proportion of CSW in the population
	unsigned long Ncsw	= census_CSW().size();
	unsigned long N		= census_alive();
	double prop_csw		= (Ncsw==0)?0:((double)(Ncsw)/N);
	
	// Effective recruitment rate
	double a = _CSW_recruitment_saturPrm_1;
	double b = _CSW_recruitment_saturPrm_2;
	
	double effRate = _CSW_maxRecruitmentRate*(1+exp(-a*b))/(1+exp(a*(prop_csw-b)));
	
	// New CSWs recruited
	unsigned int NewCSW = poisson(effRate*N*prd);
	
	// STI prevalences
	vector<unsigned long> STIprevN = census_STIinfected();
	
	// Widow prevalence
	unsigned long Nw	= census_widow();
	unsigned long Ndiv	= census_divorced();
	
	
	// DEBUG
	//cout << " CSW added: "<< NewCSW << " (total:"<<Ncsw<<")"<<endl;
	
	
	// Add all new CSWs to the population
	
	for (int i=0; i<NewCSW; i++)
	{
		unsigned long uid=0;		// doesn't matter (taken care of with addIndividual)
		Gender g = female;			// assume all CSW are female
		int maxSexPartner = 999;	// large number
		
		// Age uniformly distributed b/w pre-specified limits
		double age = _CSW_recruitment_ageMin + (_CSW_recruitment_ageMax-_CSW_recruitment_ageMin)*uniform01();
		
		// Lifetime partners
		// TO DO: something smarter??
		// (doesn't realy matter as this number is expected to be large once CSW)
		int nLifetimePartners = poisson(0.5*(age-_ageSexMin));
		
		vector<double> stiDuration(_nSTImodelled,0.0);
		vector<bool> stiSymptom(_nSTImodelled,0.0);
		
		double epsilon = 1.0/365; // if the new CSW is STI positive, assume it's just been infected
		
		// STI infections for this new CSW
		// proba of being infected with STI 's'
		// is equal to population  prevalence
		
		for (int s=0; s<_nSTImodelled; s++)
		{
			// Bernoulli trial
			double proba = (double)(STIprevN[s])/N;
			int x = binom(proba, 1);
			if (x) stiDuration[s] = epsilon;
			
			// symptom (if infected)
			stiSymptom[s] = false;
			if(x) stiSymptom[s] = (bool)(binom(_STI[s].get_proba_symptomatic_female(), 1));
		}
		
		// Prevalence of widow/divorced the same as in general population
		bool isWidow	= binom((double)(Nw)/N,1)>0?true:false;
		bool isDivorced	= binom((double)(Ndiv)/N,1)>0?true:false;
		
		// Circumcision not relevant here, all CSW are females
		bool isCircum = false;
		
		// Create individual
		Individual I(uid,g,age,maxSexPartner,
					 _CSWriskGroup,nLifetimePartners,
					 isWidow, isDivorced,isCircum,
					 stiDuration,stiSymptom,
					 _STI,_RebHIV);
		
		// Individual added to the population
		addIndividual(I);
	}
}


void Population::CSWcessation(double prd)
{
	/// SELECT CSW THAT WILL CEASE PERFOMING CSW:
	/// i) randomly select the ones that will cease
	/// ii) CSW who are older than _CSW_ageMax
	
	vector<unsigned long> csw = census_CSW();
	unsigned long Ncsw = csw.size();
	
	// Number of CSW that will drop out this period
	unsigned long Nquit = binom(_CSW_cessationRate*prd, Ncsw);
	
	if (Nquit>0)
	{
		// Choose which CSW actually drops out
		vector<long> choosen = uniformIntVectorUnique(Nquit, 0, Ncsw-1);
		
		// Set the risk group to a non-CSW value (this is the only variable that defines CSW)
		
		for (int i=0; i<Nquit; i++)
		{
			int new_riskgrp = uniformInt(0, _maxRiskGroup);
			_individual[csw[choosen[i]]].set_riskGroup(new_riskgrp);
			
			// New max number of concurrent partners
			// is set to new value ((down from 999)
			// and makes sure it's not below current number of partners
			int nmp = 1+poisson(1+2.0*new_riskgrp);
			nmp = max(nmp,_individual[csw[choosen[i]]].get_nCurrSexPartner());
			
			_individual[csw[choosen[i]]].set_nMaxCurrSexPartner(nmp); // TO DO: do not hard-code
		}
	}
	// DEBUG
	//cout << "DEBUG: number of CSW quitting:"<<Nquit<<"/"<<Ncsw<< endl;
	
	// Refresh the list of CSW (now shorter because of cessation above)
	csw.clear();
	csw = census_CSW();
	Ncsw = csw.size();
	
	// All CSWs older than _CSW_ageMax are taken out
	// of CSW activity:
	for (int i=0; i<Ncsw; i++)
	{
		unsigned long uid_i = csw[i];
		if (_individual[uid_i].get_age()>=_CSW_ageMax)
		{
			int new_riskgrp = uniformInt(0, _maxRiskGroup);
			_individual[uid_i].set_riskGroup(new_riskgrp);
			
			// New max number of concurrent partners
			// is set to new value ((down from 999)
			// and makes sure it's not below current number of partners
			int nmp = 1+poisson(1+2.0*new_riskgrp);
			nmp = max(nmp,_individual[uid_i].get_nCurrSexPartner());
			
			_individual[uid_i].set_nMaxCurrSexPartner(nmp);
		}
	}
}


void Population::kill(unsigned long uid,bool save_trace_file)
{
	// Terminates all of dead's partnerships
	int np = _individual[uid].get_nCurrSexPartner();
	
	// Record all partners UID before death
	// (this MUST be done outside the loop below
	// because the member vector containing the
	// list of partners is updated when dissolvePartnership is called)
	
	vector<unsigned long> uid_partners;
	if (np>0) uid_partners = _individual[uid].getPartnerUID();
	
	for (int p=0; p < np; p++)
	{
		unsigned long uid_p = uid_partners[p];
		
		//DEBUG
		//cout << "UID "<<uid<<" is being killed: dissolve partnership with UID "<<uid_p<<endl;
		
		dissolvePartnership(uid,uid_p,save_trace_file);
		
		// TO DO : IF IT's A SPOUSAL ONLY, THEN WIDOW !!!!
		_individual[uid_p].set_isWidow(true); // The remaining partner is now "widow"
	}
	_individual[uid].kill();
	
	// Update the list of potential pregnant females
	if(_individual[uid].get_gender()==female &&
	   _individual[uid].get_age()<=_maxPregnantAge &&
	   !(_individual[uid].get_isPregnant()))
		_UID_pot_preg = popElementValue(_UID_pot_preg, uid);
}


double weibull_hazard(double t, double shape,double scale)
{
	return shape*scale*pow(scale*t,shape-1.0);
}


double Population::death_hazard(unsigned long uid)
{
	/// RETURNS THE HAZARD FuNCTION
	/// DISTRIBUTION: WEIBULL
	
	bool debug = false;
	
	double age = _individual[uid].get_age();
	double HIVduration;
	double hazard;		// Hazard function
	
	
	// Death is forced if above max age
	if (age >= _maxAge)	hazard = 999.99;
	
	if (age < _maxAge)
	{
		HIVduration = _individual[uid].get_STIduration(HIV);
		
		hazard = weibull_hazard(age,_deathParam_Weibull_shape,_deathParam_Weibull_scale);
		
		// If HIV pos, then additional hazard due to HIV progression (not age)
		if (HIVduration>0)
			hazard += weibull_hazard(HIVduration,_deathParam_Weibull_shape_hiv,_deathParam_Weibull_scale_hiv);
	}
	
	// DEBUG
	if (hazard>1 && debug)
	{
		cout<<"HAZ="<<hazard<< " ; d-age="<<age<<" ; hivdur="<<HIVduration<<" ; hazneg="<<weibull_hazard(age,_deathParam_Weibull_shape,_deathParam_Weibull_scale);
		cout<<" ; hazpos="<<weibull_hazard(HIVduration,_deathParam_Weibull_shape_hiv,_deathParam_Weibull_scale_hiv);
		cout <<" ; sh="<<_deathParam_Weibull_shape<<" ; sc= "<< _deathParam_Weibull_scale;
		cout<< " ; shiv= "<<_deathParam_Weibull_shape_hiv<< " ; schiv = "<< _deathParam_Weibull_scale_hiv<<endl;
	}
	
	return hazard;
}


void Population::save_death_hazard(string filename)
{
	double age_hiv_acquired = 25.0;
	double L		= _deathParam_Weibull_scale;
	double K		= _deathParam_Weibull_shape;
	double L_hiv	= _deathParam_Weibull_scale_hiv;
	double K_hiv	= _deathParam_Weibull_shape_hiv;
	
	
	ofstream f(filename.c_str());
	f<<"t,HIVneg,HIVpos,ageHIVacq"<<endl;
	
	for (double t=1; t<90; t+=0.02)
	{
		double hazneg = weibull_hazard(t,K,L);
		double hazpos = hazneg + weibull_hazard(max(0.0,t-age_hiv_acquired),K_hiv,L_hiv);
		f<<t<<","<<hazneg<<","<<hazpos<<","<< age_hiv_acquired << endl;
	}
	f.close();
}


void Population::deathEvents(double prd, bool save_trace_file)
{
	/// Scan the all individuals in the population
	/// based on their age and disease stage
	/// draws a random variable that will
	/// decide death or not, during the period "prd"
	
	// LOOP ON ALL INDIVIDUALS
	
	int cnt = 0;
	for (unsigned long uid=0; uid<_size; uid++)
	{
		if(_individual[uid].isAlive())
		{
			double u = uniform01();
			double probaDeath = death_hazard(uid)*prd;
			
			// kill this individual if stochastic event says so
			if (u < probaDeath){
				kill(uid,save_trace_file);
				cnt++;
			}
		}
	}
}



unsigned long Population::census_Females()
{
	unsigned long f = 0;
	
	for (int i=0; i<_size; i++){
		if (_individual[i].get_gender()==female && _individual[i].isAlive())
			f++;
	}
	
	return f;
}

unsigned long Population::census_Partnered(int numberOfPartners)
{
	unsigned long p = 0;
	
	for (int i=0; i<_size; i++)
	{
		if ( (_individual[i].get_nCurrSexPartner() == numberOfPartners) && (_individual[i].isAlive()) )
			p++;
	}
	
	return p;
}

unsigned long Population::census_singles()
{
	// CALCULATES THE NUMBER OF SINGLES IN THE POPULATION
	
	unsigned long s = 0;
	
	for (int i=0; i<_size; i++)
	{
		if ( (_individual[i].get_nCurrSexPartner()==0)
			&& (_individual[i].isAlive()))
			s++;
	}
	return s;
}

unsigned long Population::census_singles(Gender g)
{
	// CALCULATES THE NUMBER OF SINGLES IN THE POPULATION
	// SINGLE = NO SEX PARTNER (CURRENTLY)
	
	unsigned long s = 0;
	
	for (int i=0; i<_size; i++)
	{
		if ( (_individual[i].get_nCurrSexPartner()==0)
			&& (_individual[i].isAlive())
			&& (_individual[i].get_gender()==g))
			s++;
	}
	return s;
}

double Population::census_ratioSingles()
{
	/// CALCULATES THE RATIO OF SINGLES IN THE TOTAL POPULATION
	
	return (double)(census_singles())/(double)(census_alive());
}


double Population::census_ratioSingles(Gender g)
{
	/// CALCULATES THE RATIO OF SINGLES IN THE TOTAL POPULATION
	
	return (double)(census_singles(g))/(double)(census_alive(g));
}


unsigned long Population::census_widow()
{
	unsigned long p = 0;
	
	for (unsigned long i=0; i<_size; i++)
	{
		if (_individual[i].get_isWidow() && _individual[i].isAlive())
			p++;
	}
	return p;
}


unsigned long Population::census_divorced()
{
	unsigned long p = 0;
	
	for (unsigned long i=0; i<_size; i++)
	{
		if (_individual[i].get_isDivorced() && _individual[i].isAlive())
			p++;
	}
	return p;
}


vector<double> Population::census_getAgeGap(unsigned long uid)
{
	
	// CALCULATE THE AGE GAP BETWEEN A GIVEN INDIVIDUAL ('uid')
	// AND ALL HIS/HER SEX PARTNERS
	
	if (_individual[uid].get_nCurrSexPartner()==0)
	{
		cerr <<endl<< "ERROR [census_getAgeGap]: cannot calculate Age Gap with single individual (UID:"<<uid<<")"<<endl;
		exit(1);
	}
	
	vector<double> ag;
	
	for (int i=0; i<_individual[uid].get_nCurrSexPartner(); i++)
	{
		unsigned long partnerUID = _individual[uid].getPartnerUID(i);
		ag.push_back(getAgeGap(uid, partnerUID));
	}
	
	return ag;
}


vector<double> Population::census_AgeGap()
{
	
	// *** WARNING ***
	// Age gaps will appear TWICE !!!!
	
	// TO DO: FIX ME
	
	// CALCULATES THE DISTRIBUTION OF AGE GAPS AT THE POPULATION LEVEL
	
	vector<double> ag;
	
	for (int i=0; i<_size;i++)
	{
		//cout << i<<" ## "; //DEBUG
		
		if (_individual[i].isAlive() && _individual[i].get_nCurrSexPartner() > 0)
		{
			for (int p=0; p<_individual[i].get_nCurrSexPartner(); p++)
			{
				unsigned long uid = _individual[i].get_UID();
				unsigned long partnerUID = _individual[i].getPartnerUID(p);
				ag.push_back(getAgeGap(uid, partnerUID));
				
				//DEBUG
				//cout << uid << "--" << partnerUID << " : " <<getAgeGap(uid, partnerUID)<<" ";
			}
		}
		//cout << endl; //DEBUG
	}
	
	return ag;
}


unsigned long Population::census_alive()
{
	unsigned long res = 0;
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive())
			res++;
	}
	return res;
}

unsigned long Population::census_alive(Gender g)
{
	unsigned long res = 0;
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() && (_individual[uid].get_gender()==g))
			res++;
	}
	return res;
}

unsigned long Population::census_dead()
{
	unsigned long res = 0;
	for (int uid=0; uid<_size; uid++)
	{
		if (!_individual[uid].isAlive())
			res++;
	}
	return res;
}


unsigned long Population::census_circum()
{
	unsigned long res = 0;
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() && _individual[uid].get_isCircum())
			res++;
	}
	return res;
}


vector<unsigned long> Population::census_STIinfected()
{
	/// i^th element is the number of individuals infected with STI #i
	
	vector<unsigned long> res(_nSTImodelled, 0);
	
	for (int uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive()){
			for (int sti=0; sti<_nSTImodelled; sti++){
				if (_individual[uid].get_STIduration()[sti]>0)
					res[sti]++;
			}
		}
	}
	return res;
}


unsigned long Population::census_STIcoinfected(STIname s1, STIname s2){
	
	/// Number of individuals co-infected with s1 and s2
	
	unsigned long res = 0 ;
	
	for (int uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive() &&
			_individual[uid].get_STIduration(s1)>0 &&
			_individual[uid].get_STIduration(s2)>0)
			res++;
	}
	return res;
}


unsigned long Population::census_STIcoinfected(STIname s1, STIname s2, int riskgroup){
	
	/// Number of individuals from riskgroup co-infected with s1 and s2
	
	unsigned long res = 0 ;
	
	for (int uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive() &&
			_individual[uid].get_riskGroup()==riskgroup &&
			_individual[uid].get_STIduration(s1)>0 &&
			_individual[uid].get_STIduration(s2)>0)
			res++;
	}
	return res;
}


vector<unsigned long> Population::census_CSW()
{
	/// RETURNS UIDs OF ALL CSW (alive)
	
	vector<unsigned long> uid_csw(0);
	
	for (unsigned long uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive() &&
			_individual[uid].get_riskGroup()==_CSWriskGroup){
			uid_csw.push_back(uid);
		}
	}
	//DEBUG:
	//cout<<"risk group csw: "<<_CSWriskGroup;
	//displayVector(uid_csw);
	return uid_csw;
}


double Population::census_maleVisitCSW()
{
	unsigned long n_males = 0;
	unsigned long n_males_visitCSW = 0;
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() && (_individual[uid].get_gender()==male))
		{
			n_males++;
			if (_individual[uid].get_ever_visited_CSW()) n_males_visitCSW++;
		}
	}
	
	return (double)(n_males_visitCSW)/(double)(n_males);
}


double Population::census_maleVisitCSW(double maxDurationSinceLastVisit)
{
	unsigned long n_males = 0;
	unsigned long n_males_visitCSW = 0;
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() && (_individual[uid].get_gender()==male))
		{
			n_males++;
			double d = _individual[uid].get_lastVisitCSWDuration();
			
			// 'd>0' is a test that visit ever occured
			if (d<= maxDurationSinceLastVisit && d>0)
				n_males_visitCSW++;
		}
	}
	
	return (double)(n_males_visitCSW)/(double)(n_males);
}




vector<double> Population::census_ageDistribution(vector<double> ageBreaks)
{
	vector<double> ages;
	
	// retrieve the age of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive())
		{
			ages.push_back(_individual[i].get_age());
		}
	}
	return distributionNormalized(ages, ageBreaks);
}



vector<double> Population::census_ageGapDistribution(vector<double> bins)
{
	return distributionNormalized(census_AgeGap(), bins);
}



vector<double> Population::census_ageFirstSexDistribution(vector<double> ageBreaks,Gender thegender)
{
	vector<double> ages;
	
	// retrieve the age of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_ageFirstSex()>0 &&
			_individual[i].get_gender()==thegender)
		{
			ages.push_back(_individual[i].get_ageFirstSex());
		}
	}
	return distributionNormalized(ages, ageBreaks);
}



vector<double> Population::census_ageGapFirstSexSpouseDistribution(vector<double> yearsBreaks,
																   Gender thegender)
{
	vector<double> years;
	
	// retrieve the age of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_ageFirstSex()>0 &&
			_individual[i].get_ageFirstSpouse()>0 &&
			_individual[i].get_gender()==thegender)
		{
			years.push_back(_individual[i].get_ageFirstSpouse() - _individual[i].get_ageFirstSex());
		}
	}
	return distributionNormalized(years, yearsBreaks);
}


vector<double> Population::census_nLifeSexPrtnrDistrib(vector<double> nBreaks,
													   Gender g)
{
	vector<double> num;
	
	// retrieve the lifetime number sex partners
	// of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==g )
		{
			num.push_back(_individual[i].get_nLifetimePartner());
		}
	}
	return distributionNormalized(num, nBreaks);
}




void Population::census_nLifeSexPrtnrDistrib_todelete(Gender g)
{
	string s =_DIR_OUT+"check.out";
	
	ofstream f(s.c_str());
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==g )
		{
			f<<_individual[i].get_UID()<<","<<_individual[i].get_nLifetimePartner()<<endl;
		}
	}
}





vector<unsigned long> Population::census_nMaxCurrSexPartnerCount(vector<unsigned long> nBreaks,
																 Gender g)
{
	vector<unsigned long> num;
	
	// retrieve the max number sex partners
	// of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==g )
		{
			num.push_back(_individual[i].get_nMaxCurrSexPartner());
		}
	}
	return distribution(num, nBreaks);
	
}

vector<double> Population::census_nMaxCurrSexPartnerDistrib(vector<double> nBreaks,
															Gender g)
{
	vector<double> num;
	
	// retrieve the max number sex partners
	// of all living individuals
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==g )
		{
			num.push_back(_individual[i].get_nMaxCurrSexPartner());
		}
	}
	return distributionNormalized(num, nBreaks);
	
}




vector<double> Population::census_ageMalesVisitCSWDistribution(vector<double> ageBreaks)
{
	vector<double> ages;
	
	// retrieve the age of all males who ever visited a CSW
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==male &&
			_individual[i].get_ever_visited_CSW())
		{
			ages.push_back(_individual[i].get_age());
		}
	}
	return distributionNormalized(ages, ageBreaks);
}


vector<double> Population::census_ageMalesVisitCSWDistribution(vector<double> ageBreaks,
															   double maxDurationSinceLastVisit)
{
	vector<double> ages;
	
	// retrieve the age of all males who ever visited a CSW
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive() &&
			_individual[i].get_gender()==male &&
			_individual[i].get_lastVisitCSWDuration()>0 &&
			_individual[i].get_lastVisitCSWDuration()<maxDurationSinceLastVisit)
		{
			ages.push_back(_individual[i].get_age());
		}
	}
	return distributionNormalized(ages, ageBreaks);
}



vector<unsigned long> Population::census_riskGroup()
{
	vector<unsigned long> N(_maxRiskGroup+1, 0.0);
	unsigned long Ncsw = 0;
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive())
		{
			int rg = _individual[i].get_riskGroup();
			if (rg <  _CSWriskGroup) N[rg]++;
			if (rg == _CSWriskGroup) Ncsw++;
		}
	}
	
	// Add the number of CSW as the last element
	N.push_back(Ncsw);
	
	return N;
}


unsigned long Population::census_pregnant()
{
	/// Returns total number of all pregnant women
	
	unsigned long res;
	
	for (int i=0; i<_size; i++){
		if (_individual[i].isAlive() &&
			_individual[i].get_isPregnant()){
			res++;
		}
	}
	return res;
}

unsigned long Population::census_pregnant(int riskgroup)
{
	/// Returns total number of all pregnant women
	/// of a given risk group
	
	unsigned long res = 0;
	
	for (int i=0; i<_size; i++){
		if (_individual[i].isAlive() &&
			_individual[i].get_isPregnant() &&
			(_individual[i].get_riskGroup()==riskgroup) ){
			res++;
		}
	}
	return res;
}



vector<unsigned long> Population::census_pregnant_UID()
{
	/// Returns UID of all pregnant women
	
	vector<unsigned long> res;
	
	for (int i=0; i<_size; i++){
		if (_individual[i].isAlive() &&
			_individual[i].get_isPregnant()){
			res.push_back(_individual[i].get_UID());
		}
	}
	return res;
}



/* ****************************************** */
/* ************** ANALYTICS ***************** */
/* ****************************************** */





/* ****************************************** */
/* ************** PARTNERSHIPs ************** */
/* ****************************************** */



double Population::getPartnershipDuration(unsigned long uid1, unsigned long uid2)
{
	int p = _individual[uid1].getPartnershipPosition(uid2);
	
	return _individual[uid1].getPartnershipDuration()[p];
}


bool Population::alreadyPartners(unsigned long uid1, unsigned long uid2)
{
	// ASSESS IF UID1 AND UID2 ARE ALREADY PARTNERS
	
	vector<unsigned long> puid = _individual[uid1].getPartnerUID();
	
	bool already_partnered = false;
	
	// If no existing partner, then obviously not already partnered
	if (puid.size()==0) return false;
	
	// If existing partners, check if uid2 is in the list
	for (int i=0; i<puid.size() && !already_partnered; i++)
	{
		if (puid[i]==uid2) already_partnered=true;
	}
	
	return already_partnered;
}



void Population::formPartnership(unsigned long uid1, unsigned long uid2)
{
	/// FORMS A PARTNERSHIP BETWEEN UID1 AND UID2
	/// (IT'S A CASUAL BY DEFAULT, FOR ANY NEWLY FORMED PARTNERSHIP)
	
	// Integrity checks:
	
	string errmsg_dead = "Trying to partner with a dead individual ("+ to_string(uid1) + "or" + to_string(uid2) + " is dead)";
	string errmsg_g = "Trying to partner same ("+ to_string(uid1) + "or" + to_string(uid2) + ")";
	string errmsg_prtnnum = "Maximum and current number of partners inconsistent";
	
	_individual[uid1].addPartner(uid2);
	_individual[uid2].addPartner(uid1);
	
	stopif(!_individual[uid1].isAlive() || !_individual[uid2].isAlive(),errmsg_dead);
	stopif (_individual[uid1].get_gender()==_individual[uid2].get_gender(),errmsg_g);
	stopif(uid1>=_size || uid2>=_size,"UID outside population size!");
	stopif(_individual[uid1].get_nMaxCurrSexPartner()<_individual[uid1].get_nCurrSexPartner(),errmsg_prtnnum);
	stopif(_individual[uid2].get_nMaxCurrSexPartner()<_individual[uid2].get_nCurrSexPartner(),errmsg_prtnnum);
	
	// All partnerships start with one sex act:
	bool save_trace_file = false;
	add_sexAct(uid1, uid2, 1, save_trace_file);
	unsigned long uid_male = (_individual[uid1].get_gender()==male)?uid1:uid2;
	
	// type of first sex with this partner
	// (defined only for males)
	double u = uniform01();
	vector<unsigned int> sextype(3,0.0);
	// TO DO: this vector must depend on risk group
	if (u<0.33) sextype[0]=1;
	if (u>=0.33 && u<0.66) sextype[1]=1;
	if (u>0.66) sextype[2]=1;
	_individual[uid_male].set_UID_n_sexAct_Type_period(sextype);
	
	// Update '_partnershipsMatrix'
	vector<double> tmp(2);
	tmp[0] = uid1;
	tmp[1] = uid2;
	if (_individual[uid1].get_gender()==male){
		tmp[0] = uid2;
		tmp[1] = uid1;
	}
	_partnershipsMatrix.addRowVector(tmp);
	
	// update the total number of partnerships
	_totalNumberPartnerships = _partnershipsMatrix.getNbRows();
}



void Population::dissolvePartnership(unsigned long uid1, unsigned long uid2,
									 bool save_to_file)
{
	/// Dissolve the partnership on BOTH sides
	/// (only need to run it once)
	/// and updates the partnership matrix _PartnerMatrix
	
	Individual I1 = _individual[uid1];
	Individual I2 = _individual[uid2];
	
	if (!I1.isAlive() || !I2.isAlive())
	{
		cout << "ERROR [dissolvePartnership]: trying to dissolve with a dead individual ("
		<<uid1<<"or"<<uid2<<" is dead)"<<endl;
		exit(1);
	}
	
	I1.erasePartner(uid2);
	I2.erasePartner(uid1);
	
	// Update partnership matrix
	
	Gender g1 = I1.get_gender();
	
	//unsigned long i1, i2;
	
	if (g1==female)
	{
		unsigned long pos = findPositionIn_partnershipsMatrix(uid1,uid2);
		_partnershipsMatrix.removeRow(pos);
		
		//DEBUG
		//cout << ">> DISSOLVE: ("<<uid1<<","<<uid2<<") @ pos="<< pos <<endl;
		
		// GET RID OF ALL THIS ----
		//i1 = findIndexElement(_all_UID_femaleUID, uid1);
		//i2 = findIndexElement(_all_UID_maleUID, uid2);
		//_PartnerMatrix(i1,i2) = 0;
		// -----------------
	}
	if (g1==male)
	{
		unsigned long pos = findPositionIn_partnershipsMatrix(uid2,uid1);
		_partnershipsMatrix.removeRow(pos);
	}
	
	
	_individual[uid1] = I1;
	_individual[uid2] = I2;
	
	// update the total number of partnerships
	_totalNumberPartnerships = _partnershipsMatrix.getNbRows();
	
	
	// Trace files
	if(save_to_file){
		ofstream f(_DIR_OUT + "dissolve.out",ios::app);
		f<<_simulationTime<<","<<I1.get_riskGroup()<<","<<I2.get_riskGroup()<<endl;
	}
}


double Population::meanAgePartners(unsigned long uid)		// Calculate mean age of all partners
{
	double res=0;
	
	int N = _individual[uid].get_nCurrSexPartner();
	
	assert(N>0);
	
	for (int i=0; i<N; i++)
	{
		// Retrieves UID of ith partner of this individual
		unsigned long uid_p = _individual[uid].getPartnerUID(i);
		
		// Retrieves age of the ith partner
		double age = _individual[uid_p].get_age();
		
		res += age;
		
		// DEBUG
		//cout<<"AGE P"<<i<<"/"<<N-1<<" = "<<age<<endl;
	}
	
	if (N>0) res = res/N;
	
	return res;
}



unsigned int Population::getMaxDegree()
{
	unsigned int m = 0;
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive())
		{
			int tmp = _individual[i].get_nCurrSexPartner();
			if (tmp>m) m= tmp;
		}
	}
	
	return m;
}

// ===== Matching Functions =====



void Population::update_UID_PartnershipMatrix(Individual I)
{
	// UPDATE _PartnerMatrix AFTER ONE INDIVIDUAL IS ADDED
	
	// retrieves the UIDs of all males and females
	// (including the ones that are recently added)
	
	Gender g = I.get_gender();
	unsigned long uid = I.get_UID();
	
	if (g==female)
	{
		_all_UID_femaleUID.push_back(uid);
		
		int n = _PartnerMatrix.getNbCols();
		//vector<double> v(n,0.0);
		//_PartnerMatrix.addRowVector(v); // too slow need to fix dcMatrix function...
		
		int nr = _PartnerMatrix.getNbRows()+1;
		_PartnerMatrix.resize(nr, n);
	}
	if (g==male)
	{
		_all_UID_maleUID.push_back(uid);
		
		int n = _PartnerMatrix.getNbRows();
		vector<double> v(n,0.0);
		_PartnerMatrix.addColVector(v);// too slow need to fix dcMatrix function...
		
		/*THIS DOESN't WORK
		 BECAUSE NOT AS EASY AS WITH ROW
		 NEED TO SWAP ELEMENTS OF 'val' MATRIX CLASS MEMBER
		 
		 int nc = _PartnerMatrix.getNbCols();
		 _PartnerMatrix.resize(n, nc+1);
		 _PartnerMatrix.setColumnValues(nc, v);*/
		
	}
}



double Population::form_age(double age_f)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON FEMALE AGE
	
	double tmp;
	
	if(age_f<=_formation_age_fullstart){
		// Linear component for young ages
		tmp = (age_f-_ageSexMin)/(_formation_age_fullstart-_ageSexMin);
	}
	else{
		// Logistic component for older ages
		tmp = 1.0/(1+exp(_formation_age_shape*(age_f-_formation_age_pivot)));
	}
	
	return _formation_age_fmin + (1-_formation_age_fmin)*tmp;;
}

double Population::form_age_male(double age_m)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON MALE'S AGE
	
	double tmp;
	
	if(age_m<=_formation_age_fullstart)	{
		// Linear component for young ages
		tmp = (age_m-_ageSexMin)/(_formation_age_fullstart-_ageSexMin);
	}
	else tmp = 1.0;
	
	return _formation_age_fmin + (1-_formation_age_fmin)*tmp;;
}



double Population::form_ageGap(double ageGap)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON PARTNERSHIP'S AGE GAP
	
	double gapmin	= _formation_agegap_shape[0];
	double gapavg	= _formation_agegap_mean;
	double a		= _formation_agegap_shape[1];
	double d		= _formation_agegap_shape[2];
	
	double b		= a/(d*pow(gapavg-gapmin,d));
	double fmax		= pow(a/(b*d*exp(1)),(a/d));
	double xx		= ageGap - gapmin;
	double tmp		= 0.0;
	if (xx>0) tmp = pow(xx,a)*exp(-b*pow(xx,d))/fmax;
	
	return _formation_agegap_fmin + (1-_formation_agegap_fmin)*tmp;
}

double Population::form_riskGroup(int risk_f, int risk_m)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON RISK GROUPS
	
	double a0 = _formation_RiskGroup[0];
	double a1 = _formation_RiskGroup[1];
	double rmax = _maxRiskGroup;
	return exp(-a0*(2*rmax-(risk_f+risk_m))-a1*fabs(risk_f-risk_m) );
}


double Population::form_deficit(double def_f, double def_m)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON PARTNERSHIP DEFICIT
	
	double q = _formation_PartnDeficit[0];
	return pow(def_f*def_m,q);
}

double Population::form_STIsymptom(unsigned long uid_f,unsigned long uid_m)
{
	/// RETURNS PROBABILITY LEVEL
	/// FOR PARTNERSHIP FORMATION
	/// BASED ON STI SYMPTOMS
	
	double a_f = _individual[uid_f].is_symptomatic()? _formation_STIsymptom_f:1.0;
	double a_m = _individual[uid_m].is_symptomatic()? _formation_STIsymptom_m:1.0;
	return a_f*a_m;
}


bool Population::isFormationPossible(unsigned long uid1,
									 unsigned long uid2,
									 bool save_trace_file)
{
	/// DETERMINES IF PARTNERSHIP FORMATION
	/// IS POSSIBLE BASED ON ALL FACTORS
	
//	ofstream tracefile(_DIR_OUT + "partnership_tentatives.out",ios::app);
	
	// Check both are alive
	bool alive = _individual[uid1].isAlive() && _individual[uid2].isAlive();
	string errmsg;
	errmsg = " Trying to form partnership with at least a dead individual";
	stopif(!alive, errmsg);
	
	// Check if genders are different
	Gender g1 = _individual[uid1].get_gender();
	Gender g2 = _individual[uid2].get_gender();
	bool genderCheck = (g1 != g2);
	errmsg="trying to form partnership with same gender individuals";
	stopif(!genderCheck,errmsg);

	// Check they are not already in partnership
	bool alreadyPartner = alreadyPartners(uid1, uid2);
	if (alreadyPartner) return false;

	// if individual is not open to new partnership, then impossible:
	bool open1 = _individual[uid1].isOpenToNewPartnership();
	bool open2 = _individual[uid2].isOpenToNewPartnership();
	if (!open1 || !open2){
//		if(save_trace_file) tracefile << "not_open,,,,0"<<endl;
		return false;
	}
	
	// === Age preferences ===
	
	double ageGap = getAgeGap(uid1,uid2);
	double age_female=0.0;
	double age_male=0.0;
	if (g1==female) {age_female = _individual[uid1].get_age(); age_male = _individual[uid2].get_age();}
	if (g2==female) {age_female = _individual[uid2].get_age(); age_male = _individual[uid1].get_age();}
	
	double f_age		= form_age(age_female);
	double f_age_male	= form_age_male(age_male);
	double f_ageGap		= form_ageGap(ageGap);
	
	// === Risk Group ===
	
	double f_riskgroup = form_riskGroup(_individual[uid1].get_riskGroup(),
										_individual[uid2].get_riskGroup());
	
	// === Partnership deficit ===
	
	int n1 = _individual[uid1].get_nCurrSexPartner();
	int n2 = _individual[uid2].get_nCurrSexPartner();
	int m1 = _individual[uid1].get_nMaxCurrSexPartner();
	int m2 = _individual[uid2].get_nMaxCurrSexPartner();
	
	double d1 = (double)(m1-n1)/m1;
	double d2 = (double)(m2-n2)/m2;
	double f_deficit = form_deficit(d1, d2);
	
	// === STI symptom ===
	
	unsigned long uid_f = _individual[uid1].get_gender()==female? uid1:uid2;
	unsigned long uid_m = (uid_f==uid1)?uid2:uid1;
	double f_STIsymptom = form_STIsymptom(uid_f, uid_m);

	// Bernoulli probability for formation event
	double proba = f_age_male * f_age * f_ageGap * f_riskgroup * f_deficit* f_STIsymptom;
	double u = uniform01();
	
	if(save_trace_file)
	{
//		tracefile << f_age<<","<< f_ageGap<<","<< f_riskgroup<<",";
//		tracefile << f_deficit<<","<< f_STIsymptom<<",";
//		tracefile << (u<proba)<<endl;
	}
	
	// Based on draw, decide wether formation occurs or not
	return u<proba? true : false;
}



bool Population::spousalProgression(unsigned long uid1, unsigned long uid2)
{
	/// DETERMINES IF UID1 AND UID2 CAN
	/// BECOME SPOUSES (THEY ARE ALREADY CASUAL PARTNERS)
	
	bool debug = false;
	
	// Determines who is the female
	unsigned long uid_f = uid2;
	unsigned long uid_m = uid1;
	if (_individual[uid1].get_gender()==female){
		uid_f = uid1;
		uid_m = uid2;
	}
	
	// If the female is already in a spousal union,
	// she cannot be the spouse of another male (polygyny).
	// Hence we stop here.
	if (_individual[uid_f].get_nCurrSpouse()>0)	return false;
	
	// If the female is CSW, no spousal partnership allowed
	if (_individual[uid_f].get_riskGroup()==_CSWriskGroup) return false;
	
	// Retrieves the age of the female
	double Af = _individual[uid_f].get_age();
	double gap = getAgeGap(uid1,uid2);
	
	double s1 = expGauss(Af,
						 _spousalProgress_meanAge_f,
						 _spousalProgress_varAge_f);
	double s2 = expGauss(gap,
						 _spousalProgress_meanGap,
						 _spousalProgress_varGap);
	double s = s1*s2;
	
	// Duration of this casual partnership
	int ip		= _individual[uid1].getPartnershipPosition(uid2);
	double tau	= _individual[uid1].getPartnershipDuration()[ip];
	
	double d_tau = expGauss(tau,
							_spousalProgress_durationK1,
							_spousalProgress_durationK2);
	
	// If male has already other spouses, then age gap of
	// new candidate spouse is compared to other spouses'
	int nSpouses = _individual[uid_m].get_nCurrSpouse();
	
	double K=0;
	int hasOtherSpouse = nSpouses>0 ? 1 : 0;
	
	if (hasOtherSpouse){
		// Retrieve all age gaps with spouses only
		vector<unsigned long> sp_uid = _individual[uid_m].getSpouseUID();
		unsigned long nSp = sp_uid.size();
		vector<double> ageGaps;
		
		for (int is=0; is<nSp; is++){
			double a = getAgeGap(uid_m, sp_uid[is]);
			ageGaps.push_back(a);
		}
		
		// Calculate the smallest age gaps with existing spouses
		double ageGapMin = minElementVector(ageGaps);
		
		// Define the difference b/w the new comer spouse's age gap
		// and the smallest age gap of existing spouses
		double delta = ageGapMin - getAgeGap(uid_m, uid_f);
		
		K = expGauss(delta, _spousalProgress_meanDiffAgeGap, _spousalProgress_varDiffAgeGap);
	}
	
	
	// Draw random variable that will
	// determine spousal progression or not
	
	bool spouseProgression = false;
	double proba = s*( (1-hasOtherSpouse) + hasOtherSpouse*K )*d_tau * _spousalProgress_maxRate;
	double u = uniform01();
	
	if (u<proba) spouseProgression = true;
	
	// If first time spousal partnership,
	// then record age
	if (spouseProgression){
		if (_individual[uid1].get_ageFirstSpouse()<0) _individual[uid1].set_ageFirstSpouse(_individual[uid1].get_age());
		if (_individual[uid2].get_ageFirstSpouse()<0) _individual[uid2].set_ageFirstSpouse(_individual[uid2].get_age());
	}
	
	return spouseProgression;
}


void Population::spousalScanAllCasual()
{
	/// SCAN ALL CASUAL PARTNERSHIPS AND
	/// CONSIDER IF SPOUSAL PROGRESSION HAPPENS
	
	
	// If no partneships, exits
	if (_totalNumberPartnerships==0) return;
	
	
	// Retrieves all partnerships (casual and spousal)
	dcMatrix P = getPartnershipsUID();
	
	int np = P.getNbRows();
	
	// Go through all partnerships and
	// test if spousal progression possible
	
	for (int i=0; i<np; i++)
	{
		unsigned long uid1 = P(i,0);
		unsigned long uid2 = P(i,1);
		
		int ip = _individual[uid1].getPartnershipPosition(uid2);
		
		// Check if partnership not already a spousal one
		
		// DEBUG
		//cout << endl << uid1 << " Spousal before if condition: " << _individual[uid1].get_isSpousalPartner()[ip] <<endl;
		
		if (! _individual[uid1].get_isSpousalPartner()[ip] )
		{
			bool progressToSpousal = spousalProgression(uid1,uid2);
			
			// If spousal progression successful, then update relevant variables
			if (progressToSpousal)
			{
				_individual[uid1].convertToSpouse(uid2);
				_individual[uid2].convertToSpouse(uid1);
			}
			
			// DEBUG:
			//cout << uid1 <<"<s?s>"<<uid2 << " : " << progressToSpousal <<endl;
		}
	}
	
	// Update the records of spouses in the population
	dcMatrix dummy = getSpousesUID();
}



unsigned long Population::nUnmatchedPartnerships()
{
	// RETURNS THE NUMBER OF POTENTIAL PARTNERSHIPS
	// THAT ARE NOT FULLFILLED
	
	unsigned long res = 0;
	
	for (unsigned long p=0; p<_size; p++)
	{
		if (_individual[p].isAlive())
		{
			int n_i = _individual[p].get_nCurrSexPartner();
			int m_i = _individual[p].get_nMaxCurrSexPartner();
			
			if (n_i < m_i && m_i>0)
				res += m_i-n_i;
		}
	}
	return res;
}





void Population::formOnePartnershipFromFemale_rand(unsigned long uid_fem,
												   bool save_trace_file)
{
	/// Forms a partnership with a randomly selected male given a female 'uid_fem'.
	/// Check if match is possible given each both individuals' features.
	
	bool debug = false;
	
	// Retrieve available males only
	vector<unsigned long>  males_available = getUID_available(male);
	
	unsigned long nMa = males_available.size();
	
	// If no male is available, then finished
	if (nMa==0) return;
	
	// Randomly select a male
	int i_m = uniformInt(0, nMa-1);
	unsigned long uid_male = males_available[i_m];
	
	stopif(uid_male>=_size,"problem finding an available male");
	
	// Determine if formation is possible
	bool matched = isFormationPossible(uid_fem, uid_male,save_trace_file);
	
	// --- DEBUG ---
	if (debug){
		cout << endl << uid_fem << "<?>"<<uid_male<< " : ";
		cout << matched;
		cout << " ; maxPartn=" << _individual[uid_fem].get_nMaxCurrSexPartner();
		//displayVector(females_availble);
		//displayVector(males_availble);
	}
	
	if (matched) formPartnership(uid_fem,uid_male);
}


vector<unsigned long> Population::formUIDcandidateFemale(double prd)
{
	// Retrieve available females only
	vector<unsigned long>  uid_females_available = getUID_available(female);
	
	unsigned long nFa = uid_females_available.size();
	
	// If no female is available, then finished
	vector<unsigned long> nothing(0);
	if (nFa==0) return nothing;
	
	// Draw the random number of females candidate
	// for partnership formation
	unsigned long nFc = binom(_formation_MaxRate*prd, nFa);
	
	// Select which index will be candidate
	vector<long> random_index= uniformIntVectorUnique(nFc,0,nFa-1);
	
	// Retrieve the UIDs of the selected females
	vector<unsigned long> uid_females_candidate(nFc);
	
	for (int i=0; i<nFc; i++)
	{
		uid_females_candidate[i] = uid_females_available[random_index[i]];
	}
	
	return uid_females_candidate;
}


void Population::formPartnerships(double prd, bool save_trace_file)
{
	
//	ofstream tracefile(_DIR_OUT + "nFemalesCandidate.out",ios::app);
	
	// Retrieve the UIDs of all females candidate for new partnership
	vector<unsigned long> uid_females_candidate = formUIDcandidateFemale(prd);
	
	// Number of those females
	unsigned long nFc = uid_females_candidate.size();
	
//	if(save_trace_file) tracefile<<_simulationTime<<","<<nFc<<endl;
	
	// For each of those females, pick randomly a male and try to match
	for (int i=0; i<nFc; i++){
		formOnePartnershipFromFemale_rand(uid_females_candidate[i],save_trace_file);
	}
	
	// Calculate the new total number of partnerships
	_totalNumberPartnerships = totalNumberPartnerships();
	
}




void Population::dissolvePartnerships(double prd,
									  bool save_trace_file){
	
	/// DISSOLVE PARTNERSHIPS
	
	// prd = period considered, in years (e.g. 1 day=1/365)
	
	double proportion = prd*_dissolution_MaxRate;
	
	unsigned long N = totalNumberPartnerships();
	
	// Number of partnerships candidate for dissolution
	unsigned long Dc = binom(proportion, N);
//    std::binomial_distribution<unsigned long> binom(N,proportion);
//    unsigned long Dc = binom(_RANDOM_GENERATOR);
	
	// Retrieve all partnerships
	dcMatrix P = getPartnershipsUID();
	
	// Select which partnerships will be candidate for dissolution
	vector<long> random_index = uniformIntVectorUnique(Dc,0,N-1);
	
	// Scan all selected partnerships
	// and draw a random variable for each of them
	// to determine dissolution or not:
	for (int i=0; i<Dc; i++)
	{
		// By convention, 1st col is female's UID
		unsigned long uid_f = P(random_index[i],0);
		unsigned long uid_m = P(random_index[i],1);
		
		Individual I_f = _individual[uid_f];
		Individual I_m = _individual[uid_m];
		
		double g_Spouse		= dissolve_spouse(uid_f,uid_m);
		double g_RG			= dissolve_riskGroup(uid_f,uid_m);
		double g_Duration	= dissolve_duration(uid_f,uid_m);
		double g_Age		= dissolve_age(uid_f, uid_m);
		double g_Deficit	= dissolve_partnerDeficit(uid_f,uid_m);
		double g_symptom	= dissolve_STIsymptoms(uid_f,uid_m);
		//double g_ageConc	= dissolve_ageConcPartn(uid_f, uid_m);
		
		// Probability to dissolve this partnership
		// based on its features:
		double proba = g_Spouse * g_RG * g_Duration * g_Age * g_Deficit * g_symptom ;//* g_ageConc;
		
		// Draw the Bernoulli random variable
		// to determine dissolution, or not
		bool successDissolution = false;
		
		if (uniform01() < proba){
			dissolvePartnership(uid_f,uid_m, save_trace_file);
			successDissolution = true;
		}
	}
}



unsigned long Population::totalNumberPartnerships()
{
	return _partnershipsMatrix.getNbRows();
}


unsigned long Population::findPositionIn_partnershipsMatrix(unsigned long uid_fem,
															unsigned long uid_male)
{
	unsigned long i=0;
	
	while (_partnershipsMatrix(i,0) != uid_fem ||
		   _partnershipsMatrix(i,1) != uid_male)
	{
		i++;
		
		if (i>=_partnershipsMatrix.getNbRows()){
			cout << "ERROR [Population::findPositionIn_partnershipsMatrix]: partnership not found"<<endl;
			exit(1);
		}
	}
	
	return i;
}


dcMatrix	Population::getPartnershipsUID()
{
	// GET RID OF THIS FUNCTION
	return get_partnershipsMatrix();
}

dcMatrix	Population::getSpousesUID()
{
	if (_totalNumberPartnerships==0)
	{
		cout << endl << "*** WARNING [Population::getSpousesUID]: No partnerships, so no spousal ones...";
		dcMatrix resZero(0);
		return resZero;
	}
	
	dcMatrix res(1,2);
	
	unsigned long c=0;
	
	for (int i=0; i<_totalNumberPartnerships; i++)
	{
		// Check if spousal union
		unsigned long uid_f = _partnershipsMatrix(i,0);
		unsigned long uid_m = _partnershipsMatrix(i,1);
		
		if (_individual[uid_f].isSpouse(uid_m))
		{
			if (c==0)
			{
				res(c,0) = uid_f;
				res(c,1) = uid_m;
			}
			else
			{
				vector<double> x(2);
				x[0] = uid_f;
				x[1] = uid_m;
				res.addRowVector(x);
			}
			
			c++;
		}
	}
	
	_totalNumberSpousalPartnerships = c;
	
	if (c==0)
	{
		dcMatrix dummy(0,0);
		res = dummy;
	}

	return res;
}


double Population::get_STIsusceptFactor(STIname stiname, unsigned long uid){
	return(_individual[uid].get_STIsusceptFactor(stiname));
}



double Population::dissolve_spouse(unsigned long uid1, unsigned long uid2)
{
	/// PROBABILITY TO DISSOLVE - SPOUSAL COMPONENT
	
	double proba = 1.0;
	if (_individual[uid1].isSpouse(uid2)) proba = _dissolution_spouse;
	
	return proba;
}


double Population::dissolve_riskGroup_fct(int r1, int r2){
	double shape = _dissolution_RiskGroup[0];
	return exp(shape*(double)(r1+r2 - 2*_maxRiskGroup));
}

double Population::dissolve_riskGroup(unsigned long uid1, unsigned long uid2)
{
	/// PROBABILITY TO DISSOLVE - RISK GROUP COMPONENT
	
	int r1 = _individual[uid1].get_riskGroup();
	int r2 = _individual[uid2].get_riskGroup();
	
	// If partner is CSW, then no reduction of proba to dissolve
	double res = 1.0;
	
	if (_CSWriskGroup!=r1 && _CSWriskGroup!=r2)
		res = dissolve_riskGroup_fct(r1,r2);
	
	return res;
}


double Population::dissolve_duration_fct(double prtnr_duration){
	double a = _dissolution_duration[0];
	double b = _dissolution_duration[1];
	double c = _dissolution_duration[2];
	return a+b*exp(-c*prtnr_duration);
}


double Population::dissolve_duration(unsigned long uid1, unsigned long uid2){
	/// PROBABILITY TO DISSOLVE - PARTNERSHIP DURATION COMPONENT
	return dissolve_duration_fct(getPartnershipDuration(uid1, uid2));
}


double Population::dissolve_age_fct(double age1, double age2){
	double d1	= _dissolution_age[0];
	double d2	= _dissolution_age[1];
	double pmin	= _dissolution_age[2];
	return pmin + (1-pmin)/(1+exp(d2*(age1+age2-d1)));
}


double Population::dissolve_age(unsigned long uid1, unsigned long uid2)
{
	/// PROBABILITY OF DISSOLUTION - AGE COMPONENT
	
	double a_1 = _individual[uid1].get_age();
	double a_2 = _individual[uid2].get_age();
	
	return dissolve_age_fct(a_1,a_2);
}


double Population::dissolve_partnerDeficit_fct(double d1, double d2){
	double q	= _dissolution_PartnerDeficit[0];
	double pmin	= _dissolution_PartnerDeficit[1];
	return pmin + (1-pmin)*pow( (1-d1)*(1-d2), q);
}


double Population::dissolve_partnerDeficit(unsigned long uid1, unsigned long uid2)
{
	/// PROBABILITY OF DISSOLUTION - PARTNER DEFICIT COMPONENT
	
	int n1 = _individual[uid1].get_nCurrSexPartner();
	int n2 = _individual[uid2].get_nCurrSexPartner();
	int m1 = _individual[uid1].get_nMaxCurrSexPartner();
	int m2 = _individual[uid2].get_nMaxCurrSexPartner();
	
	double d1 = (double)(m1-n1)/m1;
	double d2 = (double)(m2-n2)/m2;
	
	// Integrity check
	string errmsg =  "more concurrent partners than allowed";
	if (d1<0) {cout<<"UID="<<uid1<<"m1="<<m1<<";n1="<<n1<<endl;_individual[uid1].displayInfo();}
	if (d2<0) cout<<"UID="<<uid2<<"m2="<<m2<<";n2="<<n2<<endl;
	
	stopif(d1*d2<0, errmsg);
	// ----
	
	return dissolve_partnerDeficit_fct(d1,d2);
}



double Population::dissolve_STIsymptoms(unsigned long uid1, unsigned long uid2)
{
	/// PROBABILITY OF DISSOLUTION - STI SYMPTOMS INFLUENCE
	/// IF ANY OF THE TWO PARTNERS HAS SYMPTOMATIC STI INFECTION
	
	
	double res = 1.0;
	
	if (_individual[uid1].is_symptomatic() || _individual[uid2].is_symptomatic())
		res = _dissolution_STI_symptom;
	
	return res;
}



// ======================================================================



unsigned long Population::chooseRandomCSW()
{
	/// CHOOSE ONE CSW AMONG ALL CSWs IN THE POPULATION
	
	vector<unsigned long> uid_csw = census_CSW();
	int n_csw = uid_csw.size();
	
	unsigned long cswChosen = uniformInt(0, n_csw-1);
	return uid_csw[cswChosen];
}


void Population::updateAllDurations(double timeStep)
{
	// As time elapses, all durations increase by the time step
	
	for (int i=0; i<_size; i++)
	{
		bool isAlive = _individual[i].isAlive();
		bool isSingle = _individual[i].isSingle();
		
		// -- SingleDuration
		
		if (isAlive && isSingle) _individual[i].increaseSingleDuration(timeStep);
		
		// -- Partnerships durations
		
		if (isAlive && !isSingle)
		{
			int Np = _individual[i].get_nCurrSexPartner();
			
			for (int p=0; p<Np; p++)
			{
				_individual[i].increasePartnershipDuration(p,timeStep);
			}
		}
		
		// -- Diseases & Treatments durations
		
		if (isAlive)
		{
			_individual[i].increaseSTIdurations(timeStep);
			_individual[i].increase_STItreatDuration(timeStep);
			_individual[i].increase_lastVisitCSWDuration(timeStep);
		}
	}
}


void Population::updateAllAges(double timeStep)
{
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive()) {
			
			_individual[i].increase_age(timeStep);
			
			// If female gets too old to become pregnant, remove her from
			// the list of potential pregnant females
			double a = _individual[i].get_age();
			
			if (_individual[i].get_gender()==female &&
				!(_individual[i].get_isPregnant()) &&
				a >_maxPregnantAge &&
				a-timeStep <_maxPregnantAge) // <-- odd, but necessary to deal with the very first update
			{
				_UID_pot_preg = popElementValue(_UID_pot_preg, _individual[i].get_UID());
			}
		}
	}
}


/* **************************************************
 ************    GENETICS                ************
 *****************************************************/

void Population::update_genomes_mutations(double timestep){
    
    // loop through all STIs:
    for(int s=0; s<_STI.size(); s++){
        
        double gt       = _STI[s].get_pathogen_generationTime();
        double sng      = _STI[s].get_subst_perNucleotide_perGeneration();
        unsigned long N = _STI[s].get_genome().length();
        STIname stiname = _STI[s].get_name();
        
        double mean_n_mutations = sng * N * timestep / gt;
        std::poisson_distribution<unsigned long> poiss(mean_n_mutations);
        
        // loop through all infected individuals:
        for(unsigned long i=0; i<_size; i++){
            

            if (_individual[i].get_STIduration(stiname)>0){
                
                if(i==17){
                    double dummy=0;
                }
                
                
                // Draw the number of mutations:
                unsigned long n_mutations = poiss(_RANDOM_GENERATOR);
                if(n_mutations>0){
                    // Performs mutations:
                    string mutated = mutate_genome(_individual[i].get_genome(s), n_mutations);
                    _individual[i].set_genome(mutated, s);
                }
                // DEBUG
                //cout << STInameString(stiname) << " indiv_"<<i<<": "<<n_mutations<<" mutations."<<endl;
            }
        }
    }
}


void Population::update_genomes_clear(){
    /// Clear genome information if infection is cleared.
        
    // loop through all STIs:
    for(int s=0; s<_STI.size(); s++){
        STIname stiname = _STI[s].get_name();
        // loop through all infected individuals:
        for(unsigned long i=0; i<_size; i++){
            if (_individual[i].get_STIduration(stiname)==0){
                _individual[i].set_genome("", s);
            }
        }
    }
}

vector<string> Population::get_population_genomes(int sti_pos){
    
    vector<string> g;
    
    STIname stiname = _STI[sti_pos].get_name();
    // loop through all infected individuals:
    for(unsigned long i=0; i<_size; i++){
        if (_individual[i].get_STIduration(stiname)>0){
            g.push_back( _individual[i].get_genome(sti_pos) );
        }
    }
    return g;
}



/* **************************************************
 ************    SEXUAL ACTIVITY ************
 *****************************************************/



void Population::sexActs_number_males(unsigned long uid,
									  double period,
									  bool save_trace_file)
{
	/// Calculates the number of sex acts
	/// perfomed by this male during a given period
	
	// Males only. Number for female will be calculated from distribution cascade
	string errmsg = "individual"+ to_string(uid)+ " must be a male, but is not!";
	stopif (_individual[uid].get_gender()==female, errmsg);
	
	// Rate of sexual activity initiated with highest value
	double	Rsex		= _sexAct_maxRate_male;
	int		nSexActs	= 0;
	
	// Variables used by both single and partnered males
	
	// Reduction due to age factor:
	double age_m = _individual[uid].get_age();
	double h_age = sexAct_reduce_age(age_m);
	// Reduction due to risk group:
	double rg_m		= _individual[uid].get_riskGroup();
	double h_risk	= sexAct_reduce_riskGroup(rg_m);
	
	// Reduction due to STI symptoms
	vector<bool> stiSympt	= _individual[uid].get_STIsymptom();
	int anySTIsympt			= sumElements(stiSympt);
	bool anySTIsymptBool	= (anySTIsympt>0?true:false);
	double h_sti			= sexAct_reduce_STIsymptom(anySTIsymptBool, _individual[uid].get_gender());
	
	// Reduction due to number of partners
	int nPartners		= _individual[uid].get_nCurrSexPartner();
	double h_nPartners	= sexAct_reduce_nPartners(nPartners);
	
	// Reduction due to number of children taken care of
	double h_child = 1.0; // NOT USED --- TO DO: implement when Children implemented
	
	// Global reduction factor on male's side
	double Rm = h_age * h_risk * h_sti  * h_nPartners * h_child;
	
	double Rf=0.0;
	
	// If this male has partners
	// then also take into account their age
	// (the mean age of all partners, for simplicity):
	if (_individual[uid].get_nCurrSexPartner()>0)
	{
		double avg_age_f = meanAgePartners(uid);
		Rf = sexAct_reduce_age(avg_age_f);
		Rsex *= Rm*Rf;
	}
	
	// If males does NOT have partners => only sex workers
	// A reduction factor applies to his sex rate
	// (else, would visit CSW as if free, so too high rate)
	if (_individual[uid].get_nCurrSexPartner()==0){
		Rsex *= _sexAct_CostSexWork_reduction * Rm ;
	}
	
	// Draw the number of sex acts
	nSexActs = Rsex==0? 0 : poisson(Rsex*period);
	
//	if (save_trace_file) ff<<Rf<<","<<Rsex<<","<<nSexActs<<","<<age_m<<","<<nPartners<<endl;
	
	// update the number of sex acts
	// that will be performed by this male during this period
	_individual[uid].set_nSexActs_period(nSexActs);
	
	if(true /*save_trace_file*/){
		
		vector<double> newrow;
		vector<string> cname;
		
		newrow.push_back(uid); cname.push_back("uid");
		newrow.push_back(rg_m); cname.push_back("riskgrp");
		newrow.push_back(nSexActs); cname.push_back("nsex");

		_rec_sexact.push_back(newrow);
	}
	
}



double Population::sexAct_reduce_age(double age)
{
	/// Reduction sexual activity because of age
	
	double a_peak	= _sexAct_reduce_age_param[0];
	double sigma	= _sexAct_reduce_age_param[1];
	double q		= _sexAct_reduce_age_param[2];
	
	double res=0.0;
	if(age < a_peak)	res = pow((age-_ageSexMin)/(a_peak-_ageSexMin),2.0);
	if(age >= a_peak)	res = exp(-pow(fabs(age-a_peak)/sigma, q));
	return res;
}

double	Population::sexAct_reduce_riskGroup(int r)
{
	/// Reduce sexual activity base on individual's risk group
	
	double eps = _sexAct_reduce_risk_param[0];
	double rmax = _maxRiskGroup;
	
	return (1-eps)*r/rmax + eps;
}

double	Population::sexAct_reduce_STIsymptom(bool symptom, Gender g)
{
	/// Reduce sexual activity base on STI symptom manifestation
	
	// Reduction is gender dependent
	double eps = _sexAct_reduce_STIsymptom_param[0];
	if (g==female)
		eps= _sexAct_reduce_STIsymptom_param[1];
	
	return symptom?eps:1.0;
}


double	Population::sexAct_reduce_nPartners(int nPartners)
{
	/// Reduce sexual activity based on number of current partners
	
	double c = _sexAct_reduce_nPartner_param[0];
	return 2.0/(1.0+exp(-c*nPartners))-1.0;
}



void Population::sexAct_distribute_partnerTypes(unsigned long uid,
												bool save_trace_file)
{
	/// FOR A GIVEN NUMBER OF TOTAL SEX ACTS IN A PERIOD,
	/// DISTRIBUTE AMONG PARTNER TYPES (SPOUSE, CASUAL, CSW)
	
	// Trace file
//	ofstream tracefile(_DIR_OUT + "sexDistribPartn.out",ios::app);
	
	// Total number of sex acts
	int N = _individual[uid].get_nSexActs_period();
	
	// Total number of partners
	int nPartn = _individual[uid].get_nCurrSexPartner();
	
	// Males in partnerships
	if (N>0 && nPartn>0)
	{
		// Probability sex act with sex worker
		int rg = _individual[uid].get_riskGroup();
		double pw = probaSex_sexworker(rg);
		
		// Probability sex act with spouse
		double alpha = _sexAct_proba_distribute_partnerTypes_prefSpouse;
		
		int ns = _individual[uid].get_nCurrSpouse();
		int nc = _individual[uid].get_nCurrCasual();
		
		bool onenc = (nc>0)? 1 : 0;
		
		double ps = alpha*ns/(ns+nc)*onenc + (1-pw)*(1-onenc);
		
		// Probability sex act with casual partner
		double pc = 1-ps-pw;
		
		vector<double> proba(3);
		proba[0] = ps;
		proba[1] = pc;
		proba[2] = pw;
		
		// DEBUG:
		//cout << endl<<ns<<","<<nc<<" PROBA:"; displayVector(proba);
		
		// Now that we have the probabilities,
		// we can draw the multinomial random variable
		
		vector<unsigned int> NpartnerType = multinomial(N, proba);
		
		_individual[uid].set_nSexActs_spouse_period(NpartnerType[0]);
		_individual[uid].set_nSexActs_casual_period(NpartnerType[1]);
		_individual[uid].set_nSexActs_sexworker_period(NpartnerType[2]);
		
		// Traces
		if (save_trace_file){
//			tracefile << _simulationTime <<","<< uid <<",";
//			tracefile << _individual[uid].get_riskGroup() <<","<< _individual[uid].STI_anyInfection() <<",";
//			tracefile << NpartnerType[0] << ",";
//			tracefile << NpartnerType[1] << ",";
//			tracefile << NpartnerType[2] << ",";
//			tracefile << ns << "," << nc << endl;
		}
	}
	
	// Single males => only sex acts with sex workers
	if (N>0 && nPartn==0)
	{
		_individual[uid].set_nSexActs_sexworker_period(N);
		
		// Traces
		if (save_trace_file)
		{
//			tracefile << _simulationTime <<","<< uid <<",";
//			tracefile << _individual[uid].get_riskGroup() <<","<< _individual[uid].STI_anyInfection() <<",";
//			tracefile << -999 << ",";
//			tracefile << -999 << ",";
//			tracefile << N  << ",";
//			tracefile << 0 << "," << 0 << endl;
		}
	}
	
	// No sex acts: nothing to do (just assign value 0)
	if (N==0)
	{
		_individual[uid].set_nSexActs_spouse_period(0);
		_individual[uid].set_nSexActs_casual_period(0);
		_individual[uid].set_nSexActs_sexworker_period(0);
	}
}



void Population::sexAct_distrib_within_prtnrType(unsigned long uid,
												 string prtnrType,
												 int nPrtnr,
												 int nsexact,
												 bool saveTraceFile){
	/// Given the total number of sex acts with a partner type "prtnrType",
	/// distribute sex acts among individual partners
	
	stopif(prtnrType!="spouse" && prtnrType!="casual", "Unknown partner type.");
	
	if (nPrtnr>0 && nsexact>0)
	{
		// Define probabilities for individual partners.
		// Equal weight
		// ((think about including age, STI, HIV dependence,...))
		double ps = 1.0/nPrtnr;
		vector<double> proba(nPrtnr,ps);
		
		// Distribute sex acts among each individual partner
		vector<unsigned int> N_sexacts_indiv = multinomial(nsexact, proba);
		
		// Assign the drawn value
		for (int i=0; i<nPrtnr; i++){
			// Retrieve UID of ith parther:
			unsigned long uid_prtnr_i = 0;
			if(prtnrType=="spouse") uid_prtnr_i = _individual[uid].getSpouseUID()[i];
			if(prtnrType=="casual") uid_prtnr_i = _individual[uid].getCasualUID()[i];
			
			// Determine if female partner has any symptomatic STI.
			// If she does, then reduce the number of sex act with her.
			if(_individual[uid_prtnr_i].is_symptomatic()){
				// DEBUG
//				cout<< "symptomatic partner:";
//				displayVector(_individual[uid_prtnr_i].get_STIduration());
//				cout<< "STI[0]:"<<STInameString(_individual[uid_prtnr_i].get_STI()[0].get_name()) << endl;
//				cout<<"reduction:"<<_sexAct_reduce_STIsymptom_param[1]<<endl;
				// ------
				N_sexacts_indiv[i] = (int)(_sexAct_reduce_STIsymptom_param[1]*N_sexacts_indiv[i]);
			}
			
			add_sexAct(uid, uid_prtnr_i, N_sexacts_indiv[i],saveTraceFile);
		}
	}
}

void Population::sexAct_choose_CSW(unsigned long uid,
								   int nsexact,
								   bool saveTraceFile){
	/// Given the total number of sex acts with CSWs,
	/// assign sex acts to _ONE_ randomly chosen CSW.
	
	if (nsexact>0 && at_least_one_CSW()){
		
		// Choose the CSW:
		unsigned long uid_csw = chooseRandomCSW();
		
		// Assume all sex acts will be with THAT SAME CSW
		// during the given period.
		// (can be a problem if simulation time step is large)
		add_sexAct(uid, uid_csw, nsexact, saveTraceFile);
		
		// Although not a formal partner
		// intercourse with CSW counts in lifetime partner:
		_individual[uid].increment_nLifetimePartner();
		_individual[uid_csw].increment_nLifetimePartner();
		
		// Record that this individual has ever visited a CSW
		_individual[uid].set_ever_visited_CSW(true);
		
		// Set the duration of last visit to a tiny positive number
		// (enables incrementation)
		_individual[uid].set_lastVisitCSWDuration(0.001);
	}
}


void Population::sexAct_distribute_individualPartners(unsigned long uid,
													  bool save_trace_file)
{
	/// FOR A GIVEN NUMBER OF SEX ACTS BY PARTNER TYPES,
	/// DISTRIBUTE THOSE BETWEEN INDIVIDUAL PARTNERS
	
	stopif(_individual[uid].get_gender()!=male,
		   "Sex acts distribution implemented from male's point of view only");
	
	// Number of spouses and casual partners
	int ns = _individual[uid].get_nCurrSpouse();
	int nc = _individual[uid].get_nCurrCasual();
	
	// Number of sex acts with spouses
	int N_sexacts_s = _individual[uid].get_nSexActs_spouse_period();
	// Number of sex acts with casual partners
	int N_sexacts_c = _individual[uid].get_nSexActs_casual_period();
	// Number of sex acts with CSW
	int N_sexacts_csw = _individual[uid].get_nSexActs_sexworker_period();
	
	// DISTRIBUTE SEX ACTS AMONG INDIVIDUAL SPOUSES & CASUALS
	sexAct_distrib_within_prtnrType(uid, "spouse", ns, N_sexacts_s, save_trace_file);
	sexAct_distrib_within_prtnrType(uid, "casual", nc, N_sexacts_c, save_trace_file);
	// CHOOSE THE INDIVIDUAL CSW IF SEX ACTS WITH COMMERCIAL SEX
	sexAct_choose_CSW(uid, N_sexacts_csw, save_trace_file);
	
	// TO DO: implement check for limit on total sex acts for females
}


unsigned long Population::count_sexActs(int riskgroup){
	/// Counts the number of sex acts during a time step
	/// for a given risk group
	
	unsigned long count = 0 ;
	
	for (int uid=0; uid<_size; uid++) {
		if(_individual[uid].isAlive() &&
		   _individual[uid].get_riskGroup()==riskgroup){
			count += _individual[uid].get_nSexActs_period();
		}
	}
	return count;
}


double Population::probaSex_type0(int r1,int r2)
{
	/// Probability that a given sex act is of type 0 (condom)
	
	double rmax = (double)(_maxRiskGroup);
	double x1 =_sexAct_condom_param[0];
	double x2 =_sexAct_condom_param[1];
	
	return x1*exp(-(r1+r2)*x2/rmax);
}


double Population::probaSex_type1(int r1,int r2)
{
	/// Probability that a given sex act is of type 1 (standard risk)
	
	double probaType0 = probaSex_type0(r1,r2);
	return _sexAct_TypeLowHigh*(1-probaType0);
}



void Population::sexAct_distribute_ActTypes(unsigned long uid)
{
	/// UPDATES THE NUMBER OF SEX ACTS
	/// FOR A GIVEN INDIVIDUAL ("uid").
	/// THE TYPE OF SEX ACT IS DETERMINED (~MULTINOMIAL) HERE.
	
	stopif(_individual[uid].get_gender() != male, "must be male");
	stopif(!_individual[uid].isAlive(),"must be alive");
	
	// Retrieve the UIDs of partners with whom sex acts occurred
	vector<unsigned long> UIDpartners_sex = _individual[uid].get_UID_sexAct_period();
	int Np = UIDpartners_sex.size();
	
	// Reset for this period
	_individual[uid].reset_UID_n_sexAct_Type();
	
	// Current individual has at least one partner
	if (Np>0)
	{
		// Risk group of the male
		int rg = _individual[uid].get_riskGroup();
		
		for (int p=0; p<Np; p++)
		{
			// TO DO: segregate spouse vs casual
			
			// Risk group of the female partner
			int rg_p = _individual[UIDpartners_sex[p]].get_riskGroup();
			
			// Check if sex with CSW
			bool isCSW = (rg_p==_CSWriskGroup?true:false);
			
			// Retrieve the number of sex acts with
			// sex partner #p during given period
			unsigned int nSexActs_p = _individual[uid].get_UID_n_sexAct_period()[p];
			
			// Determine the type of sex acts (condom,low risk, high risk)
			// with partner #p. This is based on risk groups.
			if (nSexActs_p>0){
				// Calculate the proba for the multinomial random variable
				vector<double> proba(3);
				
				if (!isCSW)
				{
					proba[0] = probaSex_type0(rg, rg_p);	// sex with condom
					proba[1] = probaSex_type1(rg, rg_p);	// sex w/o condom, low risk
					proba[2] = 1-proba[0]-proba[1];			// sex w/o condom, high risk
				}
				
				if (isCSW)
				{
					// TO DO: Think about proper formula
					// initializa from a file
					proba[0] = 0.5*(_maxRiskGroup-rg)/_maxRiskGroup;	// sex with condom
					proba[1] = 0.5*(1-proba[0]);		// sex w/o condom, low risk
					proba[2] = 1-proba[0]-proba[1];		// sex w/o condom, high risk
				}
				
				vector<unsigned int> nSexType = multinomial(nSexActs_p, proba);
				
				// INCREMENTLY update the sex acts types for every partners
				_individual[uid].set_UID_n_sexAct_Type_period(nSexType);
			}
			
			if (nSexActs_p ==0){
				// Even if no sex act occurs with this partner,
				// _UID_n_sexAct_Type_period must be updated (with zeros)
				vector<unsigned int> zeros(3,0);
				_individual[uid].set_UID_n_sexAct_Type_period(zeros);
			}
		}
	}  //if (Np>0)
}



void Population::update_sexActs(double period, bool save_trace_file)
{
	/// For a given period, updates the number
	/// of sex acts:
	/// - Performed by males
	/// - Possibles for females
		
	for (int uid=0; uid<_size; uid++)
	{
		// Only for males
		if (_individual[uid].isAlive() && (_individual[uid].get_gender()==male) )
		{
			// Number of sex acts for a male
			sexActs_number_males(uid,period,save_trace_file); // force not write trace file because huge!
			
			// Distribute sex acts between partner types (spouse, casual, sex worker)
			sexAct_distribute_partnerTypes(uid,save_trace_file); // force not write trace file because huge!
			
			// If all partners are not discordant and no sex with CSW,
			// then skip sex acts distribution
			
			bool atLeastOneDiscord = STI_atLeastOneDiscordPartner(uid);
			
			if (atLeastOneDiscord || _individual[uid].get_nSexActs_sexworker_period()>0)
			{
				// TO DO: CODE CAN BE MORE OPTIMIZED
				// IF DISCORDANCE FURTHER CHECKED AT SPOUSAL AND CASUAL LEVEL
				
				// Distribute sex acts between all female partners (within types)
				sexAct_distribute_individualPartners(uid,save_trace_file);
				
				// Distribute sex act types for each individual
				sexAct_distribute_ActTypes(uid);
			}
		}
	}
}



void Population::add_sexAct(unsigned long uid1, unsigned long uid2,
							int nSexActs,
							bool save_trace_file)
{
	/// ADD SEX ACTS BETWEEN UID1 AND UID2
	
	// Must do it on both sides
	// to keep consistency
	_individual[uid1].add_sexAct(uid2, nSexActs);
	_individual[uid2].add_sexAct(uid1, nSexActs);
	
	_individual[uid1].increase_nSexActs_lifetime(nSexActs);
	_individual[uid2].increase_nSexActs_lifetime(nSexActs);
	
	// trace file
	if(save_trace_file)
	{
		ofstream tracefile(_DIR_OUT + "sexacts.out", ios::app);
		
		tracefile << _simulationTime << ",";
		tracefile << uid1 << ",";
		tracefile << _individual[uid1].get_gender() <<",";
		tracefile << _individual[uid1].get_riskGroup() << ",";
		tracefile << _individual[uid1].STI_anyInfection() << ",";
		tracefile << _individual[uid1].is_symptomatic() << ",";
		tracefile << _individual[uid1].get_age() << ",";
		tracefile << nSexActs <<endl;
		
		tracefile << _simulationTime << ",";
		tracefile << uid2 << ",";
		tracefile << _individual[uid2].get_gender() <<",";
		tracefile << _individual[uid2].get_riskGroup() << ",";
		tracefile << _individual[uid2].STI_anyInfection() << ",";
		tracefile << _individual[uid2].is_symptomatic() << ",";
		tracefile << _individual[uid2].get_age() << ",";
		tracefile << nSexActs <<endl;
	}
}




void Population::reset_sexActs()
{
	// Reset sex acts at each time step
	
	for (int i=0; i<_size; i++)
	{
		if (_individual[i].isAlive())
		{
			_individual[i].set_nSexActs_period(0);
			_individual[i].set_nSexActs_spouse_period(0);
			_individual[i].set_nSexActs_casual_period(0);
			_individual[i].set_nSexActs_sexworker_period(0);
			
			
			_individual[i].reset_UID_sexAct_period();
			_individual[i].reset_UID_n_sexAct_period();
			_individual[i].reset_UID_n_sexAct_Type();
		}
		// TO DO : check if there are other memnbers to reset
	}
}



vector<unsigned long> Population::pregnantPotentialFemales()
{
	/// Returns UID of all females that could
	/// potentially become pregnant
	/// (must have had sex with no condom)
	
	vector<unsigned long> res;
	
	for(unsigned long i=0; i<_size; i++)
	{
		if(getIndividual(i).isAlive() &&
		   getIndividual(i).get_gender()==female &&
		   //getIndividual(i).nSexActType1or2()>0 && // <-- TO DO: should be included but sex act type not recorded in females. Should retrieve all males partners and test nSexActType1or2>0 on males. For now keep it like that bc low  condom use, so very likely at least one sex act w/o condom
		   (!getIndividual(i).get_isPregnant()) &&
		   getIndividual(i).get_age()<_maxPregnantAge){
			res.push_back(getIndividual(i).get_UID());
		}
	}
	
	return res;
}



// ************ STI transmissions ************



void Population::STI_set_initial_prevalence(string filename)
{
	/// SET INITIAL STI PREVALENCES BY RISK GROUP
	///
	/// EXPECTED FILE FORMAT:
	/// FIRST COLUMN:	STI NAMES (e.g. HIV, HSV2, Ct, ...)
	/// SECOND COLUMN:	PREVALENCE FOR RISK GROUP "0"
	/// THIRD COLUMN:	PREVALENCE FOR RISK GROUP "1"
	/// FOURTH COLUMN:	PREVALENCE FOR RISK GROUP "2"
	/// LAST COLUMN:	PREVALENCE FOR CSW (HIGHEST RISK GROUP)
	

	force_seed_reset();
	
	vector<string> stiname;
	vectorFromCSVfile_string(stiname,filename.c_str(),1);
	
	int n_sti = stiname.size();
	
	vector<vector<double> > prevRiskGroup(_maxRiskGroup+2);
	
	// Store in a column vector all STI prevalence for a given risk group
	
	for (int j=0; j<_maxRiskGroup+2; j++){
		vectorFromCSVfile(prevRiskGroup[j], filename.c_str(), j+2);
		//DEBUG: cout << endl << "J="<<j;
		//displayVector(prevRiskGroup[j]);
	}
	
	for (int uid=0; uid<_size; uid++)
	{
		int rg = _individual[uid].get_riskGroup();
		if (rg==_CSWriskGroup) rg = _maxRiskGroup+1;
		
		for (int s=0; s<n_sti; s++)
		{
			string currSTI = stiname[s];
			double prev = prevRiskGroup[rg][s]; // Wanted prevalence for this STI and this risk group
			
			// Position of this STI in the (fixed) list of _STI
			int ss = STI_find_index_position(StringToSTIname(currSTI));
			
			bool isInfected = binom(prev, 1);
			
			if (isInfected){
				// Define (randomly) duration of infection
				double maxDuration = -999;
				
				if (currSTI=="HIV" || currSTI=="HSV2" || currSTI=="HPV")
					maxDuration = min(7.0 ,_individual[uid].get_age()-_ageSexMin);
				
				if (currSTI=="Ct" || currSTI=="Ng" || currSTI=="Hd"
					|| currSTI=="Bv" || currSTI=="Tv")
					maxDuration = 0.5;
				
				if(currSTI=="Tp") maxDuration = 1.5;
				
				double duration = uniform01()*maxDuration;
				_individual[uid].set_STIduration(StringToSTIname(currSTI), duration);
				
				// Symptoms
				Gender g = _individual[uid].get_gender();
				double probaSymptom=-99;
				if (g==female)	probaSymptom = _STI[ss].get_proba_symptomatic_female();
				if (g==male)	probaSymptom = _STI[ss].get_proba_symptomatic_male();
				bool isSymptomatic = false;
				double u = uniform01();
				if (u<probaSymptom) isSymptomatic = true;
				_individual[uid].set_STIsymptom(ss,isSymptomatic);
				
				//DEBUG
				//cout <<uid<< " -> RG="<<rg<<" ; STI="<<currSTI<<" ; prev="<<prev;
				//cout <<" ; DURATION="<<duration<<" ;SYMPT="<<isSymptomatic<< endl;
			}
		}
	} // end for 'uid'
}


void Population::STI_acquire(STIname stiname, unsigned long uid, double prd)
{
	/// Individual "uid" acquires STI "stiname"
	/// and STI duration for this infection is set at "prd"
	
	// DEBUG:
	//cout << endl << STInameString(stiname) << " acquired by "<<uid<<endl;
	_individual[uid].STI_acquireInfection(stiname, prd);
}


void Population::STI_transmission_indiv(STIname stiname,
										unsigned long uid1,
										unsigned long uid2,
										double timeStep,
										bool save_trace_file)
{
	/// STI transmission at the individual level.
	/// Determine who is infectious among uid1 and uid2.
	
	// Integrity checks:
	bool sti1 = _individual[uid1].get_STIduration(stiname)>0;
	bool sti2 = _individual[uid2].get_STIduration(stiname)>0;
	string errmsg = to_string(uid1)+" and " +to_string(uid2) + " are both not infectious with " +  STInameString(stiname);
	stopif((!sti1 && !sti2), errmsg);
	string errmsg2 = to_string(uid1)+" and " +to_string(uid2) + " are both infected with " +  STInameString(stiname);
	stopif((sti1 && sti2), errmsg2);
	
	unsigned long uid_inf	= (sti1?uid1:uid2);
	unsigned long uid_s		= (sti1?uid2:uid1);
	
	// "uid_s" acquires the STI:
	STI_acquire(stiname, uid_s, timeStep);
	// update list of secondary cases for "uid_inf":
	_individual[uid_inf].add_STI_secondary_cases(stiname, uid_s);
	
    // transmission of genetic material:
    _individual[uid_s].set_genome(_individual[uid_inf].get_genome(stiname), stiname);
    
    
	if (save_trace_file) {
		string fname = _DIR_OUT+"transmissions.out";
		ofstream f(fname.c_str(),ios::app);
		write_headers_if_emptyFile(fname,
								   "time,stiname,uid_from,uid_to,riskGroup_from,riskGroup_to,stiduration_from,gender_from");
		f << _simulationTime <<",";
		f <<  STInameString(stiname) <<",";
		f << uid_inf <<",";
		f << uid_s <<",";
		f << _individual[uid_inf].get_riskGroup() <<",";
		f << _individual[uid_s].get_riskGroup() <<",";
		f << _individual[uid_inf].get_STIduration(stiname) <<",";
		f << _individual[uid_inf].get_gender() << endl;
	}
}


double Population::STI_coInfection_oddsRatio(unsigned long uid, int sti)
{
	/// CALCULATES THE OVERALL ODDS RATIO
	/// FOR A INDIVIDUAL (uid) TO A SPECIFIC STI (sti)
	/// TAKING INTO ACCOUNT OTHER STI CO-INFECTION
	
	double themax = -99.9;
	
	bool anyOtherCoinfection = false;
	
	for (int s=0; s<_nSTImodelled; s++)
	{
		// Only loop through STIs infecting this individual
		if (_individual[uid].get_STIduration()[s]>0)
		{
			anyOtherCoinfection = true;
			
			// Odds-ratio for co-infection pair (sti,s)
			double tmp = _STI_SFincrease(sti,s);
			
			// if co-infections with more than 1 STI
			// then odds-ratio is the maximum of all odds-ratios
			themax = tmp>themax? tmp : themax;
		}
	}
	// If this individual has no other STI
	// then the odds-ratio is neutral (=1.0)
	double res = anyOtherCoinfection?themax:1.00;
	return res;
}


vector<double> Population::STI_CalcProbaTransmission(unsigned long uid_infect,
													 unsigned long uid_suscep)
{
	/// CALCULATE THE PROBABILITY OF TRANSMISSION
	/// FOR EVERY STI MODELLED
	/// BETWEEN TWO INDIVIDUALS
	
	// UIDs of all male's sex partners (including CSWs) during this period:
	unsigned long uid_male		= getUIDmale(uid_infect, uid_suscep);
	unsigned long uid_female	= getUIDfemale(uid_infect, uid_suscep);
	vector<unsigned long> UID_sexpartners_of_male = _individual[uid_male].get_UID_sexAct_period();
	
	
	// =====================================
	//    Multiple Proba of Transmission
	// =====================================
	
	vector<double> MPT(_nSTImodelled, 0.0);
	
	Individual I_inf = _individual[uid_infect];
	Individual I_sus = _individual[uid_suscep];
	
	for (int sti=0; sti<_nSTImodelled; sti++)
	{
		// This individual is infected
		// and the partner is suceptible
		// Also check the susceptible partner
		// is not vaccinated against this STI
		
		double dummy_debug = I_inf.get_STIduration()[sti];
		
		if (I_inf.get_STIduration()[sti]>0 &&
			I_sus.get_STIduration()[sti]==0 &&
			I_sus.get_STI_immunity()[sti]<1 ){
			
			// Infectivity for this individual and this sti
			double infectivity = I_inf.STI_IC()[sti];
			
			// Susceptibility factor of the partner for this sti
			double susceptFactor = I_sus.get_STIsusceptFactor()[sti];
						
			// Max probability of transmission per sex act
			double maxProba = I_inf.get_STI()[sti].get_probaMaxSexTransm();
			
			// Reduction factor associated with sex act type
			vector<double> SAT = I_inf.get_STI()[sti].get_probaSexTransm_SAT();
			
			// Number of sex acts of a given type
			vector<int> nSexType(3);
			
			// Sex act type is recorded only with males.
			// Find position of female in the list of
			// all male's sex partners:
			int p_position = 0;
			while (UID_sexpartners_of_male[p_position] != uid_female) p_position++;
			stopif(p_position > UID_sexpartners_of_male.size(),
				   "Can't find UID of sex partner when assigning sex act types!");
			
			// Assign the associated number of sex acts for this sex partner, by sex act type
			nSexType[0] = _individual[uid_male].get_UID_n_sexAct_Type0_period(p_position);
			nSexType[1] = _individual[uid_male].get_UID_n_sexAct_Type1_period(p_position);
			nSexType[2] = _individual[uid_male].get_UID_n_sexAct_Type2_period(p_position);
			
			// Probability of STI transmission
			// given multiple sex acts and STI co-infections:
			double MPT_sti = 1.0;
			
			// Loop through al sex acts:
			for (int T=0; T<3; T++)
			{
				if (nSexType[T]>0) // Calculate proba only if sex act occured
				{
					// Probability without any other STI co-infection
					double proba_T = infectivity * susceptFactor * SAT[T] * maxProba;
					
					// = = Co infections = =
					
					// odds-ratio given any other co-infection
					double oddsratio = STI_coInfection_oddsRatio(uid_suscep,sti);
					
					// proba taking into account the odds-ratio
					double proba_T_coinf = oddsratio * proba_T / (1.0 + (oddsratio-1.0) * proba_T);
					
					// product
					double tmp = pow(1.0 - proba_T_coinf, nSexType[T]);
					MPT_sti = MPT_sti*tmp;
				}
			}
			MPT_sti = 1.0 - MPT_sti;  // see formula of MPT in documentation
			MPT[sti] = MPT_sti;
		} // end if
	} // end for(sti)
	return MPT;
}



vector<unsigned long>  Population::STI_transmissions(double timeStep,
													 bool save_trace_file)
{
	/// - loop thru all sex acts
	/// - if any one of the 2 partners is STI+, calculate proba transmission
	/// - draw random variable and initiate transmission or not
	/// (no need to explicitly manage STI acquisition, already dealt with transmission)
	
	/// Returns a vector of incidence (incidence value for each STI)
	
	// Vector storing incidence values for each STI
	vector<unsigned long> incidence(_nSTImodelled,0);
	
	for (int uid=0; uid<_size; uid++)
	{
		// Check if this individual had sex acts
		
		// WARNING: This individual must be of one gender only (here chosen male)
		// because if gender is not checked, a partnership will be assessed TWICE
		// for transmission (once male --> female, then a 2nd time female --> male).
		// Hence, the number of sex acts is implicitly doubled
		// (actually only when first transmission fails)
		
		// WARNING 2: This (selecting one gender) would not work
		// if homosexual partnerships modelled.
		
		if ( _individual[uid].get_nSexActs_period() > 0 &&
			_individual[uid].isAlive() &&
			_individual[uid].get_gender()==male){
			
			// Retrieve UIDs of sex partners during this period
			vector<unsigned long> UIDsexpartners = _individual[uid].get_UID_sexAct_period();
			
			// Number of sex partners
			int Np = UIDsexpartners.size();
			
			// Loop through all sex partners
			for (int p=0; p<Np; p++)
			{
				// UID of the partner #p
				unsigned long uid_p = UIDsexpartners[p];
				
				// == STI transmission TO partners ==
				
				// Check first if this partnership
				// is indeed discordant
				// (else, nothing to do!)
				
				if (STI_isDiscordPartner(uid, uid_p)){
					// We don't know who's infected with which STI,
					// hence look at both partners and take the 'max'
					// of transmission probability.
					// (not really max, bc value is either 0 or strictly positive proba)
					
					vector<double> probaTransm_to	= STI_CalcProbaTransmission(uid, uid_p);
					vector<double> probaTransm_from = STI_CalcProbaTransmission(uid_p, uid);
					vector<double> probaTransm		= max_vector(probaTransm_to, probaTransm_from);
					
					// DEBUG ===
					if (_individual[uid].get_STIduration(HIV)>0 &&
						_individual[uid_p].get_STIduration(Tp)>0) {
//						cout << endl <<_simulationTime<<  " Attempt to transmit: "<<uid<<"->"<<uid_p<<endl;
//						
//						cout<<"HIV durations:"<<_individual[uid].get_STIduration(HIV)<<";"<<_individual[uid_p].get_STIduration(HIV)<<endl;
//						cout<<"Tp durations:"<<_individual[uid].get_STIduration(Tp)<<";"<<_individual[uid_p].get_STIduration(Tp)<<endl;
//						displayVector(probaTransm_to);
//						displayVector(probaTransm_from);
//						displayVector(probaTransm);
//						
//						vector<double> dummy	= STI_CalcProbaTransmission(uid, uid_p);
						
					}
					// =========
					 
					
					
					for (int sti=0; sti<_nSTImodelled; sti++)
					{
						// Loop only through STIs that
						// have a chance to be transmissible
						if (probaTransm[sti]>0){
							// Actual transmission if random variable
							// translates into 'success' (of transmission)
							
							double u = uniform01();
							bool successTransmission = (u<probaTransm[sti]);
							
							if (successTransmission){
								STI_transmission_indiv(_STI[sti].get_name(),
													   uid, uid_p,
													   timeStep,
													   save_trace_file);
								incidence[sti]++;
								// DEBUG
								//cout << "** TRANSMISSION of "<<STInameString(sti);
								//cout << " from "<<uid<< " to " <<uid_p<<endl;
							}
						}
					} // end_loop on all STIs
				} // end_if has any infection
				
				// ****** IMPORTANT ******
				// STI acquisition FROM partners is dealt with
				// when the loop on individual goes through
				// the infected partner(s)
				// ***********************
			}
		} // end_if had any sex act during this period
	}// end_loop on all individuals
	
	return incidence;
}



double Population::STI_proba_MTCT(STIname sti, double stiduration)
{
	/// Calculate the actual probability of MTCT
	double proba = _STI[positionSTIinVector(sti, _STI)].get_proba_MTCT();
	// special case for syphilis
	// where mtct is believed to
	// decrease with time.
	if (sti==Tp) proba = proba/(1.0+exp(3.0*(stiduration-2.0)));
	
	// DEBUG
	//cout << STInameString(sti)<<"_"<< stiduration << " DEBUG MTCT PROBA: "<<proba<<endl;
	// =====
	
	return proba;
}



void Population::STI_update_naturalClearance()
{
	/// ALL STIS THAT ARE NATURALLY CLEARED
	/// BY IMMUNE SYSTEM ARE TERMINATED HERE
	/// AFTER A PRE-SPECIFIED PERIOD
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() && _individual[uid].STI_anyInfection())
		{
			for (int s=0; s<_nSTImodelled; s++) // loop thru all STIs of this individual
			{
				if (_individual[uid].get_STIduration()[s]>0)
				{
					// Retrieve the clearence duration of STI #s
					double maxSTIduration = _individual[uid].get_STI()[s].get_naturalClearanceDuration();
					
					if (_individual[uid].get_STIduration()[s] > maxSTIduration)
					{
						_individual[uid].set_STIduration(s,0.0);    // STI is cleared (_STIduration=0)
						_individual[uid].set_STIsymptom(s, false);  // STI symptoms disappear
					}
				}
			}
		}
	}
}


void Population::STI_update_templateToAllIndividuals()
{
	/// UPDATE STI TEMPLATE TO ALL INDIVIDUAL
	
	/// THIS FUNCTION MUST BE CALLED WHEN THE STI PARAMETERS
	/// ARE CHANGED AT THE POPULATION LEVEL (E.G. DURING CALIBRATION)
	/// AND HAVE TO BE APPLIED TO ALL INDIVIDUALS
	
	for (int uid=0; uid<_size; uid++)
		if (_individual[uid].isAlive())	_individual[uid].set_STI(_STI);
}




double Population::STI_prevalence(STIname stiname)
{
	/// GLOBAL PREVALENCE OF ONE STI
	
	unsigned long N=0;	// total number of individuals alive
	unsigned long infected=0; // total number of infected individuals with this STI
	
	unsigned int sti_i = positionSTIinVector(stiname, _STI);
	
	for (int uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive()){
			N++;
			if(_individual[uid].get_STIduration()[sti_i]>0)   // slow code:(_individual[uid].get_STIduration(stiname)>0)
				infected++;
		}
	}
	return (double)(infected)/N;
}



vector<double>	Population::STI_prevalences(vector<STIname> stinames)
{
	/// GLOBAL PREVALENCE OF SELECTED STIs
	
	vector<double> prev;
	
	for (int s=0; s<stinames.size(); s++)
		prev.push_back(STI_prevalence(stinames[s]));
	
	return prev;
}




double Population::STI_prevalence(STIname s, int riskGroup)
{
	/// GLOBAL PREVALENCE OF ONE STI
	/// WITHIN A GIVEN RISK GROUP
	
	unsigned long N=0;			// total number of individuals alive, given risk group
	unsigned long infected=0;	// total number of infected individuals with this STI, given risk group
	
	
	for (int uid=0; uid<_size; uid++){
		if (_individual[uid].isAlive() &&
			_individual[uid].get_riskGroup()==riskGroup){
			N++;
			if (_individual[uid].get_STIduration(s)>0)	infected++;
		}
	}
	return (double)(infected)/N;
}


dcMatrix Population::STI_prevalence_by_riskGroup(vector<STIname> stinames)
{
	/// MATRIX WHERE:
	/// ROWS = STIS
	/// COLUMNS = RISK GROUPS, THE LAST ONE BEING CSW
	
	
	dcMatrix res(0,0);
	
	for (int s=0; s<stinames.size(); s++)
	{
		// Row for the matrix
		vector<double> row_sti;
		
		for (int rg=0; rg<=_maxRiskGroup; rg++)
		{
			row_sti.push_back(STI_prevalence(_STI[s].get_name(),rg));
		}
		// last column for CSW
		row_sti.push_back(STI_prevalence(_STI[s].get_name(),_CSWriskGroup));
		
		res.addRowVector(row_sti);
		row_sti.clear();
	}
	
	return res;
}




vector<double> Population::STI_prevalence_by_age(STIname s,
												 vector<double> agebreaks)
{
	/// GLOBAL PREVALENCE BY AGE OF ONE STI
	
	vector<double> x_nosti;
	vector<double> x_sti;
	
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() )
		{
			if (_individual[uid].get_STIduration(s)>0){
				x_sti.push_back(_individual[uid].get_age());
			}
			else{
				x_nosti.push_back(_individual[uid].get_age());
			}
		}
	}
	
	// Distribution of counts
	vector<unsigned long> d_nosti = distribution(x_nosti, agebreaks);
	vector<unsigned long> d_sti = distribution(x_sti, agebreaks);
	
	// Prevalenve with an age bucket
	vector<double> res;
	for (int i=0; i<d_nosti.size(); i++)
	{
		double tmp = 0;
		
		if (d_sti[i]+d_nosti[i]>0)
			tmp =(double)d_sti[i]/(double)(d_sti[i]+d_nosti[i]);
		
		res.push_back(tmp);
	}
	return res;
}


vector<double> Population::STI_prevalence_by_age(STIname s,
												 Gender g,
												 vector<double> agebreaks)
{
	/// GLOBAL PREVALENCE BY AGE
	/// FOR A GIVEN GENDER
	/// OF ONE STI
	
	
	vector<double> x_nosti;
	vector<double> x_sti;
	
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() &&
			_individual[uid].get_gender()==g)
		{
			if (_individual[uid].get_STIduration(s)>0)	{
				x_sti.push_back(_individual[uid].get_age());
			}
			else{
				x_nosti.push_back(_individual[uid].get_age());
			}
		}
	}
	
	// Distribution of counts
	vector<unsigned long> d_nosti = distribution(x_nosti, agebreaks);
	vector<unsigned long> d_sti = distribution(x_sti, agebreaks);
	
	// Prevalence within an age bucket
	vector<double> res;
	
	for (int i=0; i<d_nosti.size(); i++)
	{
		double tmp = 0;
		
		if (d_sti[i]+d_nosti[i]>0)
			tmp =(double)d_sti[i]/(double)(d_sti[i]+d_nosti[i]);
		
		res.push_back(tmp);
	}
	return res;
}




vector<double> Population::STI_prevalence_by_age(STIname s,
												 int riskGroup,
												 vector<double> agebreaks)
{
	
	/// GLOBAL PREVALENCE BY AGE
	/// FOR A GIVEN GENDER
	/// OF ONE STI
	
	
	vector<double> x_nosti;
	vector<double> x_sti;
	
	
	for (int uid=0; uid<_size; uid++)
	{
		if (_individual[uid].isAlive() &&
			_individual[uid].get_riskGroup()==riskGroup)
		{
			//			if (_individual[uid].get_STIduration()[s]>0)
			if (_individual[uid].get_STIduration(s)>0)
			{
				x_sti.push_back(_individual[uid].get_age());
			}
			else
			{
				x_nosti.push_back(_individual[uid].get_age());
			}
		}
		
	}
	
	// Distribution of counts
	
	vector<unsigned long> d_nosti = distribution(x_nosti, agebreaks);
	vector<unsigned long> d_sti = distribution(x_sti, agebreaks);
	
	// Prevalenve with an age bucket
	
	vector<double> res;
	
	for (int i=0; i<d_nosti.size(); i++)
	{
		double tmp = 0;
		
		if (d_sti[i]+d_nosti[i]>0)
			tmp =(double)d_sti[i]/(double)(d_sti[i]+d_nosti[i]);
		
		res.push_back(tmp);
	}
	
	return res;
}


void Population::update_STI_mtct_cumcount(vector<bool> mtct){
	/// Increase the cumulative count of
	/// MTCT events
	
	for(int i=0; i<mtct.size();i++)
		if(mtct[i]) _STI_mtct_cumcount[i]++;
}



vector<unsigned long> Population::Reff_inst_all_infectious(STIname stiname){
	
	/// "Instantaneous" Effective reproductive number for all infectious individuals for a given STI
	/// (only for the _timeStep period)
	
	vector<double> Reff_dist;
	vector<unsigned long> cnt;
	
	for(unsigned long uid=0; uid<_size; uid++){
		if(_individual[uid].isAlive() &&
		   _individual[uid].get_STIduration(stiname)>0)
		{
			unsigned long n = _individual[uid].get_STI_secondary_cases(stiname).size();
			cnt.push_back(n);
		}
	}
	return cnt;
}


double Population::Reff_inst_mean(STIname stiname){
	
	/// "Instantaneous" Mean effective reproductive number for a STI
	/// (only for the _timeStep period)
	
	vector<unsigned long> cnt = Reff_inst_all_infectious(stiname);
	double Reff = 0.0;
	
	if(cnt.size()>0) Reff = meanElements(cnt);
	return Reff;
}


void Population::update_secondary_cases(STIname stiname){
	
	/// "Cumulative" Effective reproductive number for all infectious individuals for a given STI.
	/// Assume simulation is run (else this function makes no sense)
	
	// Retrieve secondary cases for the current time step:
	vector<unsigned long> cnt = Reff_inst_all_infectious(stiname);
	// Determine STI index:
	int i_sti = STI_find_index_position(stiname);
	// Append those new cases to existing ones:
	for(int i=0; i<cnt.size(); i++) _secondary_cases[i_sti].push_back(cnt[i]);
}

double Population::Reff_cum_mean(STIname stiname){
	/// Calculate Reff cumulative (since simulation started)
	/// for a given STI
	
	int i_sti = STI_find_index_position(stiname);
	return(meanElements(_secondary_cases[i_sti]));
}




int Population::STI_find_index_position(STIname s)
{
	/// FOR A GIVEN STI NAME
	/// RETURNS ITS POSITION
	/// (FROM THE VECTOR '_STI')
	
	int	res = -99;
	bool found = false;
	
	for (int i=0;i<_nSTImodelled;i++){
		if (_STI[i].get_name()==s){
			res=i;
			found=true;
		}
	}
	string errmsg = "STI name not found: " + STInameString(s);
	stopif(!found, errmsg);
	return res;
}



void Population::STI_save_all_infectivityCurves(string filerootname)
{
	/// SAVE INFECTIVITY CURVES
	/// FOR EVERY STI MODELLED
	
	for (int i=0; i<_STI.size(); i++){
		string name_i = filerootname + "_" + STInameString(_STI[i].get_name()) + ".out";
		_STI[i].check_infectivityCurve(name_i);
	}
}

vector< vector<double> > Population::get_infectivityCurve(STIname stiname, Gender g){
	
	/// Get infectivity curve for a given STI.
	///
	/// first vector returned: times (at which IC is calculated)
	/// second vector: infectivity curve values

	
	vector< vector<double> > x(2);
	
	// Time unit is years
	double tmax = 80;
	double dt = 2.0/365.0;
	
	int sti_i = positionSTIinVector(stiname, _STI);
	
	for(double t=0.0; t<=tmax; t+=dt){
		x[0].push_back(t);
		x[1].push_back(_STI[sti_i].infectivityCurve(t, g));
	}
	return x;
}



bool  Population::STI_isDiscordPartner(unsigned long uid1,unsigned long uid2)
{
	/// Check if this pair of indivudals is discordant
	/// for any STI modelelled
	
	// Integrity checks:
	stopif(_individual[uid1].get_gender()==_individual[uid2].get_gender(),"Individuals must be of different gender!");
	
	bool isDiscordant = false;
	bool sti1 =_individual[uid1].STI_anyInfection();
	bool sti2 =_individual[uid2].STI_anyInfection();
	
	// Concordant negative
	if (!sti1 && !sti2) isDiscordant = false;
	
	// If one partner has no STI but the other has at least one then discordant:
	if ( (sti1 && !sti2) || (!sti1 && sti2) ) isDiscordant = true;
	
	// Both partners have at least an STI,
	// need to check which ones
	if ( sti1 && sti2){
		vector<string> stilist1 = _individual[uid1].STI_listInfection();
		vector<string> stilist2 = _individual[uid2].STI_listInfection();
		
		for (int i=0; i<stilist1.size() && !isDiscordant; i++){
			if (!isElementPresent(stilist2, stilist1[i])) isDiscordant = true;
		}
	}
	return isDiscordant;
}



bool Population::STI_atLeastOneDiscordPartner(unsigned long uid)
{
	/// IS THERE AT LEAST ONE PARTNERSHIP WHICH IS DISCORDANT?
	
	bool res = false;
	
	// Retrieve all partners
	vector <unsigned long> uidp = _individual[uid].getPartnerUID();
	
	for (int i=0; i<uidp.size() && !res ; i++)
		res = STI_isDiscordPartner(uid, uidp[i]);
	
	return res;
}


bool Population::at_least_one_CSW(){
	
	/// IS THERE AT LEAST ONE CSW IN THE POPULATION?
	
	bool res = false;
	for(int i=0; i<_size; i++)
		if(_individual[i].isAlive() && _individual[i].get_riskGroup()==_CSWriskGroup){
			res = true;
			break;
		}
	return res;
}


/* ************************************* */
/* ************** HELPERS ************** */
/* ************************************* */


void Population::displayInfo(bool IndividualDetails)
{
	// May not be used, but initialized the correct count
	// of number of partnerships
	
	dcMatrix P	= getPartnershipsUID();
	dcMatrix SP	= getSpousesUID();
	
	
	// Display info
	
	cout << endl <<endl;
	cout << "================================"<<endl;
	cout << "====== Info on Population ======"<<endl;
	cout << "================================"<<endl <<endl;
	
	cout << "Size: "						<< _size<<endl;
	
	unsigned long nfem =census_Females();
	unsigned long nalive =census_alive();
	
	cout << "Alive: "						<< nalive<<endl;
	cout << "Number females: "				<< nfem << " (" << 100*nfem/nalive<< "%)"<< endl;
	cout << "Number females available: "	<< getUID_available(female).size() <<endl;
	cout << "Number males available: "		<< getUID_available(male).size() <<endl;
	cout << "Number of widowed: "			<< census_widow() <<endl;
	cout << "Number of males circumcised: "		<< census_circum() <<" ("<< 100*census_circum()/nalive << "%)" << endl;
	cout << "Target circumcision proportion: "		<< 100*_proportion_circum << "%" <<endl;
	
	coutline(40);
	
	cout << "Demographics parameters: "<<endl;
	
	cout << "Birth rate = "<< _birthRate <<endl;
	cout << "Infant Mortality rate = "<< _infantMortality <<endl;
	cout << "Child Mortality rate = "<< _childMortality <<endl;
	cout << "_deathParam_Weibull_scale = "<< _deathParam_Weibull_scale <<endl;
	cout << "_deathParam_Weibull_shape = "<< _deathParam_Weibull_shape <<endl;
	cout << "_deathParam_Weibull_scale_hiv = "<< _deathParam_Weibull_scale_hiv <<endl;
	cout << "_deathParam_Weibull_shape_hiv = "<< _deathParam_Weibull_shape_hiv <<endl;
	
	coutline(40);
	
	cout << "Partnerships formation parameters: "<<endl;
	
	//	cout << "mean age female = "<<_formation_meanAge_female;
	//	cout << "\t  variance age female = "<<_formation_varAge_female<<endl;
	
	cout << "age start full = "<< _formation_age_fullstart<<endl;
	cout << "age pivot = "<< _formation_age_pivot<<endl;
	cout << "shape age formation = "<<_formation_age_shape<<endl;
	cout << "formation min proba component = " << _formation_age_fmin <<endl;
	
	cout << "mean age gap = "<<_formation_agegap_mean<<endl;
	//cout << "variance age gap = "<<_formation_agegap_var<<endl;
	cout << " min proba component age gap = "<< _formation_agegap_fmin << endl;
	
	
	cout << "Max rate spousal formation: " << _spousalProgress_maxRate << endl ;
	cout << "mean age female new spouse: " << _spousalProgress_meanAge_f ;
	cout << "\t variance: "<< _spousalProgress_varAge_f << endl;
	cout << "mean age gap new spouse: " << _spousalProgress_meanGap ;
	cout << "\t variance: "<< _spousalProgress_varGap << endl;
	cout << "duration shape new spouse k1: " << _spousalProgress_durationK1 ;
	cout << "\t k2: "<< _spousalProgress_durationK2 << endl;
	cout << "Mean diff age gap multi-spousal: " << _spousalProgress_meanDiffAgeGap;
	cout << "\t variance: "<< _spousalProgress_varDiffAgeGap << endl;
	
	coutline(40);
	cout << "Number of risk groups: "<<_maxRiskGroup+1<<endl;
	coutline(40);
	
	cout << "Number of partnerships: " << _totalNumberPartnerships << endl;
	cout << "Number of spousal partnerships: " << _totalNumberSpousalPartnerships << endl;
	
	vector<unsigned long> pp(6);
	for (int p=0; p<pp.size(); p++)
	{
		pp[p] = census_Partnered(p);
		cout << "Indiv with "<< p <<" partner"; if (p>1) cout<<"s";
		cout << ": "<<pp[p]<<" ("<< (double)(pp[p])/_size<<")"<<endl;
	}
	
	
	coutline(40);
	cout << "Condom use parameters: "; displayVector(_sexAct_condom_param);
	cout << "Ratio Low v.s high risk sex type: "<<_sexAct_TypeLowHigh<<endl;
	cout << "Number of CSW : "<< census_CSW().size() << endl;
	
	coutline(40);
	cout << endl << "STI infections: " <<endl;
	for (int sti=0; sti<_nSTImodelled; sti++)
	{
		cout << STInameString(_STI[sti].get_name())<<" : " << census_STIinfected()[sti] << endl;
	}
	
	bool detailSTI = false;
	
	if (detailSTI)
	{
		for (int i=0; i<_STI.size(); i++) {
			_STI[i].displayInfo();
			
		}
	}
	
	coutline(40);
	
	// Increase in susceptibility when STI co-infections
	displayInfo_STI_SFincrease();
	
	cout << endl << " Rebound of HIV infectivity:"<<endl;
	cout << "(fraction of peak (early acute phase) HIV infectivity)"<<endl;
	
	for (int i=0;i<_nSTImodelled;i++)
	{
		if (STInameString(_STI[i].get_name())!="HIV")
		{
			cout<<STInameString(_STI[i].get_name())<<" --> "<<_RebHIV[i] <<endl;
		}
		
	}
	coutline(40);
	
	
	// Circumcision and reduced susceptibility
	
	cout << endl << " Circumcision and reduced susceptibility:"<<endl;
	cout << "(susceptibility reduction factor when circumcised)"<<endl;
	for (int i=0;i<_nSTImodelled;i++)
	{
		cout<<STInameString(_STI[i].get_name())
		<<" --> "<< _STI[i].get_circum_SF_reduction()<<endl;
	}
	
	
	coutline(40);
	
	
	
	// ==== --- If Individuals details asked --- ====
	
	if (IndividualDetails)
	{
		cout <<endl<< "UID\t alv\t Gdr\t Age\t RskG\t nSx\t nSp\t nMx\t ";
		cout << "nLf\t nLfS\t Div\t Wid\t";
		cout << "P1\t P2\t P3\t P4\t";
		cout << "HIV";
		cout << endl;
		
		for (int p=0;p< 120; p++) cout<<"=";
		cout<<endl;
		
		for (int p=0; p<_size; p++)
		{
			cout <<_individual[p].get_UID() << "\t"
			<< _individual[p].isAlive()<< "\t"
			<< _individual[p].get_gender()<< "\t"
			<< _individual[p].get_age()<< "\t"
			<< _individual[p].get_riskGroup()<< "\t"
			<< _individual[p].get_nCurrSexPartner()<< "\t"
			<< _individual[p].get_nCurrSpouse()<< "\t"
			<< _individual[p].get_nMaxCurrSexPartner()<< "\t"
			<< _individual[p].get_nLifetimePartner()<< "\t"
			<< _individual[p].get_nLifetimeSpouse()<< "\t"
			<< _individual[p].get_isDivorced()<< "\t"
			<< _individual[p].get_isWidow()<< "\t";
			
			// partners' UID (if any)
			if (_individual[p].get_nCurrSexPartner()>0) cout << _individual[p].getPartnerUID(0)<< "\t"; else cout<<"\t";
			if (_individual[p].get_nCurrSexPartner()>1) cout << _individual[p].getPartnerUID(1)<< "\t"; else cout<<"\t";
			if (_individual[p].get_nCurrSexPartner()>2) cout << _individual[p].getPartnerUID(2)<< "\t"; else cout<<"\t";
			if (_individual[p].get_nCurrSexPartner()>3) cout << _individual[p].getPartnerUID(3)<< "\t"; else cout<<"\t";
			
			
			cout << _individual[p].get_STIduration(HIV);
			cout << endl;
		}
	}
	cout<<endl;
}


void Population::displayInfo_SexActs_lastPeriod()
{
	int nPartnerMax = 4; // maximum partners displayed
	
	cout << "UID\t Gnd\t Rsk\t" ;
	for (int i=0; i<nPartnerMax; i++) {
		cout << "P"<<i<<"\t";
	}
	
	cout << "Nsex\t NsexSp\t NsexC\t NsexW";
	
	cout << endl;
	
	for (int i=0; i<_size; i++)
	{
		cout << i << "\t";
		cout << _individual[i].get_gender()<<"\t";
		cout << _individual[i].get_riskGroup()<<"\t";
		
		for (int p=0; p<nPartnerMax; p++)
		{
			if (_individual[i].get_nCurrSexPartner()>p)
			{
				cout << _individual[i].getPartnerUID(p);
				if (_individual[i].isSpouse(_individual[i].getPartnerUID(p))) cout<<"s";
				cout << "\t";
			}
			else cout<<"\t";
		}
		
		cout << _individual[i].get_nSexActs_period()<<"\t";
		cout << _individual[i].get_nSexActs_spouse_period()<<"\t";
		cout << _individual[i].get_nSexActs_casual_period()<<"\t";
		cout << _individual[i].get_nSexActs_sexworker_period()<<"\t";
		cout << endl;
	}
}


void Population::displayInfo_SexActs_Pairs_lastPeriod()
{
	// DISPLAY ALL SEX ACTS BY SEX ACT TYPES AND UID PAIRS
	
	for (int uid=0; uid<_size; uid++)
	{
		// Retrieve UIDs of sex partner during this period
		vector<unsigned long> UIDsexpartners = _individual[uid].get_UID_sexAct_period();
		
		// Retrieve number of sex acts
		vector<int> nSexActsPerUID = _individual[uid].get_UID_n_sexAct_period();
		
		int Np = UIDsexpartners.size(); // number of partners
		if (Np>0)
		{
			// Loop through all sex partners
			for (int p=0; p<Np; p++)
			{
				int nSexActsPerUID_p = nSexActsPerUID[p];
				
				string cswFlag = "";
				if (isCSW(UIDsexpartners[p])) cswFlag="$";
				
				cout << uid << " had "<<nSexActsPerUID_p<<" sex act(s) with ";
				cout << UIDsexpartners[p]<<cswFlag<<" ";
				
				if (nSexActsPerUID_p>0)
				{
					cout << "(sex type condom/lo/hi: "<< _individual[uid].get_UID_n_sexAct_Type0_period(p)<<"/";
					cout << _individual[uid].get_UID_n_sexAct_Type1_period(p)<<"/";
					cout << _individual[uid].get_UID_n_sexAct_Type2_period(p)<<")";
				}
				cout << endl;
			}
		}
	}
}



void Population::displayInfo_STI_SFincrease()
{
	int nrow = _STI_SFincrease.getNbRows();
	int ncol = _STI_SFincrease.getNbCols();
	
	string errmsg = "matrix _STI_SFincrease is not of correct dimension!";
	stopif(nrow!= _nSTImodelled || ncol!= _nSTImodelled, errmsg);
	
	cout << endl << " SUSCEPTIBILITY INCREASE CO-INFECTION MATRIX:"<<endl<<endl;
	cout << "(column:STI already infecting; row:STI exposed to)"<<endl<<endl;
	
	string colgap="\t";
	
	for (int i=0; i<nrow+1; i++)
	{
		for (int j=0; j<ncol+1; j++)
		{
			if (i==0 && j==0) cout << colgap;
			if (i==0 && j>0) cout << STInameString(_STI[j-1].get_name()) <<colgap;
			
			if (i>0 && j==0) cout <<  STInameString(_STI[i-1].get_name()) <<colgap;
			if (i>0 && j>0) cout <<  _STI_SFincrease(i-1,j-1) <<colgap;
		}
		cout << endl;
	}
	
}

vector<unsigned long> Population::getUIDs(){
	
	vector<unsigned long> x;
	
	for(unsigned long i=0; i<_size; i++){
		unsigned long u =_individual[i].get_UID();
		x.push_back(u);
	}
	return x;
}



dcDataFrame Population::export_to_dataframe(){
	
	
	// Initialize data frame:
	vector<unsigned long> uid = getUIDs();
	dcDataFrame df(uid,"uid");
	
	// Retrieve all requested data:
	
	vector<double> isalive;
	for(unsigned long i=0; i<_size; i++) isalive.push_back(_individual[i].isAlive());

	vector<double> gender;
	for(unsigned long i=0; i<_size; i++) gender.push_back(_individual[i].get_gender());

	vector<double> dateinpop;
	for(unsigned long i=0; i<_size; i++) dateinpop.push_back(_individual[i].get_dateInPopulation());

	vector<double> age;
	for(unsigned long i=0; i<_size; i++) age.push_back(_individual[i].get_age());
	
	vector<double> riskgrp;
	for(unsigned long i=0; i<_size; i++) riskgrp.push_back(_individual[i].get_riskGroup());
	
	vector<double> nCurrPartn;
	for(unsigned long i=0; i<_size; i++) nCurrPartn.push_back(_individual[i].get_nCurrSexPartner());
	
	vector<double> nCurrMaxPartn;
	for(unsigned long i=0; i<_size; i++) nCurrMaxPartn.push_back(_individual[i].get_nMaxCurrSexPartner());
	
	vector<double> nLifetimePartner;
	for(unsigned long i=0; i<_size; i++) nLifetimePartner.push_back(_individual[i].get_nLifetimePartner());
	
	vector<double> nCurrSpouse;
	for(unsigned long i=0; i<_size; i++) nCurrSpouse.push_back(_individual[i].get_nCurrSpouse());
	
	vector<double> nLifetimeSpouse;
	for(unsigned long i=0; i<_size; i++) nLifetimeSpouse.push_back(_individual[i].get_nLifetimeSpouse());
	
	vector<double> circum;
	for(unsigned long i=0; i<_size; i++) circum.push_back(_individual[i].get_isCircum());
	
	vector<double> singleDur;
	for(unsigned long i=0; i<_size; i++) singleDur.push_back(_individual[i].get_singleDuration());
	
	vector<double> age1sex;
	for(unsigned long i=0; i<_size; i++) age1sex.push_back(_individual[i].get_ageFirstSex());
	
	vector<double> age1partner;
	for(unsigned long i=0; i<_size; i++) age1partner.push_back(_individual[i].get_ageFirstPartner());
	
	vector<double> age1spouse;
	for(unsigned long i=0; i<_size; i++) age1spouse.push_back(_individual[i].get_ageFirstSpouse());
	
	vector<double> everVisitCSW;
	for(unsigned long i=0; i<_size; i++) everVisitCSW.push_back(_individual[i].get_ever_visited_CSW());
	
	vector<double> isPregnant;
	for(unsigned long i=0; i<_size; i++) isPregnant.push_back(_individual[i].get_isPregnant());
	
	vector<double> gestDur;
	for(unsigned long i=0; i<_size; i++) gestDur.push_back(_individual[i].get_gestationDuration());
	
	vector<double> nChild;
	for(unsigned long i=0; i<_size; i++) nChild.push_back(_individual[i].get_nChildBorn());
	
	vector<double> nSexActs_lifetime;
	for(unsigned long i=0; i<_size; i++) nSexActs_lifetime.push_back((double)(_individual[i].get_nSexActs_lifetime()));
	
	// Construct the data frame:
	
	df.addcol("dateinpop", dateinpop);
	df.addcol("isalive", isalive);
	df.addcol("gender", gender);
	df.addcol("age", age);
	df.addcol("riskgroup", riskgrp);
	df.addcol("nCurrPartn", nCurrPartn);
	df.addcol("nCurrMaxPartn", nCurrMaxPartn);
	df.addcol("nCurrSpouse", nCurrSpouse);
	df.addcol("nLifetimePartner", nLifetimePartner);
	df.addcol("nLifetimeSpouse", nLifetimeSpouse);
	df.addcol("circum", circum);
	df.addcol("singleDur", singleDur);
	df.addcol("age1sex", age1sex);
	df.addcol("age1partner", age1partner);
	df.addcol("age1spouse", age1spouse);
	df.addcol("everVisitCSW", everVisitCSW);
	df.addcol("isPregnant", isPregnant);
	df.addcol("gestDur", gestDur);
	df.addcol("nChild", nChild);
	df.addcol("nSexActs_lifetime", nSexActs_lifetime);
	
	// Partnerships:
	
	int maxNumberPartners = 5;
	
	for (int p=0; p<maxNumberPartners; p++)
	{
		vector<double> uid_p;
		vector<double> prtn_dur;
		//
		
		for(unsigned long i=0; i<_size; i++){
			
			unsigned long tmp = 0;
			double tmp_dur = 0;
			if (_individual[i].get_nCurrSexPartner()>p){
				tmp		= _individual[i].getPartnerUID(p);
				tmp_dur	= _individual[i].getPartnershipDuration()[p];
			}
			uid_p.push_back(tmp);
			prtn_dur.push_back(tmp_dur);
		}
		string header = "UIDpartner" + to_string(p+1);
		df.addcol(header, uid_p);
		header = "durPrtn" + to_string(p+1);
		df.addcol(header, prtn_dur);
	}
	
	// STI
	
	for (int sti=0; sti<_STI.size(); sti++){
		
		vector<double> sti_dur;
		vector<double> sti_sympt;
		vector<double> sti_treat;
		vector<double> sti_immun;
		vector<double> sti_vacc_date;
		vector<double> sti_IC;
		
		string stiname = STInameString(_STI[sti].get_name());
		
		for(unsigned long i=0; i<_size; i++){
			sti_dur.push_back(_individual[i].get_STIduration()[sti]);
			sti_sympt.push_back(_individual[i].get_STIsymptom()[sti]);
			sti_immun.push_back(_individual[i].get_STI_immunity()[sti]);
			sti_treat.push_back(_individual[i].get_STItreatDuration()[sti]);
			sti_vacc_date.push_back(_individual[i].get_STI_vacc_time()[sti]);
			sti_IC.push_back(_individual[i].STI_IC()[sti]);
		}
		string header = stiname + "duration";
		df.addcol(header, sti_dur);
		header = stiname + "sympt";
		df.addcol(header, sti_sympt);
		header = stiname + "treat";
		df.addcol(header, sti_treat);
		header = stiname + "immun";
		df.addcol(header, sti_immun);
		header = stiname + "vaccTime";
		df.addcol(header, sti_vacc_date);
		header = stiname + "_IC";
		df.addcol(header, sti_IC);
	}
	return df;
}





void Population::saveToCSVFile(string pathFile)
{
	/// SAVE THE WHOLE POPULATION (INCLUDING DEADS)
	/// TO A FILE FOR ANALYSIS OUTSIDE C++
	
	ofstream f(pathFile.c_str());
	
	int maxNumberPartners = 5;
	
	// ==== Headers ====
	
	f << "UID,dateInPopulation,isAlive,gender,age,";
	f << "riskGroup,nCurrSexPartner,nMaxCurrSexPartner,";
	f << "nLifetimeSexPartner,nCurrSpouse,nLifeTimeSpouse,";
	f << "isDivorced,isWidow,isCircum,singleDuration,";
	f << "ageFirstSex,ageFirstPartner,ageFirstSpouse,";
	f << "everVisitCSW,isPregnant,gestationDuration,nChildBorn,";
	
	for (int p=0; p<maxNumberPartners; p++)
		f << "partnerUID"<< p <<",";
	
	for (int p=0; p<maxNumberPartners; p++)
	{
		f << "pDuration"<< p;
		f <<",";
	}
	
	// STIs headers
	
	for (int i=0; i<_STI.size(); i++)
	{
		f << STInameString(_STI[i].get_name()) <<"duration,";
	}
	
	for (int i=0; i<_STI.size(); i++)
	{
		f << STInameString(_STI[i].get_name()) <<"symptom,";
	}
	
	// Vaccination headers
	
	for (int i=0; i<_STI.size(); i++)
	{
		f << STInameString(_STI[i].get_name()) <<"vacc";
		if (i<_STI.size()-1) f<<",";
	}
	
	
	// ==== (end of headers) ====
	
	f << endl;
	
	
	// ==== DATA ====
	
	for (int i=0; i<_size; i++)
	{
		f <<_individual[i].get_UID() << ","
		<< _individual[i].get_dateInPopulation() << ","
		<< _individual[i].isAlive()<< ","
		<< _individual[i].get_gender()<< ","
		<< _individual[i].get_age()<< ","
		<< _individual[i].get_riskGroup()<< ","
		<< _individual[i].get_nCurrSexPartner()<< ","
		<< _individual[i].get_nMaxCurrSexPartner()<< ","
		<< _individual[i].get_nLifetimePartner()<< ","
		<< _individual[i].get_nCurrSpouse()<< ","
		<< _individual[i].get_nLifetimeSpouse()<< ","
		<< _individual[i].get_isDivorced()<<","
		<< _individual[i].get_isWidow()<<","
		<< _individual[i].get_isCircum() <<","
		<< _individual[i].get_singleDuration()<<","
		<< _individual[i].get_ageFirstSex() <<","
		<< _individual[i].get_ageFirstPartner()<<","
		<< _individual[i].get_ageFirstSpouse()<<","
		<< _individual[i].get_ever_visited_CSW() <<","
		<< _individual[i].get_isPregnant() <<","
		<< _individual[i].get_gestationDuration() <<","
		<< _individual[i].get_nChildBorn() <<","
		;
		
		
		for (int p=0; p<maxNumberPartners; p++)
		{
			if (_individual[i].get_nCurrSexPartner()>p)
			{
				f << _individual[i].getPartnerUID(p);
			}
			f << ",";
		}
		
		for (int p=0; p<maxNumberPartners; p++)
		{
			if (_individual[i].get_nCurrSexPartner()>p)
			{
				f << _individual[i].getPartnershipDuration()[p];
			}
			
			f<<",";
		}
		
		// STIs data
		for (int sti=0; sti<_STI.size(); sti++)
			f << _individual[i].get_STIduration()[sti] <<",";
		
		for (int sti=0; sti<_STI.size(); sti++)
			f << _individual[i].get_STIsymptom()[sti] <<",";
		
		// Vaccination data
		for (int sti=0; sti<_STI.size(); sti++){
			f << _individual[i].get_STI_immunity()[sti];
			if (sti<_STI.size()-1) f << ",";
		}
		
		f << endl;
	}
}




void Population::save_outputs_demog(string pathFolder)
{
	/// SAVE DEMOGRAPHIC OUTPUTS FROM CURRENT POPULATION
	
	
	// === Age distribution ===
	
	vector<double> ageBreaks;
	for (int i=0; i<99; i++)  ageBreaks.push_back(1.0*i);
	
	vector<double> AD = census_ageDistribution(ageBreaks);
	
	ageBreaks.pop_back();
	dcMatrix M(ageBreaks);
	M.addColVector(AD);
	
	M.WriteToFileCSV(pathFolder+"census_ageDistribution.out");
	
	//DEBUG
	//exit(1);
	
	
	// (add more here - if needed)
	
}



void Population::save_outputs_demog(string pathFolder,
									vector<double> ageBreaks,
									unsigned int iMC,
									unsigned int idate)
{
	/// SAVE DEMOGRAPHIC OUTPUTS FROM CURRENT POPULATION
	
	
	// === Age distribution ===
	
	// makes sure the breaks are broad enough
	ageBreaks.push_back(999.9);
	
	vector<double> AD = census_ageDistribution(ageBreaks);
	ageBreaks.pop_back();
	dcMatrix M(ageBreaks);
	M.addColVector(AD);
	
	string fileout = pathFolder+"ageDistrib_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	
	M.WriteToFileCSV(fileout);
	
	// (add more here - if needed)
	
}



void Population::save_outputs_prtnr(string pathFolder)
{
	/// SAVE PARTNERSHIP OUTPUTS FROM CURRENT POPULATION
	
	
	// === Age gap distribution ===
	
	vector<double> ageBreaks;
	for (int i=-80; i<80; i++)  ageBreaks.push_back(1.0*i);
	
	vector<double> AGD = census_ageGapDistribution(ageBreaks);
	ageBreaks.pop_back();
	dcMatrix M(ageBreaks);
	M.addColVector(AGD);
	
	M.WriteToFileCSV(pathFolder+"census_ageGapDistrib.out");
	
	
	
	// === Single Ratio ===
	
	
	double singleratio_f = census_ratioSingles(female);
	double singleratio_m = census_ratioSingles(male);
	
	
	string ff = pathFolder+"census_ratioSingles_f.out";
	string fm = pathFolder+"census_ratioSingles_m.out";
	
	ofstream sr_f(ff.c_str());
	ofstream sr_m(fm.c_str());
	
	sr_f<<singleratio_f<<endl;
	sr_m<<singleratio_m<<endl;
	
	
	
	// (add more here - if needed)
	
}



void Population::save_outputs_prtnr(string pathFolder,
									vector<double> ageBreaks,
									unsigned int iMC,
									unsigned int idate)
{
	/// SAVE PARTNERSHIP OUTPUTS FROM CURRENT POPULATION
	
	
	// makes sure the breaks are broad enough
	ageBreaks.insert(ageBreaks.begin(), -999.9);
	ageBreaks.push_back(999.9);
	
	
	// === Age gap distribution ===
	
	vector<double> AGD = census_ageGapDistribution(ageBreaks);
	ageBreaks.pop_back();
	dcMatrix M(ageBreaks);
	M.addColVector(AGD);
	
	string fileout = pathFolder+"ageGapDistrib_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	M.WriteToFileCSV(fileout);
	
	// === Single Ratio ===
	
	double singleratio_f = census_ratioSingles(female);
	double singleratio_m = census_ratioSingles(male);
	
	string filename = pathFolder+"singleRatio_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	ofstream thefile(filename.c_str());
	thefile<<"female,"<<singleratio_f<<endl;
	thefile<<"male,"<<singleratio_m<<endl;
	
	// old stuff - delete when sure it's not used anymore:
	//	string ff = pathFolder+"singleRatio_f_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	//	string fm = pathFolder+"singleRatio_m_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	//
	//	ofstream sr_f(ff.c_str());
	//	ofstream sr_m(fm.c_str());
	//
	//	sr_f<<singleratio_f<<endl;
	//	sr_m<<singleratio_m<<endl;
	
	
	// ============================
	// (add more here - if needed)
	
}



void Population::save_outputs_sex(string pathFolder)
{
	/// SAVE SEXUAL BEHAVIOUR OUTPUTS FROM CURRENT POPULATION
	
	// === Age First sex distribution ===
	
	vector<double> ageBreaks;
	for (int i=0; i<90; i++)  ageBreaks.push_back(1.0*i);
	
	vector<double> AFSD_f = census_ageFirstSexDistribution(ageBreaks,female);
	vector<double> AFSD_m = census_ageFirstSexDistribution(ageBreaks,male);
	
	ageBreaks.pop_back();
	
	dcMatrix M_f(ageBreaks);
	M_f.addColVector(AFSD_f);
	M_f.WriteToFileCSV(pathFolder+"census_ageFirstSexDistribution_f.out");
	
	
	dcMatrix M_m(ageBreaks);
	M_m.addColVector(AFSD_m);
	M_m.WriteToFileCSV(pathFolder+"census_ageFirstSexDistribution_m.out");
	
	
	
	// === Age Gap First sex and First spouse distribution ===
	
	vector<double> ageGapBreaks;
	for (int i=-80; i<80; i++)  ageGapBreaks.push_back(1.0*i);
	
	vector<double> AGFSD_f = census_ageGapFirstSexSpouseDistribution(ageGapBreaks,female);
	vector<double> AGFSD_m = census_ageGapFirstSexSpouseDistribution(ageGapBreaks,male);
	
	ageGapBreaks.pop_back();
	
	dcMatrix MG_f(ageGapBreaks);
	MG_f.addColVector(AGFSD_f);
	MG_f.WriteToFileCSV(pathFolder+"census_ageGapFirstSexSpouseDistribution_f.out");
	
	dcMatrix MG_m(ageGapBreaks);
	MG_m.addColVector(AGFSD_m);
	MG_m.WriteToFileCSV(pathFolder+"census_ageGapFirstSexSpouseDistribution_m.out");
	
	
	
	// === Lifetime number of sex partners distribution ===
	
	vector<double> nBreaks;
	for (int i=0; i<999; i++)  nBreaks.push_back(1.0*i);
	
	vector<double> LSPD_f = census_nLifeSexPrtnrDistrib(nBreaks, female);
	vector<double> LSPD_m = census_nLifeSexPrtnrDistrib(nBreaks, male);
	
	nBreaks.pop_back();
	
	dcMatrix ML_f(nBreaks);
	ML_f.addColVector(LSPD_f);
	ML_f.WriteToFileCSV(pathFolder+"census_nLifeSexPrtnrDistrib_f.out");
	
	dcMatrix ML_m(nBreaks);
	ML_m.addColVector(LSPD_m);
	ML_m.WriteToFileCSV(pathFolder+"census_nLifeSexPrtnrDistrib_m.out");
	
}



void Population::save_outputs_sex(string pathFolder,
								  vector<double> ageBreaks_f,
								  vector<double> ageBreaks_m,
								  vector<double> ageGapBreaks_f,
								  vector<double> ageGapBreaks_m,
								  vector<double> nBreaks_f,
								  vector<double> nBreaks_m,
								  vector<double> ageBreaks_malesVisitCSW,
								  unsigned int iMC, unsigned int idate)
{
	/// SAVE SEXUAL BEHAVIOUR OUTPUTS FROM CURRENT POPULATION
	
	// makes sure the breaks are broad enough
	
	ageBreaks_f.push_back(999.9);
	ageBreaks_m.push_back(999.9);
	ageGapBreaks_f.insert(ageGapBreaks_f.begin(), -999.9);
	ageGapBreaks_m.insert(ageGapBreaks_m.begin(), -999.9);
	ageGapBreaks_f.push_back(999.9);
	ageGapBreaks_m.push_back(999.9);
	ageBreaks_malesVisitCSW.push_back(999.9);
	
	nBreaks_f.push_back(99999.9);
	nBreaks_m.push_back(99999.9);
	
	
	string fileout_f;
	string fileout_m;
	
	// === Age First sex distribution ===
	
	
	
	vector<double> AFSD_f = census_ageFirstSexDistribution(ageBreaks_f,female);
	vector<double> AFSD_m = census_ageFirstSexDistribution(ageBreaks_m,male);
	
	ageBreaks_f.pop_back();
	ageBreaks_m.pop_back();
	
	fileout_f = pathFolder+"age1sex_f_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	fileout_m = pathFolder+"age1sex_m_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	
	dcMatrix M_f(ageBreaks_f);
	M_f.addColVector(AFSD_f);
	M_f.WriteToFileCSV(fileout_f);
	
	dcMatrix M_m(ageBreaks_m);
	M_m.addColVector(AFSD_m);
	M_m.WriteToFileCSV(fileout_m);
	
	
	
	// === Age Gap First sex and First spouse distribution ===
	
	if (ageGapBreaks_f.size()>0) {
		vector<double> AGFSD_f = census_ageGapFirstSexSpouseDistribution(ageGapBreaks_f,female);
		vector<double> AGFSD_m = census_ageGapFirstSexSpouseDistribution(ageGapBreaks_m,male);
		
		ageGapBreaks_f.pop_back();
		ageGapBreaks_m.pop_back();
		
		fileout_f = pathFolder+"ageGap1SexSpouseDistrib_f_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
		fileout_m = pathFolder+"ageGap1SexSpouseDistrib_m_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
		
		dcMatrix MG_f(ageGapBreaks_f);
		MG_f.addColVector(AGFSD_f);
		MG_f.WriteToFileCSV(fileout_f);
		
		dcMatrix MG_m(ageGapBreaks_m);
		MG_m.addColVector(AGFSD_m);
		MG_m.WriteToFileCSV(fileout_m);
	}
	
	
	// === Lifetime number of sex partners distribution ===
	
	
	vector<double> LSPD_f = census_nLifeSexPrtnrDistrib(nBreaks_f, female);
	vector<double> LSPD_m = census_nLifeSexPrtnrDistrib(nBreaks_m, male);
	
	nBreaks_f.pop_back();
	nBreaks_m.pop_back();
	
	fileout_f = pathFolder+"lftNP_f_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	fileout_m = pathFolder+"lftNP_m_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	
	
	dcMatrix ML_f(nBreaks_f);
	ML_f.addColVector(LSPD_f);
	ML_f.WriteToFileCSV(fileout_f);
	
	dcMatrix ML_m(nBreaks_m);
	ML_m.addColVector(LSPD_m);
	ML_m.WriteToFileCSV(fileout_m);
	
	
	
	// === Proportion males visiting CSW ===
	
	double MVC = census_maleVisitCSW();
	
	string fMVC = pathFolder+"visitCSW_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	ofstream MVC_f(fMVC.c_str());
	
	MVC_f<<MVC<<endl;
	
	
	// === Age distribution males visiting CSW ===
	
	double maxDurationSinceLastVisit = 1.0;
	
	vector<double> AMVC = census_ageMalesVisitCSWDistribution(ageBreaks_malesVisitCSW,
															  maxDurationSinceLastVisit);
	
	ageBreaks_malesVisitCSW.pop_back();
	
	dcMatrix MAMVC(ageBreaks_malesVisitCSW);
	MAMVC.addColVector(AMVC);
	
	string fileout =pathFolder+"visitCSWage_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
	
	MAMVC.WriteToFileCSV(fileout);
	
}




void Population::save_outputs_epi(string pathFolder,
								  vector<STIname> STIname_targeted,
								  vector<double> agebreaks_HIVprev_f,
								  vector<double> agebreaks_HIVprev_m,
								  unsigned int iMC, unsigned int idate)
{
	/// Save epidemiological outputs from model
	
	// STI global prevalences
	
	if (STIname_targeted.size()>0){
		
		vector<double> STIprev;
		for (int i=0;i<STIname_targeted.size(); i++)
		{
			STIprev.push_back(STI_prevalence(STIname_targeted[i]));
		}
		string fileout_sti = pathFolder+"STI_prev_global_mc"+int2string(iMC)+".out";
		vectorToCSVFile_Row(STIprev, fileout_sti);
	}
	
	// HIV prevalence by age and gender
	
	if(agebreaks_HIVprev_f.size()>0)
	{
		vector<double> hiv_f = STI_prevalence_by_age(HIV, female, agebreaks_HIVprev_f);
		agebreaks_HIVprev_f.pop_back();
		dcMatrix HIVprev_f(agebreaks_HIVprev_f);
		HIVprev_f.addColVector(hiv_f);
		string fileout_f = pathFolder+"HIV_prev_age_f_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
		HIVprev_f.WriteToFileCSV(fileout_f);
	}
	
	if(agebreaks_HIVprev_m.size()>0)
	{
		vector<double> hiv_m = STI_prevalence_by_age(HIV, male, agebreaks_HIVprev_m);
		agebreaks_HIVprev_m.pop_back();
		dcMatrix HIVprev_m(agebreaks_HIVprev_m);
		HIVprev_m.addColVector(hiv_m);
		string fileout_m = pathFolder+"HIV_prev_age_m_D"+to_string(idate)+"mc"+int2string(iMC)+".out";
		HIVprev_m.WriteToFileCSV(fileout_m);
	}
	
}



void Population::FileOutput(string pathFile)
{
	/// CSV FILE FOR FURTHER ANALYSIS OUTSIDE C++
	
	ofstream f(pathFile.c_str());
	
	// Headers
	f<<"p1,p2,isPartner,risk1,risk2,anySTI"<<endl;
	
	
	for (int i=0; i<_size; i++){
		// Retrieves gender of left-hand partner
		string g_i="F";
		if (_individual[i].get_gender()==male) g_i="M";
		
		unsigned long uid = _individual[i].get_UID();
		
		int     nConc = _individual[i].get_nCurrSexPartner();
		
		int     rg1 = _individual[i].get_riskGroup();
		bool    anySTI = _individual[i].STI_anyInfection();
		
		if (nConc>0){
			for (int p=0; p<nConc; p++){
				unsigned long uid_p = _individual[i].getPartnerUID(p);
				int rg2 = _individual[uid_p].get_riskGroup();
				
				string g_p="F";
				if (_individual[uid_p].get_gender()==male) g_p="M";
				
				// Edge between partners
				f << g_i << uid << "," << g_p << uid_p <<",1,"<<rg1<<","<<rg2<<",";
				f << anySTI << endl;
			}
		}
	}
}



unsigned long Population::getUIDfemale(unsigned long uid1, unsigned long uid2)
{
	Gender g1 = _individual[uid1].get_gender();
	Gender g2 = _individual[uid2].get_gender();
	
	string errmsg = to_string(uid1) + " and " + to_string(uid2) + " same gender (must be different)!";
	stopif(g1==g2,errmsg);
	
	unsigned long uid_f = uid1;
	if (g1==male) uid_f = uid2;
	
	return uid_f;
}


unsigned long Population::getUIDmale(unsigned long uid1, unsigned long uid2)
{
	unsigned long uid_m = uid1;
	unsigned long uid_f = getUIDfemale(uid1, uid2);
	
	if (uid_f==uid1) uid_m = uid2;
	return uid_m;
}


void Population::treat_indiv(unsigned long uid, STIname sti)
{
	_individual[uid].treat(sti);
}


void Population::cure_indiv(unsigned long uid, STIname sti)
{
	/// Completely cure STI from this individual
	
	_individual[uid].set_STIduration(sti, 0.00);
	_individual[uid].set_STIsymptom(sti, false);
	_individual[uid].set_STItreatDuration(0.00, sti);
	_individual[uid].clear_STI_secondary_cases(sti);
}


void Population::vaccinate_indiv(unsigned long uid, STIname stiname)
{
	/// Vaccinate an individual against a given STI
	/// (not necessarily successful though)
	
	// -- Checks --
	stopif(!_individual[uid].isAlive(), "Cannot vaccinate a dead individual!");
	// ------------
	
	int i_sti = positionSTIinVector(stiname, _STI);
	
	// Record this individual has been vaccinated (irrespective of success)
	_individual[uid].set_STI_vacc(i_sti, true);
	_individual[uid].set_STI_vacc_time(i_sti, _simulationTime);
	
	// Determine if vaccine will be succesfull for that individual:
	double p = _STI[i_sti].get_proba_vaccineFailure();
	double vaxSuccess = 1.0;
	double u = uniform01();
	if(u<p) vaxSuccess = 0.0;
	_individual[uid].set_STI_immunity(i_sti, vaxSuccess);
		
	//DEBUG
	//cout<<endl<<"VAX DEBUG::: UID "<<uid<< STInameString(stiname)<<"vax: "<<vaxSuccess<<endl;
	// ----
}



void Population::set_STIsusceptFactor(unsigned long uid, STIname stiname, double x){

	int sti_i = positionSTIinVector(stiname, _STI);
	_individual[uid].set_STIsusceptFactor(sti_i, x);
}


void Population::set_STI_immunity(unsigned long uid, STIname stiname, double immunity){
	
	/// Set the immunitiy level of a given individual for a STI
	
	int sti_i = positionSTIinVector(stiname, _STI);
	_individual[uid].set_STI_immunity(sti_i, immunity);
}


void Population::set_STI_MTCT(unsigned long uid, vector<bool> mtct){
	_individual[uid].set_STI_MTCT(mtct);
}

void Population::set_STI_MTCT(unsigned long uid, STIname stiname, bool mtct){
	_individual[uid].set_STI_MTCT(stiname, mtct);
}


void Population::debug_check_partnerUID()
{
	//cout <<endl<<"debug_check_partnerUID"<<endl;
	
	for(int i=0;i<_size; i++){
		if (_individual[i].isAlive()){
			for(int j=0;j<_individual[i].getPartnerUID().size();j++){
				unsigned long a =_individual[i].getPartnerUID(j);
				stopif(a>=_size, "DEBUG");
			}
		}
	}
}



/* *****************************************
 **********   OUTSIDE CLASS ****************
 ****************************************** */



vector<unsigned long> convertTo_unsignedLong(vector<double> x, bool deleteFirst)
{
	vector<unsigned long> x_converted;
	
	int i_start = 0;
	if (deleteFirst) i_start=1;
	
	for	(int p=i_start;p<x.size();p++) x_converted.push_back((unsigned long)(x[p]));
	
	return(x_converted);
}

vector<int> convertTo_int(vector<double> x, bool deleteFirst)
{
	vector<int> x_converted;
	
	int i_start = 0;
	if (deleteFirst) i_start=1;
	
	for	(int p=i_start;p<x.size();p++) x_converted.push_back((int)(x[p]));
	
	return(x_converted);
}

vector<Gender> convertTo_Gender(vector<double> x, bool deleteFirst)
{
	vector<Gender> x_converted;
	
	int i_start = 0;
	if (deleteFirst) i_start=1;
	
	for	(int p=i_start;p<x.size();p++) x_converted.push_back((Gender)(x[p]));
	
	return(x_converted);
}

double		expGauss(double x, double mean, double var)
{
	// Gaussian density NOT normalized
	return exp(-(x-mean)*(x-mean)/2/var/var);
}

double		expGauss2(double x, double k1, double k2)
{
	return k1*exp(-k2*x*x);
}

