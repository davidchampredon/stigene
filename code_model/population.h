/*
 *  population.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-07-20.
 *  Copyright 2013 David Champredon. All rights reserved.
 *
 */

#ifndef population_h
#define population_h

#include "individuals.h"
#include "dcTools.h"
#include "dcDataFrame.h"


class Population
{
	// ==== Fundamental members ====
	
	unsigned long		_size;			// size of the population (all individuals, including the deads)
	vector<Individual>	_individual;	// Individual, atomic component of a population
	
	double				_maxAge;		// Maximum age of individuals (force death if this age is reached)
	
	double				_simulationTime;// keeps track of time when embedded in a simulation
	double				_timeStep;		// time step when embedded in a simulation
	
	// ==== Demographics ====
	
	double				_birthRate;	
	double				_infantMortality;	// Mortality annual rate of children <1 yrs-old
	double				_childMortality;	// Mortality annual rate of children <5 yrs-old and >1 yrs-old
	
	double				_probaPregnantPerSexAct; // probability to become pregnant after one sex act
	double				_maxPregnantAge;
	
	
	
	double				_deathParam_Weibull_shape;		// Natural death Weibull param hazard rate
	double				_deathParam_Weibull_scale;
	double				_deathParam_Weibull_shape_hiv;	// HIV-induced Weibull param
	double				_deathParam_Weibull_scale_hiv;
	
	double				_proportion_circum;	// Proportion of males circumcised
	
	
	
	
	// ==== PARTNERSHIPS ====
	
	/* -- GET RID OF THIS -- */
	dcMatrix				_PartnerMatrix;				// dcMatrix (Nf,Nm) with all partnerships
	/* -- */

	
	dcMatrix				_partnershipsMatrix;		// dcMatrix of UIDs of all partnerships; Nx2 matrix; N: number of partnerships; line = pair of UID of a partnership; 1st column: UID female, 2nd column: UID male

	vector<unsigned long>	_all_UID_femaleUID;		// UID of females (rows) associated with _PartnerMatrix
	vector<unsigned long>	_all_UID_maleUID;		// UID of males (columns) associated with _PartnerMatrix
	
	unsigned long		_totalNumberPartnerships;
	unsigned long		_totalNumberSpousalPartnerships;
	

	// -- Matching probabilities parameter--
	
	double				_formation_MaxRate;			// Maximum consideration rate for partnership formation by females
	vector<double>		_formation_RiskGroup;		// Parameters defining proba good match, riskGroup component
	vector<double>		_formation_PartnDeficit;	// Parameters defining proba good match, partner deficit component
	
	// Age at new partnership formation
	
	double				_formation_age_fullstart;		// Age when female is fully age-attractive
	double				_formation_age_pivot;			// Parameter representing age where female is half age-attractive (logistic fct)
	double				_formation_age_shape;			// Shape parameter for female age attractiveness
	double				_formation_age_fmin;			// Shape parameter minimum effect of age proba component
	// ----
	
	// Age Gap = AgeMale - AgeFemale
	double				_formation_agegap_mean;			// mean age gap of the population
	double				_formation_agegap_var;			// variance of age gap at the population level
	
	vector<double>		_formation_agegap_shape;		// shape parameters
	double				_formation_agegap_fmin;			// Shape parameter minimum effect of age gap proba component
	
	//double				_formation_correl_Age_AgeGap;	// Correlation b/w age at new formation and age gap for females
	//double				_formation_shapeAge;			// ("alpha") shape parameter for age dependant new formation
	
	// STI symptom impact
	double				_formation_STIsymptom_f;	// Relative reduction of proba formation if female has any STI symptom
	double				_formation_STIsymptom_m;	// Relative reduction of proba formation if male has any STI symptom

		// -- Spousal progression
	double				_spousalProgress_maxRate;		// Max rate (unit=year^-1) of spousal progression
	double				_spousalProgress_meanAge_f;		// Mean age of female entering a new spousal relationship
	double				_spousalProgress_varAge_f;		// Variance of age of female entering a new spousal relationship
	double				_spousalProgress_meanGap;		// Mean age gap of a new spousal relationship
	double				_spousalProgress_varGap;		// Variance of age gap of a new spousal relationship
	
	double				_spousalProgress_durationK1;	// Shape parameter for duration based rate
	double				_spousalProgress_durationK2;	// Shape parameter for duration based rate	
	
	double				_spousalProgress_meanDiffAgeGap;	// Mean difference b/w existing spouses' age gap and age gap ofa new spouse candidate
	double				_spousalProgress_varDiffAgeGap;		// Var of diff b/w existing spouses' age gap and age gap ofa new spouse candidate
	
		// -- Dissolution

	double				_dissolution_MaxRate;	// Maximum consideration rate for partnership dissolution
	double				_dissolution_spouse;	// Reduction factor of risk to dissolve when in spousal union
	vector<double>		_dissolution_RiskGroup;	// Parameters for risk-based dissolution component
	vector<double>		_dissolution_duration;	// Parameters for duration-based dissolution component
	vector<double>		_dissolution_age;		// Parameters for age-based (both female and male) dissolution component
	vector<double>		_dissolution_PartnerDeficit;// Parameters for partner deficit dissolution component
	double				_dissolution_STI_symptom;	// Relative reduction of proba to dissolve if any partner is symptomatic of any STI
	double				_dissolution_ageConcPartn;	// shape parameter for link b/w age and # of concurrent partnership
		
	
	// ==== CSW ====
	
	double				_CSW_maxRecruitmentRate;		// Maximum recruitment rate of new CSWs 
	double				_CSW_recruitment_saturPrm_1;	// Saturation parameter for CSW recruitment
	double				_CSW_recruitment_saturPrm_2;
	double				_CSW_recruitment_ageMin;		// Minimum age to start commercial sex
	double				_CSW_recruitment_ageMax;		// Maximum age to start commercial sex
	double				_CSW_ageMax;					// Maximum age for practicing CSW
	double				_CSW_cessationRate;				// Drop out rate of commercial sex
		
	

	// ==== Sexual Behaviour ====
	
	double				_ageSexMin;	// Minimum age sexually active
	double				_ageSexMax;	// Maximum age sexually active
	
	int					_maxRiskGroup;	// Maximum integer value for risk group (except CSW)
	int					_CSWriskGroup;	// Risk group value for commercial sex workers (should be distinctively high compared to other)
	vector<double>		_propRiskGroup;	// Proportion of the whole population in each risk group (except CSW)
	
	vector<double>		_nMaxCurrSexPartner_param_f; // Parameters defining the maximum number of concurrent partner for a given female
	vector<double>		_nMaxCurrSexPartner_param_m; // Parameters defining the maximum number of concurrent partner for a given male
	
	
	double				_sexAct_maxRate_male;		// Maximum rate of sexual activity for males (unit: sex acts per year)
	double				_sexAct_maxRate_female;		// Maximum rate of sexual activity for females (unit: sex acts per year)

	vector<double>		_sexAct_reduce_age_param;			// parameters involved in reducing sexual activity according to age
	vector<double>		_sexAct_reduce_risk_param;			// parameters involved in reducing sexual activity according to risk group
	vector<double>		_sexAct_reduce_STIsymptom_param;	// parameters involved in reducing sexual activity according to STI symptoms
	vector<double>		_sexAct_reduce_nPartner_param;		// parameters involved in reducing sexual activity according number of current partners

	double				_sexAct_proba_distribute_partnerTypes_prefSpouse; // Preference parameter for sex acts with spouses (rather than casual)
	vector<double>		_sexAct_proba_distribute_partnerTypes_sexWorkParam; // Preference parameter for sex acts with sex workers
	
	vector<double>		_sexAct_proba_sexWorker_param;	// parameters defining the probability of intercourse with sex workers
	
	double				_sexAct_CostSexWork_reduction;	// reduction factor of number of sexual acts due to commercial sex costs
	
	vector<double>		_sexAct_condom_param;		// Parameters for condom use
	double				_sexAct_TypeLowHigh;		// Ratio of low risk sex vs high risk, for all sex acts WITHOUT condom
	
	// ===== STI =====
	
	vector<STI>			_STI;				// STIs modelled ; template for each Individual's STI
	unsigned long		_nSTImodelled;		// Number of STIs modelled in the population
	dcMatrix			_STI_SFincrease;	// dcMatrix defining the increase factor of the susceptibility factor when infected with another STI
	vector<double>		_RebHIV;			// HIV infectivity rebound due to STIs co-infections
	
	vector< vector<unsigned long> > _secondary_cases; // Number of all secondary cases for every STIs (size of vector increases as simulation runs)
	
	vector<unsigned long>	_STI_mtct_cumcount;	// Cumulative count of MTCT events, for each STI
	
	
	// Book keepers
	
	vector<unsigned long>	_UID_pot_preg; // UIDs of all potential pregenant females
	
	
	
	// Record calculation details
	vector< vector<double> > _rec_sexact;

	
	// -----------------
	// Private functions
	// -----------------
		
	void		sexAct_distrib_within_prtnrType(unsigned long uid,
												string prtnrType,
												int nPrtnr, int nsexact,
												bool saveTraceFile);

	void		sexAct_choose_CSW(unsigned long uid,
								   int nsexact,
								   bool saveTraceFile);
	
public:
	
	
	// === Constructors ===
	
	Population(){}
	Population(unsigned long size);	
	
	
	// === Pseudo-Constructors ===
	
	
	void initFromFile(string pathFile, string STI_filename,
					  string _STI_SFincrease_filename, string _RebHIV_filename);
	
	void initFromFile2(string pathFile, string STI_filename,
					   string _STI_SFincrease_filename, string _RebHIV_filename,
					   string STI_treatment_filename,
					   string STI_vaccine_filename);
	
	void initPartnershipDuration();
	void initSingleDuration();
	
	
	void setup_for_simulation(unsigned long founder_size,
							  double founder_female_ratio,
							  double founder_prop_csw,
							  string folder_inputs,
							  string file_STI_features,
							  string file_STI_SFincrease,
							  string file_STI_HIVrebound,
							  string file_STI_treatment,
							  string file_STI_vaccine,
                              string file_STI_genome,
							  bool debugInfo);
	
	void set_and_check_UID();
	
	
	void create_founder_population(unsigned long size, double female_ratio, double prop_csw);
	
	
	
	// =====================
	// === Get Functions ===
	// =====================
	
	
	Individual		getIndividual(unsigned long uid) {return _individual[uid];}
	unsigned long	get_nSTImodelled() {return _nSTImodelled;}
	double			get_ageSexMin() {return _ageSexMin;}
	dcMatrix		get_partnershipsMatrix() { return _partnershipsMatrix;}
	
	double			get_simulationTime() { return _simulationTime;}
	
	unsigned long	get_size() {return _size;}
		
	vector<double>	get_age();
	double			getAge(int i) {return _individual[i].get_age();}
	double			getAgeGap(unsigned long uid1, unsigned long uid2);
	
	vector<unsigned long> getUID(Gender g); // returns all UIDs of Gender g
	vector<unsigned long> getUID_available(Gender g); // returns all UIDs of Gender g that are available for a new partnership
	
	int				get_maxRiskGroup() {return _maxRiskGroup;}
	
	double			get_birthRate() {return _birthRate;}
	double			get_infantMortality() {return _infantMortality;}
	double			get_childMortality() {return _childMortality;}
	
	double			get_probaPregnantPerSexAct() {return _probaPregnantPerSexAct;}
	
	double			get_deathParam_Weibull_shape() {return _deathParam_Weibull_shape;}
	double			get_deathParam_Weibull_scale() {return _deathParam_Weibull_scale;}
	double			get_deathParam_Weibull_shape_hiv() {return _deathParam_Weibull_shape_hiv;}
	double			get_deathParam_Weibull_scale_hiv() {return _deathParam_Weibull_scale_hiv;}
	
	double			get_formation_MaxRate() { return _formation_MaxRate;}
	vector<double>	get_formation_RiskGroup() {return _formation_RiskGroup;}
	vector<double>	get_formation_PartnDeficit() {return _formation_PartnDeficit;}
	
	double			get_formation_age_fullstart() {return _formation_age_fullstart;}
	double			get_formation_age_pivot() {return _formation_age_pivot;}
	double			get_formation_age_shape() {return _formation_age_shape;}
	double			get_formation_age_fmin() {return _formation_age_fmin;}
	
//	double			get_formation_meanAge_female() {return _formation_meanAge_female;}
//	double			get_formation_varAge_female() {return _formation_varAge_female;}
	double			get_formation_agegap_mean() {return _formation_agegap_mean;}
	double			get_formation_agegap_var() {return _formation_agegap_var;}
//	double			get_formation_correl_Age_AgeGap() {return _formation_correl_Age_AgeGap;}
//	double			get_formation_shapeAge() {return _formation_shapeAge;}
	
	double			get_spousalProgress_maxRate() {return _spousalProgress_maxRate;}
	double			get_spousalProgress_meanAge_f() {return _spousalProgress_meanAge_f;}	
	double			get_spousalProgress_varAge_f() {return _spousalProgress_varAge_f;}	
	double			get_spousalProgress_meanGap() {return _spousalProgress_meanGap;}			
	double			get_spousalProgress_varGap() {return _spousalProgress_varGap;}			
	
	double			get_spousalProgress_durationK1() {return _spousalProgress_durationK1;}		
	double			get_spousalProgress_durationK2() {return _spousalProgress_durationK2;}			
	
	double			get_spousalProgress_meanDiffAgeGap() {return _spousalProgress_meanDiffAgeGap;}		
	double			get_spousalProgress_varDiffAgeGap() {return _spousalProgress_varDiffAgeGap;}
	
	double			get_dissolution_MaxRate() { return _dissolution_MaxRate;}
	vector<double>	get_dissolution_duration() { return _dissolution_duration;}
	vector<double>	get_dissolution_RiskGroup() { return _dissolution_RiskGroup;}
	vector<double>	get_dissolution_age() { return _dissolution_age;}
	vector<double>	get_dissolution_PartnerDeficit() { return _dissolution_PartnerDeficit;}
	double			get_dissolution_ageConcPartn() { return _dissolution_ageConcPartn;}
	
	vector<double>	get_propRiskGroup(){return _propRiskGroup;}
	

	vector<STI>		get_STI() {return _STI;}
	
	double			get_STIsusceptFactor(STIname stiname, unsigned long uid);
	
	vector<unsigned long>	get_UID_pot_preg() {return _UID_pot_preg;}
	vector<unsigned long>	get_STI_mtct_cumcount() { return _STI_mtct_cumcount;}
	
	vector< vector<double> >	get_rec_sexact() {return _rec_sexact ;}
		
	
	// === Subset Functions ===
	
	vector<Individual>	subsetIndividual(Gender g); // return all individuals of gender g
	
	
	// ===============================
	// ===     Set functions       ===
	// ===============================

	void			set_simulationTime(double t) {_simulationTime = t;}
	void			set_timeStep(double t) {_timeStep = t;}
	
	void			set_maxRiskGroup(int r) {_maxRiskGroup = r;}
	
	// -- Demographics
	
	void			set_birthRate(double rate) {_birthRate=rate;}
	void			set_infantMortality(double rate) {_infantMortality=rate;}
	void			set_childMortality(double rate) {_childMortality=rate;}
	
	void			set_deathParam_Weibull_shape(double x) {_deathParam_Weibull_shape=x;}
	void			set_deathParam_Weibull_scale(double x) {_deathParam_Weibull_scale=x;}
	void			set_deathParam_Weibull_shape_hiv(double x) {_deathParam_Weibull_shape_hiv=x;}
	void			set_deathParam_Weibull_scale_hiv(double x) {_deathParam_Weibull_scale_hiv=x;}
	
	void			set_probaPregnantPerSexAct(double x) {_probaPregnantPerSexAct=x;}
	void			set_maxPregnantAge(double x) {_maxPregnantAge=x;}
	void			set_pregnant(unsigned long uid) {_UID_pot_preg=popElementValue(_UID_pot_preg, uid);	_individual[uid].set_isPregnant(true);}
	void			set_gestationDuration(unsigned long uid, double x){_individual[uid].set_gestationDuration(x);}
	void			increment_gestationDuration(unsigned long uid, double x)
					{_individual[uid].set_gestationDuration(x+_individual[uid].get_gestationDuration());}
	void			increment_nChildBorn(unsigned long uid);

	void			set_UID_pot_preg(vector<unsigned long> x){_UID_pot_preg = x;}
	
	
	
	// -- Formation
	
	void			set_formation_MaxRate(double rate) {_formation_MaxRate = rate;}
	
	void			set_formation_age_fullstart(double x) {_formation_age_fullstart=x;}
	void			set_formation_age_pivot(double x) {_formation_age_pivot=x;}
	void			set_formation_age_shape(double x) {_formation_age_shape=x;}
	void			set_formation_age_fmin(double x) {_formation_age_fmin=x;}
	
	void			set_formation_agegap_mean(double x) {_formation_agegap_mean=x;}
	void			set_formation_agegap_var(double x) {_formation_agegap_var=x;}
	void			set_formation_agegap_fmin(double x) {_formation_agegap_fmin=x;}
	
	void			set_formation_agegap_shape(vector<double> x) {_formation_agegap_shape=x;}
	void			set_formation_agegap_shape(double x, int i) {_formation_agegap_shape[i]=x;}
	
	void			set_formation_STIsymptom_f(double x) {_formation_STIsymptom_f=x;}
	void			set_formation_STIsymptom_m(double x) {_formation_STIsymptom_m=x;}

//	void			set_formation_AgeGap_female(double m_g, double v_g) {_formation_agegap_mean=m_g ; _formation_agegap_var=v_g;}
	
	void			set_formation_RiskGroup(vector<double> x) {_formation_RiskGroup=x;}
	void			set_formation_RiskGroup(double x, int i) {_formation_RiskGroup[i]=x;}
	
	void			set_formation_PartnDeficit(vector<double> x) {_formation_PartnDeficit=x;}
	void			set_formation_PartnDeficit(double x, int i) {_formation_PartnDeficit[i]=x;}

	
	
	
	// -- Spousal progress
	
	
	void			set_spousalProgress_maxRate(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_meanAge_f(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_varAge_f(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_meanGap(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_varGap(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_durationK1(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_durationK2(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_meanDiffAgeGap(double x){_spousalProgress_maxRate=x;}
	void			set_spousalProgress_varDiffAgeGap(double x){_spousalProgress_maxRate=x;}
	

	
	// -- Dissolution
	
	void			set_dissolution_MaxRate(double x) {_dissolution_MaxRate=x;}
	void			set_dissolution_spouse(double x) {_dissolution_spouse=x;}
	void			set_dissolution_STI_symptom(double x){_dissolution_STI_symptom=x;}
	void			set_dissolution_ageConcPartn(double x){_dissolution_ageConcPartn=x;}
	
	void			set_dissolution_RiskGroup(vector<double>x){_dissolution_RiskGroup=x;}
	void			set_dissolution_duration(vector<double>x){_dissolution_duration=x;}
	void			set_dissolution_age(vector<double>x){_dissolution_age=x;}
	void			set_dissolution_PartnerDeficit(vector<double>x){_dissolution_PartnerDeficit=x;}

	void			set_dissolution_RiskGroup(double x, int i){_dissolution_RiskGroup[i]=x;}
	void			set_dissolution_duration(double x, int i){_dissolution_duration[i]=x;}
	void			set_dissolution_age(double x, int i){_dissolution_age[i]=x;}
	void			set_dissolution_PartnerDeficit(double x, int i){_dissolution_PartnerDeficit[i]=x;}

	void			setSpousalParameters(string spousalFile);
	void			set_spousalProgress_Age_param(double m_af, double v_af, double m_g, double v_g);
	void			set_spousalProgress_Duration_param(double k1, double k2);
	void			set_spousalProgress_DiffAgeGap_param(double meanDelta, double varDelta);
	
	void			set_dissolution_risk_param(vector<double> prm) {_dissolution_RiskGroup = prm;}
	void			set_dissolution_duration_param(vector<double> prm) {_dissolution_duration = prm;}
	void			set_dissolution_age_param(vector<double> prm) {_dissolution_age = prm;}
	
	// -- Miscellaneous
	
	
	void			setIndividual_UID(vector<unsigned long> uid);
	void			setIndividual_gender(vector<Gender> gender);
	void			setIndividual_age(vector<double> age);
	void			setIndividual_nMaxCurrSexPartner(vector<int> nMaxCurrSexPartner);
	
	void			setIndividual_nCurrSexPartner();
	
	//void			setAgeGaps(double meanAgeGap,double varAgeGap) {_formation_agegap_mean=meanAgeGap; _formation_agegap_var=varAgeGap;}
	void			setAgeSex(double ageSexMin, double ageSexMax) {_ageSexMin=ageSexMin;_ageSexMax=ageSexMax;}
	
	
	// -- Sexual activity
	
	void			set_propRiskGroup(vector<double> x) {_propRiskGroup=x;}
	void			set_propRiskGroup(double x, int i) {_propRiskGroup[i]=x;}
	
	void			set_nMaxCurrSexPartner_param_f(vector<double> x) {_nMaxCurrSexPartner_param_f=x;}
	void			set_nMaxCurrSexPartner_param_m(vector<double> x) {_nMaxCurrSexPartner_param_m=x;}
	void			set_nMaxCurrSexPartner_param_f(double x, int i) {_nMaxCurrSexPartner_param_f[i]=x;}
	void			set_nMaxCurrSexPartner_param_m(double x, int i) {_nMaxCurrSexPartner_param_m[i]=x;}
	
	void			set_sexAct_maxRate_female(double x) {_sexAct_maxRate_female=x;}
	void			set_sexAct_maxRate_male(double x) {_sexAct_maxRate_male=x;}
	
	void			set_sexActMaxRate(double maxRate_male, double maxRate_female) {_sexAct_maxRate_male=maxRate_male; _sexAct_maxRate_female=maxRate_female;}

	void			set_sexAct_reduce_age_param(vector<double> prm) {_sexAct_reduce_age_param = prm;}
	void			set_sexAct_reduce_age_param(double x, int i) {_sexAct_reduce_age_param[i] = x;}
	
	void			set_sexAct_reduce_risk_param(vector<double> prm) {_sexAct_reduce_risk_param = prm;}
	void			set_sexAct_reduce_risk_param(double x, int i) {_sexAct_reduce_risk_param[i] = x;}
	
	void			set_sexAct_reduce_STIsymptom_param(vector<double> prm) {_sexAct_reduce_STIsymptom_param = prm;}
	void			set_sexAct_reduce_STIsymptom_param(double x, int i) {_sexAct_reduce_STIsymptom_param[i] = x;}
	
	void			set_sexAct_reduce_nPartner_param(vector<double> prm) {_sexAct_reduce_nPartner_param = prm;}
	void			set_sexAct_reduce_nPartner_param(double x, int i) {_sexAct_reduce_nPartner_param[i] = x;}
	
	void			set_sexAct_proba_distribute_partnerTypes_prefSpouse(double a) {_sexAct_proba_distribute_partnerTypes_prefSpouse=a;}

	void			set_sexAct_proba_distribute_partnerTypes_sexWorkParam(vector<double> x) {_sexAct_proba_distribute_partnerTypes_sexWorkParam=x;}
	void			set_sexAct_proba_distribute_partnerTypes_sexWorkParam(double x, int i) {_sexAct_proba_distribute_partnerTypes_sexWorkParam[i]=x;}
	
	void			set_sexAct_proba_sexWorker_param(vector<double> x) {_sexAct_proba_sexWorker_param = x;}
	void			set_sexAct_proba_sexWorker_param(double x, int i) {_sexAct_proba_sexWorker_param[i] = x;}
	
	void			set_sexAct_CostSexWork_reduction(double x) {_sexAct_CostSexWork_reduction=x;}
	
	void			set_sexAct_condom_param(vector<double> x) {_sexAct_condom_param=x;}
	void			set_sexAct_condom_param(double x, int i) {_sexAct_condom_param[i]=x;}
	void			set_sexAct_TypeLowHigh(double x) {_sexAct_TypeLowHigh=x;}
	
	
	// -- STI transmission
	
	void			set_STI_SFincrease(string filename);
	
	void			set_RebHIV(vector<double> x) {_RebHIV=x;}
	void			set_RebHIV_fromFile(string filename);
	
	void			set_HIVparam(vector<double> param);
	void			set_STI_probaMaxSexTransm(STIname name, double proba);
	void			set_STI_probaSexTransm_SAT(STIname name, vector<double> param);
	void			set_STI_naturalClearanceDuration(STIname name, double x);
	void			set_STI_proba_symptomatic_female(STIname name, double x);
	void			set_STI_proba_symptomatic_male(STIname name, double x);
	void			set_STI_latentDuration(STIname name, double x);
	void			set_STI_infectiousDuration(STIname name, double x);
	void			set_STI_proba_recurrence(STIname name, double x);
	void			set_STI_recurrence_freq(STIname name, double x);
	void			set_STI_recurrence_duration(STIname name, double x);
	void			set_STI_minInfectiousness(STIname name, double x);
	
	// CSW
	
	void			set_CSW_maxRecruitmentRate(double x) {_CSW_maxRecruitmentRate=x;}
	void			set_CSW_recruitment_saturPrm_1(double x) {_CSW_recruitment_saturPrm_1=x;}
	void			set_CSW_recruitment_saturPrm_2(double x) {_CSW_recruitment_saturPrm_2=x;}
	void			set_CSW_recruitment_ageMin(double x) {_CSW_recruitment_ageMin=x;}
	void			set_CSW_recruitment_ageMax(double x) {_CSW_recruitment_ageMax=x;}	
	void			set_CSW_cessationRate(double x) {_CSW_cessationRate=x;}
	void			set_CSW_ageMax(double x) {_CSW_ageMax = x;}
	
	
	
	// === Set several parameters ===
	
	void			setAllParameters(string folder_inputs);		// Set all parameters (file names are hard coded)
	
	void			set_Population_features(string filename);
	
	void			set_Formation_parameters(string filename);
	void			set_Demographic_parameters(string filename);
	void			set_SpousalProgress_parameters(string filename);
	void			set_Dissolution_parameters(string filename);
	
	void			set_SexualActivity_parameters(string filename);
	
	void			set_CSW_demographics_parameters(string filename);
	

	
	// ===============================
	// ===    Parameters Update    ===
	// ===============================
	
	
	// Update selected parameters that are defined, through their given names
	// (convention: given name must be same as C name (e.g. "_birthRate" <--> _birthRate))
	void		UpdateSelectedParameter(string paramName, double value);
	void		UpdateSelectedParameter(vector<string> paramName, vector<double> value);
	void		UpdateSelectedParameter_file(string filename, unsigned int paramSet = 1);
	
	
	
	
	// ===============================
	// ===       Demographics      ===
	// ===============================

	void			addIndividual(Individual I);
	
	void			youthArrivals(double prd,bool save_trace_file);		// Temporary function: introduce young in the population
	
	void			CSWrecruitment(double prd);		// Recruitment of new CSW
	void			CSWcessation(double prd);		// CSW droping out commercial sex
	
	void			kill(unsigned long uid, bool save_trace_file);		// Kills individual

	void			save_death_hazard(string filename); // For external check
	double			death_hazard(unsigned long i);	// Annual death probability of ith individual
	void			deathEvents(double prd,
								bool save_trace_file);	// Select individuals that die withn the next period "prd" given an age,[[disease stage]]
	
	unsigned long	getMaxUID();
	
	
	unsigned long	census_alive();			// Number of individuals alive
	unsigned long	census_alive(Gender);	// Number of individuals alive for a given gender
	unsigned long	census_dead();			// Number of dead individuals
	
	unsigned long	census_Females();		// Number of females alived in the whole population
	
	unsigned long	census_Partnered(int numberOfPartners); // Number of individuals having "numberOfPartners" partners
	unsigned long	census_singles(); // Number of single individuals (no sex partners other than CSW)
	unsigned long	census_singles(Gender g); // Number of single individuals (no sex partners other than CSW)
	double			census_ratioSingles();	// Ratio of single individual to the total population
	double			census_ratioSingles(Gender g);	// Ratio of single individual to the total population
	
	unsigned long	census_widow();			// Number of individuals who are widowed
	unsigned long	census_divorced();		// Number of individuals who are divorced
	
	unsigned long	census_circum();		// Number of male individuals who are circumcised
	
	vector<double>	census_getAgeGap(unsigned long uid); // Get age gaps of individual with all his/her partners
	vector<double>	census_AgeGap();		// Get age gaps of all individuals with all his/her partners 
											// WARNING: WILL APPEAR TWICE (one for each partner)
	
	vector<double>	census_ageGapDistribution(vector<double> bins);	// Distribution of age gap for the whole population
	vector<double>	census_ageDistribution(vector<double> ageBreaks);	// Distribution of ages in the whole population
	vector<double>	census_ageFirstSexDistribution(vector<double> ageBreaks, Gender thegender);	// Distribution of ages at 1st sex for the whole population
	vector<double>	census_ageGapFirstSexSpouseDistribution(vector<double> ageBreaks, Gender thegender);	// Distribution of ages at 1st sex for the whole population
	
	vector<double>	census_ageMalesVisitCSWDistribution(vector<double> ageBreaks);	// Distribution ages males ever visited CSW
	vector<double>	census_ageMalesVisitCSWDistribution(vector<double> ageBreaks,
														double maxDurationSinceLastVisit);
	
	vector<unsigned long>	census_nMaxCurrSexPartnerCount(vector<unsigned long> nBreaks, Gender g);	// Distribution of maximum number of concurrent sex partners
	vector<double>			census_nMaxCurrSexPartnerDistrib(vector<double> nBreaks, Gender g);	// Distribution of maximum number of concurrent sex partners
	vector<double>			census_nLifeSexPrtnrDistrib(vector<double> nBreaks, Gender g);	// Distribution of lifetime number of sex partners
	
	unsigned long			census_pregnant();		// total number of pregnant females
	unsigned long			census_pregnant(int riskgroup);		// total number of pregnant females of given risk group
	vector<unsigned long>	census_pregnant_UID();	// UID of all pregnant women

	
	//DEBUG ----
	void census_nLifeSexPrtnrDistrib_todelete(Gender g);
	// ----
	
	vector<unsigned long>	census_CSW();			// Return all UIDs of commercial sex workers (CSW)
	double					census_maleVisitCSW();	// Return proportion of maleswho have ever visited a CSW
	double					census_maleVisitCSW(double maxDurationSinceLastVisit);
	
	vector<unsigned long>	census_STIinfected();	// i^th element is the number of individuals infected with STI #i
	unsigned long			census_STIcoinfected(STIname s1, STIname s2) ; // Number of individuals co-infected with s1 and s2
	unsigned long			census_STIcoinfected(STIname s1, STIname s2, int riskgroup) ;

	vector<unsigned long>	census_riskGroup();		// i^th elementis number of individual in risk group #i. Last element CSW
	
	
	// ===============================
	// ===    Sexual Activity      ===
	// ===============================
	
	void			reset_sexActs();			// Resets the number of sex acts and partners' UID for this period
	void			update_sexActs(double period, bool save_trace_file);
	
	void			sexActs_number_males(unsigned long uid, double period, bool save_trace_file);
	void			sexAct_distribute_partnerTypes(unsigned long uid, bool save_trace_file);
	void			sexAct_distribute_individualPartners(unsigned long uid, bool save_trace_file);
	void			sexAct_distribute_ActTypes(unsigned long uid);

	unsigned long	count_sexActs(int riskgroup);
	
	vector<unsigned long>	pregnantPotentialFemales(); // count number of females who can become pregnant (must had sex)
	
	double			proba_nMaxCurrSexPartner(Gender g, int riskgroup);
	double			probaSex_sexworker(int riskgroup);
	double			probaSex_type0(int r1, int r2);
	double			probaSex_type1(int r1, int r2);
	
	// === Analytics ===
	
	double			averageAge() {return averageElements(get_age());}
	
	
	
	// === Partnerships ===
	
	dcMatrix		getPartnershipsUID();		// Retrieves the pair of UIDs for each partnerships
	dcMatrix		getSpousesUID();			// Retrieves the pair of UIDs for each spousal partnerships
	
	double			getPartnershipDuration(unsigned long uid1, unsigned long uid2);
	
	bool			alreadyPartners(unsigned long uid1, unsigned long uid2); // 'true' if already in partnership, false else
	
	void			formPartnership(unsigned long uid1, unsigned long uid2);
	void			dissolvePartnership(unsigned long uid1, unsigned long uid2,
										bool save_trace_file);
	
	double			meanAgePartners(unsigned long uid);		// Calculate mean age of all partners of a given individual
	
	unsigned long	findPositionIn_partnershipsMatrix(unsigned long uid_fem, unsigned long uid_male); // find the position (row) of the partnership (uidf,uidm) in _partnershipsMatrix
	
	unsigned int	getMaxDegree();		// Retrieve maximum degree of network (=max number of partners)

	void			update_UID_PartnershipMatrix(Individual I);	// Update "_PartnerMatrix_(fe)maleUID" when individuals are added to the population
		
	// -- Match functions for partnership formation
	
	double			form_age(double age_f);						// Impact of age on formation probability
	double			form_age_male(double age_m);				// Impact of male's age on formation probability
	double			form_ageGap(double ageGap);					// Impact of age on formation probability
	double			form_riskGroup(int risk_f, int risk_m);		// Impact of risk group on formation probability
	double			form_deficit(double def_f, double def_m);	// Impact of Spartnership deficit on formation probability
	double			form_STIsymptom(unsigned long uid_f, unsigned long uid_m);// Impact of STI symptom on formation probability
	
	bool			isFormationPossible(unsigned long uid1, unsigned long uid2, bool save_to_file);
	
	
	// -- Spousal progression
	
	bool			spousalProgression(unsigned long uid1, unsigned long uid2);
	void			spousalScanAllCasual();
	
	
	// -- Dissolution functions of partnerships
	
	bool			isDissolutionPossible(unsigned long uid1, unsigned long uid2); // True is candidate partnership will dissolve
	
	double			dissolve_spouse(unsigned long uid1, unsigned long uid2);
	double			dissolve_riskGroup(unsigned long uid1, unsigned long uid2);
	double			dissolve_duration(unsigned long uid1, unsigned long uid2);
	double			dissolve_age(unsigned long uid1, unsigned long uid2);
	double			dissolve_partnerDeficit(unsigned long uid1, unsigned long uid2);
	
	double			dissolve_STIsymptoms(unsigned long uid1, unsigned long uid2);
	
	double			dissolve_riskGroup_fct(int r1, int r2);
	double			dissolve_duration_fct(double prtnr_duration);
	double			dissolve_age_fct(double age1, double age2);
	double			dissolve_partnerDeficit_fct(double d1, double d2);

	
	// --- counting
	
	unsigned long	totalNumberPartnerships();
	unsigned long	nUnmatchedPartnerships();	// Number of unmatched potential partnerships in whole population
	unsigned long	totalNumberSpousalPartneships() {return _totalNumberSpousalPartnerships;}

	
		// --- Formation & Dissolution processes
	
	void			formOnePartnershipFromFemale_rand(unsigned long uid_fem, bool save_trace_file);// Attempt to forms 1 partnership with specified female and a random available male

	vector<unsigned long>	formUIDcandidateFemale(double prd); // Returns UIDs of all females candidate for new partnership
	
	void			formPartnerships(double prd,bool save_trace_file);
	void			dissolvePartnerships(double prd, bool save_trace_file);	
	
	
		// --- CSW
	
	unsigned long	chooseRandomCSW();	// Randomly pick a CSW and returns her UID
	bool			isCSW(unsigned long uid) {return (_individual[uid].get_riskGroup()==_CSWriskGroup?true:false);}
	bool			at_least_one_CSW();
	
	// === Sexual Activity ===
	
	double			sexAct_reduce_age(double age);				// Reduce sexual activity based on individual's age
	double			sexAct_reduce_riskGroup(int r);				// Reduce sexual activity based on individual's risk group
	double			sexAct_reduce_STIsymptom(bool s,Gender g);	// Reduce sexual activity based on individual's presence of STI symptoms
	double			sexAct_reduce_AIDS(double aidsDuration);	// Reduce sexual activity based on individual's AIDS progression
	double			sexAct_reduce_nPartners(int nPartners);		// Reduce sexual activity based on individual's number of current partners
	
	void			add_sexAct(unsigned long uid1, unsigned long uid2,
							   int nSexActs,
							   bool save_trace_file);	// Add sex acts b/w 2 individuals
	
	
	// ===============================
	// ===    STI TRANSMISSIONS    ===
	// ===============================
	
	void			STI_setAllParameters(string filename);
    void            STI_set_genomes(string filename);
	void			STI_set_initial_prevalence(string filename); // Set inital prevalence for the Population
	
	
	vector<double>	STI_CalcProbaTransmission(unsigned long uid_infect, 
												   unsigned long uid_suscep);
	
	void			STI_acquire(STIname stiname, unsigned long uid,
								double timeStep);	// "uid" acquires STI "stiname" and sets the STI duration to "prd"
	
	void			STI_transmission_indiv(STIname stiname, 
										   unsigned long uid1, unsigned long uid2,
										   double timeStep, bool save_trace_file);
	
	vector<unsigned long> STI_transmissions(double timeStep, bool save_trace_file);	// creates all transmission events following sex acts
	
	double			STI_coInfection_oddsRatio(unsigned long uid, int sti); // odds ratio for a given individual, to STI 'sti'

	double			STI_proba_MTCT(STIname sti, double stiduration);		// probability of mother-to-child transmission
	void			set_STI_MTCT(unsigned long uid, vector<bool> mtct);
	void			set_STI_MTCT(unsigned long uid, STIname stiname, bool mtct);
	
	// Terminate STIs that have reached the max duration
	// before immune system clearance, accros the whole population
	
	void			STI_update_naturalClearance();
	
	
	// Update Population-level STIs features to the all individuals
	void			STI_update_templateToAllIndividuals();
	
	
	int				STI_find_index_position(STIname);
	void			STI_save_all_infectivityCurves(string filerootname);
	
	vector< vector<double> > get_infectivityCurve(STIname sti, Gender g);
    
    // Discordant partnership
    
    bool            STI_isDiscordPartner(unsigned long uid1, unsigned long uid2); // is the partnership discordant for at least one STI
    bool            STI_atLeastOneDiscordPartner(unsigned long uid); // is there at least one partner who is discordant
	
	// === STI Prevalences/Incidences ===
	
	vector<double>	STI_prevalences();							// Global prevalence for every STIs
	vector<double>	STI_prevalences(vector<STIname> stinames);	// Global prevalence for selected STIs
	double			STI_prevalence(STIname s);					// Global prevalence for just one STI
	double			STI_prevalence(STIname s, int riskGroup);	// Prevalence for just one STI, within a given risk group

	dcMatrix		STI_prevalence_by_riskGroup(vector<STIname> stinames);							// Prevalence per risk group, for every STIs
	vector<double>	STI_prevalence_by_age(STIname s, vector<double> agebreaks);						// Prevalence by age for just one STI
	vector<double>	STI_prevalence_by_age(STIname s, Gender g, vector<double> agebreaks);			// Prevalence by age for just one STI, within a given gender
	vector<double>	STI_prevalence_by_age(STIname s, int riskGroup,vector<double> agebreaks);		// Prevalence by age for just one STI, within a given risk group
	
	
	void			update_STI_mtct_cumcount(vector<bool> mtct);

	
	// --- Reproductive numbers ---
	
	vector<unsigned long>	Reff_inst_all_infectious(STIname stiname);	// Instantaneous effective reproductive number for all infectious individuals for a given STI
	double					Reff_inst_mean(STIname stiname);			// Instantaneous Mean effective reproductive number for a STI
	void					update_secondary_cases(STIname stiname);	// Cumulative effective reproductive number for all infectious individuals for a given STI
	double					Reff_cum_mean(STIname stiname);				// Cumulative Mean effective reproductive number for a STI
    
    // ===============================
	// ===       TREATMENTS        ===
	// ===============================
	
	void			treat_indiv(unsigned long uid, STIname sti);
	void			cure_indiv(unsigned long uid, STIname sti);
	
	
	// ===============================
	// ===       VACCINATION       ===
	// ===============================

	void		vaccinate_indiv(unsigned long uid, STIname stiname);
	void		set_STI_immunity(unsigned long uid, STIname stiname, double immunity);
	void		set_STIsusceptFactor(unsigned long uid, STIname stiname, double x);
    
    
    // ===============================
    // ===       GENETICS     ===
    // ===============================

    void    update_genomes_mutations(double timestep);
    void    update_genomes_clear();
    
    vector<string> get_population_genomes(int sti_pos);
    
    
	// ===============================
	// ===       MISCELLENAOUS     ===
	// ===============================

	
	// === Time, Ages, Durations updates ===
	
	void	updateAllDurations(double timeStep);
	void	updateAllAges(double timestep);
	
	// === Helpers ===	
	
	void			displayInfo(bool IndividualDetails);
	void			displayInfo_SexActs_lastPeriod();
	void			displayInfo_SexActs_Pairs_lastPeriod();
	void			displayInfo_STI_SFincrease();
	
	// retrieve all UIDs:
	vector<unsigned long> getUIDs();
	
	unsigned long	getUIDfemale(unsigned long uid1, unsigned long uid2); // retrieves the UID female of a given pair
	unsigned long	getUIDmale(unsigned long uid1, unsigned long uid2); // retrieves the UID male of a given pair
	
	dcDataFrame		export_to_dataframe();
	
	// === Output Files
	
	void			saveToCSVFile(string pathFile);

	void			save_outputs_demog(string pathFolder);		// Save to file demographic outputs (e.g. age distribution)
	void			save_outputs_demog(string pathFolder,
									   vector<double> ageBreaks,
									   unsigned int iMC=0,
									   unsigned int idate=0);

	void			save_outputs_prtnr(string pathFolder);		// Save to file partnership outputs (e.g. age gap distribution, single ratio,...)
	void			save_outputs_prtnr(string pathFolder,
									   vector<double> ageGapBreaks,
									   unsigned int iMC=0,
									   unsigned int idate=0);
	
	void			save_outputs_sex(string pathFolder);		// Save to file sexual behavior outputs (e.g. number of lifetime sex partners distribution)
	void			save_outputs_sex(string pathFolder,
									 vector<double> ageBreaks_f,
								  vector<double> ageBreaks_m,
								  vector<double> ageGapBreaks_f,
								  vector<double> ageGapBreaks_m,
								  vector<double> nBreaks_f,
								  vector<double> nBreaks_m,
								  vector<double> ageBreaks_malesVisitCSW,
									 unsigned int iMC=0,
									 unsigned int idate=0);

	
	void			save_outputs_epi(string pathFolder,
									 vector<STIname> STIname_targeted,
									 vector<double> agebreaks_HIVprev_f,
									 vector<double> agebreaks_HIVprev_m,
									 unsigned int iMC=0,
									 unsigned int idate=0);		// Save to file epidemiological outputs (e.g. STI prevalences)

	
    void			FileOutput(string pathFile);
	void			debug_check_partnerUID();
};




vector<unsigned long> convertTo_unsignedLong(vector<double> x, bool deleteFirst);
vector<int> convertTo_int(vector<double> x, bool deleteFirst);
vector<Gender> convertTo_Gender(vector<double> x, bool deleteFirst);

double		expGauss(double x, double mean, double var);
double		expGauss2(double x, double k1, double k2);

double		probaSex_sexworker(int riskgroup,int maxriskgroup, double c1, double c2);
double		weibull_hazard(double t, double shape,double scale);

#endif