/*
 *  individuals.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-07-20.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef individuals_h
#define individuals_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <cstdlib>

#include "dcTools.h"
#include "STI.h"
#include "globalVar.h"
#include "gender.h"

using namespace std;


/** Class representing the individual agents of the population
 Each individual had a unique identification number (UID)
 All features relevant to a given individual are stored
 in this class. For example:
 - gender
 - age
 - number of current partners
 - STI infection status
 - treatment of a given STI
 
 */

class Individual
{
	
	// == Fundamental features ==
	
	bool			_isAlive;
	unsigned long	_UID;		// Unique ID: this is the position
	// of this Individual in the "Population"
	Gender			_gender;
	double			_age;
	int				_riskGroup;
	
	bool			_isWidow;	// Tracks if the individual is widowed
	bool			_isDivorced;// Tracks if the individual is divorced
	
	bool            _isCircum;  // Circumcision status (male only!)
	
	bool			_isPregnant; // Pregnancy status (for female only!)
	double			_gestationDuration;
	int				_nChildBorn;			// number of child born (for female only)
	double			_dateInPopulation;		// first date when entered the Population
	
	// == Partnerships ==
	
	int						_nCurrSexPartner;		// number of current sex partners
	int						_nMaxCurrSexPartner;	// max number of current sex partners
	int						_nCurrSpouse;			// number of current spouses
	
	vector<unsigned long>	_partnerUID;			// UID of all partners (including spouses)
	vector<bool>			_isSpousalPartner;		// _isSPousalPartner[i]=1 if ith partner a spouse, 0 if casual
	
	int						_nLifetimePartner;		// lifetime number of partners
	int						_nLifetimeSpouse;		// lifetime number of spouses
	vector<double>			_partnershipDuration;	// Duration of each current partnerships
	double					_singleDuration;		// Duration of current single status
	
	double					_ageFirstPartner;		// Age at first partnership formation (must be equal to _ageFirstSex; internal check)
	double					_ageFirstSpouse;		// Age at first spousal partnership formation
	
	
	// == Sexual activity ==
	
	int				_nSexActs_period;			// Number of sex acts for a given period (with all partners)
	int				_nSexActs_spouse_period;	// Number of sex acts for a given period with all spouses only (for males only)
	int				_nSexActs_casual_period;	// Number of sex acts for a given period with all casual only (for males only)
	int				_nSexActs_sexworker_period;	// Number of sex acts for a given period with sex workers only (for males only)
	
	unsigned long	_nSexActs_lifetime;			// Number of sex act during lifetime.
	
	double			_ageFirstSex;				// Age at first sex
	
	
	vector<unsigned long>	_UID_sexAct_period;		// UIDs of individual having sex with this individual during a given period
	vector<int>				_UID_n_sexAct_period;	// number of sex acts with each UID WHERE SEX ACT PERFORMED during the period
	vector<int>				_UID_n_sexAct_Type0_period;	// number of sex acts of type 0 with each partner within the period
	vector<int>				_UID_n_sexAct_Type1_period;	// number of sex acts of type 1 with each partner within the period
	vector<int>				_UID_n_sexAct_Type2_period;	// number of sex acts of type 2 with each partner within the period
	
	// CSW are not considered as formal partners
	// so the count of sex act types is segregated
	// (not in vector because assumed same CSW for a given period)
	int						_CSW_n_sexAct_Type0_period;	// number of sex acts of type 0 with the CSW during the period
	int						_CSW_n_sexAct_Type1_period;
	int						_CSW_n_sexAct_Type2_period;
	
	bool					_ever_visited_CSW;			// If individual (male) has ever visited a CSW
	double					_lastVisitCSWDuration;		// Duration since last visit to CSW (for males only)
	
	
	// == Diseases ==
	
	vector<STI>		_STI;					// STI modelled - Copy of template from Population
	
	vector<double>	_RebHIV;				// Rebound of HIV infectivity - Copy of template from Population
	
	vector<double>	_STIsusceptFactor;		// susceptibility factors to every STIs modelled for this individual (comes after the population level susceptibility defined in matrix _STI_SFincrease)
	vector<double>	_STIduration;			// Duration since STI infection. Not infected with ith STI: _STIduration[i]=0
	vector<bool>	_STIsymptom;			// ith element=1 if ith STI is symptomatic (0 else)
	
	vector< vector<unsigned long> > _STI_secondary_cases; // UID of all secondary cases infected by this individual, for every STI
	vector<bool>	_STI_MTCT;				// whether mother to child transmission occurs (relevant only for pregnant females)
	
	// == Treatment ==
	
	vector<double>	_STItreatDuration;		// duration since treatment started
	vector<int>		_STItreatTMS;			// Microbiological success of treatment in this individual
	vector<double>	_STItreatAdherence;		// Adherence to STI treatments
	
	
	// == Vaccination ==
	
	vector<double>	_STI_vacc_time;		// time when vaccinated
	vector<bool>	_STI_vacc;			// STI vaccination status
	vector<double>	_STI_immunity;		// STI immunity status: 1=full immunity, 0=no immunity

	
public:
	
	// === Constructors ===
	
	Individual(){}
	
	// core construction
	void	construct(unsigned long uid, Gender g,
					  double age, int maxSexPartner,
					  int riskGroup, int nLifetimePartner,
					  bool isWidow, bool isDivorced,
					  bool isCircum);
	
	
	
	Individual (unsigned long uid, Gender g, double age,
				int maxSexPartner,
				int riskGroup, int nLifetimePartner,
				bool isWidow, bool isDivorced, bool isCircum,
				vector<double> STIdurations,
				vector<bool> STIsymptoms,
				vector<STI> templateSTI,
				vector<double> RebHIV);
	
	
	
	// === Get Functions ===
	
	bool						isAlive() {return _isAlive;}
	unsigned long				get_UID() {return _UID;}
	double						get_age() {return _age;}
	Gender						get_gender() {return _gender;}
	int							get_riskGroup() {return _riskGroup;}
	
	bool						get_isPregnant() {return _isPregnant;}
	double						get_gestationDuration() {return _gestationDuration;}
	int							get_nChildBorn() {return _nChildBorn;}
	double						get_dateInPopulation() {return _dateInPopulation;}
	
	
	int							get_nCurrSexPartner() {return _nCurrSexPartner;}
	int							get_nMaxCurrSexPartner() {return _nMaxCurrSexPartner;}
	int							get_nCurrSpouse() {return _nCurrSpouse;}
	int							get_nCurrCasual() {return _nCurrSexPartner-_nCurrSpouse;}
	
	int							get_nLifetimePartner() {return _nLifetimePartner;}
	int							get_nLifetimeSpouse() {return _nLifetimeSpouse;}
	
	vector<unsigned long>		getPartnerUID();				// Retrieve vector of all partners' UID
	unsigned long				getPartnerUID(int nPartner);	// Retrieve UID on nth partner
	
	vector<unsigned long>		getSpouseUID();					// Retrieve vector of all spouses' UID
	vector<unsigned long>		getCasualUID();					// Retrieve vector of all casuals' UID
	
	int							getPartnershipPosition(unsigned long uid_partner);
	vector<bool>				get_isSpousalPartner() {return _isSpousalPartner;}
	
	double						get_singleDuration() {return _singleDuration;}
	
	double						get_ageFirstPartner() {return _ageFirstPartner;}
	double						get_ageFirstSpouse() {return _ageFirstSpouse;}
	
	vector<double>				getPartnershipDuration() {return _partnershipDuration;}
	
	bool						get_isWidow() {return _isWidow;}
	bool						get_isDivorced() {return _isDivorced;}
	
	bool						get_isCircum() {return _isCircum;}
	
	int							get_nSexActs_period() {return _nSexActs_period;}
	int							get_nSexActs_spouse_period() {return _nSexActs_spouse_period;}
	int							get_nSexActs_casual_period() {return _nSexActs_casual_period;}
	int							get_nSexActs_sexworker_period() {return _nSexActs_sexworker_period;}
	unsigned long				get_nSexActs_lifetime() {return _nSexActs_lifetime;}
	
	bool						get_ever_visited_CSW() {return _ever_visited_CSW;}
	double						get_lastVisitCSWDuration() {return _lastVisitCSWDuration;}
	
	double						get_ageFirstSex() {return _ageFirstSex;}
	
	vector<unsigned long>		get_UID_sexAct_period(){return _UID_sexAct_period;}
	vector<int>					get_UID_n_sexAct_period(){return _UID_n_sexAct_period;}
	
	// Retrieve number of sex acts by sex act type for a given partner (uid_partner)
	// and (implicitly) during the last period (if simulation run)
	int							get_UID_n_sexAct_Type0_period(int uid_partner) {return _UID_n_sexAct_Type0_period[uid_partner];}
	int							get_UID_n_sexAct_Type1_period(int uid_partner) {return _UID_n_sexAct_Type1_period[uid_partner];}
	int							get_UID_n_sexAct_Type2_period(int uid_partner) {return _UID_n_sexAct_Type2_period[uid_partner];}
	
	// -- STI --
	
	vector<STI>			get_STI() {return _STI;}
	
	vector<double>		get_STIduration() {return _STIduration;}
	double				get_STIduration(STIname stiname);
	
	vector<double>		get_STIsusceptFactor() {return _STIsusceptFactor;}
	double				get_STIsusceptFactor(STIname stiname);
	
	vector<bool>		get_STIsymptom() {return _STIsymptom;}
	bool				get_STIsymptom(STIname);
	
	vector<double>		get_STItreatDuration() {return _STItreatDuration;}
	double				get_STItreatDuration(STIname sti);
	int					get_STItreatTMS(STIname sti);
	double				get_STItreatAdherence(STIname sti);
	
	
	vector<bool>		get_STI_vacc(){return _STI_vacc;}
	vector<double>		get_STI_vacc_time(){return _STI_vacc_time;}
	
	vector<double>		get_STI_immunity() {return _STI_immunity;}
	double				get_STI_immunity(STIname stiname);

	vector<bool>		get_STI_MTCT() {return _STI_MTCT;}

	
	vector<unsigned long> get_STI_secondary_cases(STIname stiname);
    
    // -- Genetics --
    string      get_genome(int i) {return _STI[i].get_genome();}
    string      get_genome(STIname stiname) {return _STI[positionSTIinVector(stiname, _STI)].get_genome();}
    
	
	// ==============================
	// ====== Set Functions ======
	// ==============================
	
	void		set_UID(unsigned long uid) {_UID = uid;}
	void		kill() {_isAlive=0;}
	
	void		set_age(double age) {_age = age;}
	void		increase_age(double timeStep) {_age += timeStep;}
	
	void		set_gender(Gender g) {_gender = g;}
	
	void		set_dateInPopulation(double x) {_dateInPopulation = x;}
	
	void		set_partnerUID(vector<unsigned long> pUID) {_partnerUID = pUID;}
	
	void		set_isWidow(bool x) {_isWidow = x;}
	void		set_isDivorced(bool x) {_isDivorced = x;}
	
	void		set_isPregnant(bool x) {_isPregnant = x;}
	void		set_gestationDuration(double x) {_gestationDuration=x;}
	
	void		set_nChildBorn(int n) {_nChildBorn = n;}
	void		increment_nChildBorn() {_nChildBorn++;}
	
	void		set_nMaxCurrSexPartner(int n) {_nMaxCurrSexPartner = n;}
	void		set_nCurrSexPartner(int n);
	void		decrease_nCurrSexPartner(); //_nCurrSexPartner--;
	
	void		set_nLifetimePartner(int n) {_nLifetimePartner = n;}
	void		increment_nLifetimePartner() {_nLifetimePartner++;}
	
	void		initPartnershipDuration(); // Initialize _partnershiDuration when population newly created
	void		setPartnershipDuration(int n, double d); // sets duration=d to partnership #n
	void		increasePartnershipDuration(int n, double timeStep); // increase duration by timeStep of partnership #n
	
	void		set_singleDuration(double d) {_singleDuration = d;}
	void		increaseSingleDuration(double timestep) {_singleDuration += timestep;}
	
	void		setIsSpousalPartner(int p, bool x);
	void		increase_nCurrSpouse() {_nCurrSpouse++; _nLifetimeSpouse++;}
	void		decrease_nCurrSpouse() {_nCurrSpouse--;}
	
	void		set_nSexActs_period(int n) {_nSexActs_period=n;}
	void		set_nSexActs_spouse_period(int n) {_nSexActs_spouse_period=n;}
	void		set_nSexActs_casual_period(int n) {_nSexActs_casual_period=n;}
	void		set_nSexActs_sexworker_period(int n) {_nSexActs_sexworker_period=n;}
	
	void		set_ever_visited_CSW(bool x) {_ever_visited_CSW=x;}
	void		set_lastVisitCSWDuration(double x) {_lastVisitCSWDuration=x;}
	
	void		increase_lastVisitCSWDuration(double x) {if(_lastVisitCSWDuration>0) _lastVisitCSWDuration+=x;}
	
	void		set_ageFirstSex(double x) {_ageFirstSex=x;}
	void		set_ageFirstPartner(double x) {_ageFirstPartner=x;}
	void		set_ageFirstSpouse(double x) {_ageFirstSpouse=x;}
	
	void		increase_nSexActs_period(int n) {_nSexActs_period+=n;}
	void		increase_nSexActs_spouse_period(int n) {_nSexActs_spouse_period+=n;}
	void		increase_nSexActs_casual_period(int n) {_nSexActs_casual_period+=n;}
	void		increase_nSexActs_lifetime(int n) {_nSexActs_lifetime += n;}
	
	void		add_sexAct(unsigned long uid_partner, int n);
	
	void		set_UID_n_sexAct_Type_period(vector<unsigned int> nSexType);
	
	void		set_riskGroup(int r) {_riskGroup=r;}
	
	void		init_STIsusceptFactor();
	void		set_STIsusceptFactor(int i, double x) {_STIsusceptFactor[i]=x;}
	
	void		set_STI(vector<STI> s) {_STI = s; _STIduration.resize(s.size(),0.0);}
	
	void		set_STIduration(int sti_index, double duration);
	void		set_STIduration(STIname	stiname, double duration);
	
	void		set_STIsymptom(int sti_index, bool isSymptomatic);
	
	void		set_STI_MTCT(vector<bool> mtct) {_STI_MTCT = mtct;}
	void		set_STI_MTCT(STIname stiname, bool mtct);
	
	void		set_RebHIV(vector<double> x) {_RebHIV=x;}
	
	
	
	void		set_STItreatDuration(vector<double> x){_STItreatDuration=x;}
	void		set_STItreatTMS(vector<int> x) {_STItreatTMS=x;}
	void		set_STItreatAdherence(vector<double> x){_STItreatAdherence=x;}
	
	void		set_STItreatDuration(double x, STIname sti);
	void		set_STItreatTMS(int x, STIname sti) ;
	void		set_STItreatAdherence(double x, STIname sti);
	
	void		set_STI_vacc_time(vector<double> x) {_STI_vacc_time = x;}
	void		set_STI_vacc_time(int i,double x) {_STI_vacc_time[i] = x;}
	
	void		set_STI_vacc(vector<bool> x) {_STI_vacc=x;}
	void		set_STI_vacc(int i, bool x) {_STI_vacc[i]=x;}
	
	void		set_STI_immunity(vector<double> x) {_STI_immunity=x;}
	void		set_STI_immunity(int i, double x) {_STI_immunity[i]=x;}

    void        set_genome(string genome, int i) {_STI[i].set_genome(genome);}
    void        set_genome(string genome, STIname stiname) { _STI[positionSTIinVector(stiname, _STI)].set_genome(genome);}
	
	// === Partnerships ===
	
	void		addPartner(unsigned long uid);	// Main function to add a partnership
	void		erasePartner(unsigned long uid);
	
	void		convertToSpouse(unsigned long uid_partner);
	
	bool		isOpenToNewPartnership();
	
	
	bool		isSpouse(unsigned long uid);		// 1 if "uid" is a spouse of this individual
	bool		isSingle() {return (_nCurrSexPartner==0);}
	
	
	
	// === Sexual activity ===
	
	double		sexAct_reduce_age();	// Reduce the rate of sexual activity according to age
	
	void		reset_UID_sexAct_period() {_UID_sexAct_period.clear();}
	void		reset_UID_n_sexAct_Type(); // Reset all vectors _UID_n_sexAct_TypeX_period
	void		reset_UID_n_sexAct_period();
	
	int			nSexActType1or2();	// Total number (across all partners) of sex acts of type 1 or 2 (not type0[:condom])
	
	// === STI ===
	
	void			STI_initializeAll(vector<STI> templateSTI);		// Initialize all STI when this individual is "constructed"
	
	void			STI_resetAllDurations();						// Reset all STI durations to 0
	
	void			set_STIsymptom(STIname stiname, bool isSymptomatic);
	
	bool			is_symptomatic();	// If this individual has symptoms from any STIs
	bool			is_symptomatic(STIname stiname);	// If this individual has symptoms from that given STIs
	
	void			STI_acquireInfection(STIname name, double durationSinceInfection);	// STI infection is acquired by this individual
	
	bool			STI_anyInfection();		// Check if there is any STI infection
	bool			STI_infected(STIname stiname);	// If this individual is infected with a given STI
	
	
	int				STI_nModelled();		// Number of STIs modelled
	int				STI_find_index_position(STIname s);
	
	vector<string>	STI_listInfection();    // Returns all STI infections of this individual
	
	vector<double>	STI_IC();               // Infectivity curve for every STIs modelled
	void			display_STI_IC();		// Display values of all infectivity curves at the current simulation time
	
	void			STI_updateSusceptFactor();	// update the susceptibility factors, given this individual's features
	
	void			increaseSTIdurations(double prd);
	
	void			add_STI_secondary_cases(STIname stiname, unsigned long uid_infected);
	void			clear_STI_secondary_cases(STIname stiname);
	
	
	// === STI Treatment ===
	
	void		treat(STIname);     // Attempts to treat an infecting STI
	bool		STI_treated(STIname stiname);	// If this individual is treated against a given STI
	
	void		increase_STItreatDuration(double prd);
	
	double		TRE(STIname);		// Treatment effect on infectivity curve
	
	
	
	// === STI vaccination ===
	
	double		VRE(STIname);		// Vaccination effect on infectivity curve
	
	
	// === Miscellenaous ===
	
	
	
	
	// === Helper Functions ===
	
	void	displayInfo();
	
};


#endif