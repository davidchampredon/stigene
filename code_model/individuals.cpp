/*
 *  individuals.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-07-20.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "individuals.h"


/* ************* CONSTRUCTORS ************* */


void Individual::construct(unsigned long uid, Gender g, double age,
						   int maxSexPartner, int riskGroup, int nLifetimePartner,
						   bool isWidow, bool isDivorced,
						   bool isCircum)
{
	_isAlive = 1;			// By default, a new Individual is alive (!)
	
	// By default, no sex partner (need to be assigned a partner from specific process)
	_nCurrSexPartner	= 0;
	_nCurrSpouse		= 0;
	_partnerUID.clear();
	
	_nMaxCurrSexPartner	= maxSexPartner;
	stopif(maxSexPartner==0,"_nMaxCurrSexPartner must be greater than 1 for any individual");

	
	// TO DO: set lifetime partner to 0
	_nLifetimePartner	= nLifetimePartner;
	
	_nLifetimeSpouse	= 0; // TO DO: generate this like nLifetimePartner
	_singleDuration		= 0; // TO DO: generate this
	
	
	_ageFirstSex		= -99.0;		// By default, no sex ever occured
	_ageFirstSpouse		= -99.0;		// By default, no spouses ever
	_ageFirstPartner	= -99.0;		// By default, no partners ever
	
	_UID		= uid;
	_gender		= g;
	_age		= age;
	
	_riskGroup	= riskGroup;
	
	_isWidow	= isWidow;
	_isDivorced	= isDivorced;
	
	_nSexActs_period		= 0;	// Number of sex acts for a given period (with all partners)
	_nSexActs_spouse_period = 0;	// Number of sex acts for a given period with all spouses only (for males only)
	_nSexActs_casual_period = 0;	// Number of sex acts for a given period with all casual only (for males only)
	_nSexActs_lifetime		= 0;

	_nSexActs_sexworker_period	= 0;
	_ever_visited_CSW			= false;
	_lastVisitCSWDuration		= 0.0;
	
	_isCircum = isCircum;	// Circumcision status (relevant for males only)
	
	_isPregnant = false;
	_gestationDuration = 0.0;
	_nChildBorn = 0;
	
}



Individual::Individual(unsigned long uid, Gender g, double age, int maxSexPartner,
					   int riskGroup, int nLifetimePartner,
					   bool isWidow, bool isDivorced, bool isCircum,
					   vector<double> STIdurations,
					   vector<bool> STIsymptoms,
					   vector<STI> templateSTI,
					   vector<double> templateRebHIV)
{
	construct(uid, g, age, maxSexPartner, riskGroup,
			  nLifetimePartner, isWidow, isDivorced, isCircum);
	
	// === STIs ===
	
	STI_initializeAll(templateSTI);
	
	// Initialize durations of infections
	int nSTI = templateSTI.size() ;
	
	
	string errmsg = "STI initialization not coherent!";
	stopif(nSTI != STIdurations.size(),errmsg);
	
	for (int i=0; i<nSTI; i++){
		set_STIduration(i,STIdurations[i]);
		set_STIsymptom(i, STIsymptoms[i]);
		_STI_MTCT.push_back(false);
	}
	set_RebHIV(templateRebHIV);
}





/* ************* PARTNERSHIPS ************* */

void Individual::addPartner(unsigned long uid)
{
	/// Main function to add a partnership
	
	// If this is the first partner ever,
	// then record Individual's age
	
	if (_ageFirstPartner<0){
		_ageFirstPartner = _age;
	}
	
	stopif(_nCurrSexPartner==_nMaxCurrSexPartner,
		   "Cannot add more partner than maximum allowed");
	
	_nCurrSexPartner++;			// One more concurrent sex partner
	_partnerUID.push_back(uid);	// Vector of partner increased
	_nLifetimePartner++;		// One more sex partner in lifetime
	_partnershipDuration.push_back(0); // New partnership just began => duration=0
	_isSpousalPartner.push_back(0);		// New partnership is not spousal at begining
}


void Individual::erasePartner(unsigned long uid2)
{
	/// ERASE A PARTNER
	/// _ONLY_ FROM THIS INDIVIDUAL'S POINT OF VIEW
	
	// Check if they are indeed coupled
	bool check=false;
	unsigned long p2=0; // index of partner uid2
	
	for (int p=0; p<_nCurrSexPartner; p++)	{
		if (_partnerUID[p]==uid2){
			check=true;
			p2 = p;
		}
	}
	
	if (!check)	{
		cout << "ERROR [erasePartner]: trying to dissolve a partnership that is not formed! ("
		<< _UID <<"-X-"<<uid2<<")"<<endl;
		displayInfo();
		exit(1);
	}
	
	bool isAspouse = isSpouse(uid2);

	// Cumbersome thing to code to access "erase" method in vectors ...
	// Remove all elements of vectors that contained info on this dissolved partnership
	
	vector<unsigned long>::iterator p2_uid = _partnerUID.begin()+p2;
	vector<double>::iterator p2_pd = _partnershipDuration.begin()+p2;
	vector<bool>::iterator p2_sp = _isSpousalPartner.begin()+p2;
	
	// Erase relevant element in the list of partners
	_partnerUID.erase(p2_uid);			// Erase partner's UID from the list of partners
	_partnershipDuration.erase(p2_pd);	// Erase duration with this partner
	_isSpousalPartner.erase(p2_sp);		// Erase spousal info with this partner
	
	_nCurrSexPartner-- ; // Number of current partners decreased by 1
	
	if (isAspouse) {
		_nCurrSpouse --;	// If it was a spouse, then decrease count by one
		_isDivorced = true;
	}
}


void Individual::convertToSpouse(unsigned long uid_partner)
{
	int ip = getPartnershipPosition(uid_partner);
	setIsSpousalPartner(ip, true);
	increase_nCurrSpouse();
}


// =========================================
// ============= SET FUNCTIONS =============
// =========================================

void Individual::setIsSpousalPartner(int p, bool x)
{
	if (p >= _nCurrSexPartner)
	{
		cout << endl << "ERROR [Individual::setIsSpousalPartner]";
		cout << "cannot set spousal status to partnership that doesn't exist!"<<endl;
		exit(1);
	}
	// DEBUG:
	//cout << endl <<	"Setting partner #" << p << " spouse="<< x << " for UID " << _UID<<endl;

	vector<bool> tmp = _isSpousalPartner;
	
	tmp[p] = x;
	
	_isSpousalPartner = tmp;
}	

void Individual::set_UID_n_sexAct_Type_period(vector<unsigned int> nSexType)
{
	string errmsg = "vector size of sex types must be 3";
	stopif(nSexType.size()!=3,errmsg);
	
	_UID_n_sexAct_Type0_period.push_back(nSexType[0]);
	_UID_n_sexAct_Type1_period.push_back(nSexType[1]);
	_UID_n_sexAct_Type2_period.push_back(nSexType[2]);
}




// =========================================
// ============= GET FUNCTIONS =============
// =========================================


int Individual::getPartnershipPosition(unsigned long uid_partner)
{
	// Retrieve the position of a partnership, given partner's uid
	
	int n = _nCurrSexPartner;
	
	int i;
	for (i=0; i<n; i++) 
	{
		if (getPartnerUID(i)==uid_partner) break;
	}
	
	string errmsg="Can't find position of partnership b/w "+to_string(_UID)+" and "+to_string(uid_partner);
	stopif(i==n,errmsg);
	
	return i;
}


vector<unsigned long> Individual::getPartnerUID() 
{
	return _partnerUID;
}


unsigned long Individual::getPartnerUID(int nPartner) 
{
	unsigned long res = 0;
	
	if (_nCurrSexPartner>0) res = _partnerUID[nPartner];
	return res;
}


vector<unsigned long> Individual::getSpouseUID() 
{
	vector<unsigned long> res;
	
	for (int i=0; i<_nCurrSexPartner; i++)
	{
		if (_isSpousalPartner[i]) 
		{
			res.push_back(_partnerUID[i]);
		}
	}
	
	if (res.size() != _nCurrSpouse)
	{
		cout << endl << "ERROR [Individual::getSpouseUID]: incoherence between '_nCurrSpouse' and '_isSpousalPartner'"<<endl;
		exit(1);
	}
	
	return res;
}


vector<unsigned long> Individual::getCasualUID() 
{
	vector<unsigned long> res;
	
	for (int i=0; i<_nCurrSexPartner; i++)
	{
		if (!_isSpousalPartner[i]) 
		{
			res.push_back(_partnerUID[i]);
		}
	}
	
	if (res.size() != _nCurrSexPartner - _nCurrSpouse)
	{
		cout << endl << "ERROR [Individual::getCasualUID]: incoherence between '_nCurrSpouse' and '_isSpousalPartner'"<<endl;
		exit(1);
	}
	
	return res;
}


double Individual::get_STIduration(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STIduration[i];
}


double Individual::get_STIsusceptFactor(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STIsusceptFactor[i];
}

bool Individual::get_STIsymptom(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STIsymptom[i];
}

double Individual::get_STItreatDuration(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STItreatDuration[i];
}

int Individual::get_STItreatTMS(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STItreatTMS[i];
}

double Individual::get_STItreatAdherence(STIname stiname)
{
	int i = positionSTIinVector(stiname, _STI);
	return _STItreatAdherence[i];
}

double Individual::get_STI_immunity(STIname stiname)
{
	/// Immunity to 'stiname'
	int i = positionSTIinVector(stiname, _STI);
	return _STI_immunity[i];
}


vector<unsigned long> Individual::get_STI_secondary_cases(STIname stiname){
	
	int i_sti = STI_find_index_position(stiname);
	return _STI_secondary_cases[i_sti];
}




// =====================================
// ========= SET FUNCTIONS =============
// =====================================



void Individual::set_nCurrSexPartner(int n) 
{	
	_nCurrSexPartner = n; 
	
	// Partnerships duration initialized at 0
	vector<double> pd(_nCurrSexPartner,0);
	_partnershipDuration = pd;
}

void Individual::decrease_nCurrSexPartner()
{
	string errmsg = "cannot decrease _nCurrSexPartner because = 0 (uid:"+to_string(_UID) + ")";
	stopif (_nCurrSexPartner==0, errmsg);
	_nCurrSexPartner--;
}


void Individual::initPartnershipDuration()
{
	/// Initiate vector when population is created
	
	vector<double> pd(_nCurrSexPartner,0);
	_partnershipDuration = pd;
}


void Individual::setPartnershipDuration(int n, double d)
{
	if (n >= _nCurrSexPartner)
	{
		cout << endl << "ERROR [setPartnershipDuration]: cannot assign a partnership duration that does not exist"<< endl;
		exit(1);
	}
	
	_partnershipDuration[n] = d;
}


void Individual::increasePartnershipDuration(int n, double timeStep)
{
	if (n >= _nCurrSexPartner)
	{
		cout << endl << "ERROR [increasePartnershipDuration]: cannot increase a partnership duration that does not exist"<< endl;
		exit(1);
	}
	
	_partnershipDuration[n] += timeStep;
}


bool Individual::isSpouse(unsigned long uid)
{
	if (_nCurrSexPartner==0) return false;
	
	int ip = getPartnershipPosition(uid);
	return get_isSpousalPartner()[ip];
}




bool Individual::isOpenToNewPartnership()
{
	return (get_nCurrSexPartner() < get_nMaxCurrSexPartner());
}


void Individual::add_sexAct(unsigned long uid_partner, int n)
{
	/// ADD n SEX ACTS BETWEEN THIS INDIVIDUAL AND uid_partner
	
	// add the UID of the partner this individual has sex with
	_UID_sexAct_period.push_back(uid_partner);
	_UID_n_sexAct_period.push_back(n);
	
	// The increase count of sex acts is already 
	// updated in males (needed at the start,
	// before distributing sex acts to partner types etc)
	// Hence, increase count only for females here
	if (_gender==female) _nSexActs_period += n;
	
	// If it's the first time ever have sex, then record age
	// (by default, _ageFirstSex=-99)
	if (_ageFirstSex<0) _ageFirstSex = _age;
	
	
	// DEBUG:
//	if (_ageFirstPartner+1<_ageFirstSex && _ageFirstPartner>0 && _UID>2000)
//	{
//		cout << endl<<"UID="<<_UID<< " ; gender=" <<_gender
//		<< " ; _ageFirstPartner="<<_ageFirstPartner
//		<< " ; _ageFirstSex="<<_ageFirstSex<<endl;
//	}
	//cout << "add "<< n <<" sexact: " << _UID << " +++ " << uid_partner;
	//cout <<" ; total: "<< _nSexActs_period <<endl;
}



void Individual::reset_UID_n_sexAct_Type()
{
	_UID_n_sexAct_Type0_period.clear();	
	_UID_n_sexAct_Type1_period.clear();	
	_UID_n_sexAct_Type2_period.clear();
}



void Individual::reset_UID_n_sexAct_period()
{
	_UID_n_sexAct_period.clear();
}



void Individual::init_STIsusceptFactor()
{
	_STIsusceptFactor.clear();
	int n = _STI.size();
	
	// susceptibility factor = 1.00
	// by default (multiplicative factor)
	_STIsusceptFactor.resize(n,1.00);
	
	// If circumcised, then susceptibility
	// decreased for some STIs
	if (_isCircum)	{
		for (int i=0;i<n;i++)
			_STIsusceptFactor[i]= _STI[i].get_circum_SF_reduction();
	}
	
	// If already vaccinated but
	// not immunized
	for(int i=0; i<n; i++){
		if(_STI_vacc[i]){
			double susc_vacc = 1-_STI_immunity[i];
			_STIsusceptFactor[i]=_STIsusceptFactor[i] * susc_vacc;
		}
	}
}



int Individual::nSexActType1or2()
{
	/// Total number (across all partners) of
	/// sex acts of type 1 or 2 (not type0[:condom])
	
	int res = 0;
	
	int n_type1 = sumElements(_UID_n_sexAct_Type1_period);
	int n_type2 = sumElements(_UID_n_sexAct_Type2_period);
	
	res = n_type1+n_type2;
	
	return res;
}


/* ************* STI **************/


void Individual::STI_initializeAll(vector<STI> STItemplate)
{
	/// INITIALIZE ALL STIs FEATURES
	/// FOR A GIVEN INDIVIDUAL
	
	_STI = STItemplate;
	
	int n = _STI.size();
	
	_STIduration.clear();
	_STIduration.resize(n,0);
	
	// -- Vaccination (none)
	
	_STI_vacc.clear();
	_STI_vacc.resize(_STI.size(),false);
	
	_STI_vacc_time.clear();
	_STI_vacc_time.resize(_STI.size(),99e9);
	
	_STI_immunity.clear();
	_STI_immunity.resize(_STI.size(),0.0);

	init_STIsusceptFactor();
	
	_STIsymptom.clear();
	_STIsymptom.resize(n,false);
	
	// Resize the vector holding the vectors of secondary cases
	// (one vector for each STI):
	_STI_secondary_cases.resize(n);
	
	// Initialize (to none) STI treatments
	
	vector<double> zero(STItemplate.size(), 0.0);
	vector<int> zeroint(STItemplate.size(), 0);
	set_STItreatDuration(zero);
	set_STItreatAdherence(zero);
	set_STItreatTMS(zeroint);
	
	
	
	// DEBUG:
	// cout << _STI.size() << " STI INITIALIZED"<<endl;
}



void Individual::set_STIduration(int sti, double duration)
{
	if (sti >= _STI.size()){
		cout << " ERROR [STI_setDuration]: number of STIs ("<< _STI.size()<<") defined is smaller than "<<sti<<endl;
		exit(1);
	}
	_STIduration[sti] = duration;
}



void Individual::set_STIduration(STIname stiname, double duration)
{
	bool found = false;
	int ss=-99;
	
	// Search for the STI name
	for (int s=0; s<_STI.size(); s++){
		if (_STI[s].get_name()==stiname){
			found = true;
			ss = s;
		}
	}
	string errmsg = "ERROR [STI_setDuration]: STI name not found: " + STInameString(stiname);
	stopif(!found,errmsg);
	
	// Set the duration value to the relevant STI
	_STIduration[ss] = duration;
}


void Individual::set_STIsymptom(int sti_index, bool isSymptomatic)
{
	_STIsymptom[sti_index] = isSymptomatic;
}


void Individual::set_STIsymptom(STIname stiname, bool isSymptomatic)
{
	int i_sti = positionSTIinVector(stiname, _STI);
	_STIsymptom[i_sti] = isSymptomatic;
}

void Individual::set_STI_MTCT(STIname stiname, bool mtct)
{
	int i_sti = positionSTIinVector(stiname, _STI);
	_STI_MTCT[i_sti] = mtct;
}


bool Individual::is_symptomatic()
{
	/// RETURNS TRUE IF THIS INDIVIDUAL
	/// HAS SYMPTOM TO _ANY_ STIs
	
	bool isSymptomatic = false;
	
	for (int i=0; i<_STIsymptom.size() && !isSymptomatic; i++){
		if (_STIsymptom[i]) {
			isSymptomatic=true;
			break;
		}
	}
	return isSymptomatic;
}

bool Individual::is_symptomatic(STIname stiname)
{
	/// RETURNS TRUE IF THIS INDIVIDUAL
	/// HAS SYMPTOM TO _THAT_ STIs
	
	int i_sti = positionSTIinVector(stiname, _STI);
	return _STIsymptom[i_sti];
}



void Individual::STI_resetAllDurations()
{
	for (int i=0; i<_STI.size(); i++) {
		set_STIduration(i, 0.0);
	}
}


void Individual::increaseSTIdurations(double prd)
{
	for (int i=0; i<_STIduration.size(); i++) 
	{
		if (_STIduration[i]>0) // condition needed, in order not to increase disease not acquired (when _STIduration=0)
			_STIduration[i] += prd;
	}
}


void Individual::STI_acquireInfection(STIname name, double durationSinceInfection)
{
	/// SET THE PARAMETERS VALUES (BOTH DETERMINISTIC AND STOCHASTIC)
	/// OF THIS STI INFECTING THIS INDIVIDUAL UPON INFECTION
	
	// Calculates the index corresponding to the STI name given
	int i = positionSTIinVector(name, _STI);
	
	string errmsg = to_string(_UID)+ " is already infected with "+ STInameString(name);
	stopif(_STIduration[i]>0, errmsg);
	
	// This STI has started infection
	_STIduration[i] = durationSinceInfection;
	
	// === SYMPTOMATIC STATUS ===
	
	// Determine if this STI will be symptomatic
	// (gender-based)
	
	// TO DO: try to redesign such that _STIsymptom is not independent
	
	double u		= uniform01();
	double pf		= _STI[i].get_proba_symptomatic_female();
	double pm		= _STI[i].get_proba_symptomatic_male();
	double proba	= (_gender==female)?pf:pm;
	
	bool symptom	= false;
	if (u<proba) symptom = true;
	_STI[i].set_is_symptomatic(symptom);
	_STIsymptom[i] = symptom;
	
	// DEBUG
	//cout << _UID << " Newly aqcuired "<<STInameString(_STI[sti].get_name())<<" symptomatic:"<<tmp<<endl;
	

	// === RECURRENCE/PERSISTENCE STATUS ===
	
	// Determine if recurrent symptoms will occur with this STI
	// By default, STI is not recurent
	
	_STI[i].set_is_recurrent(false);
	
	// Only for recurrent/persistent STIs

	if ((name==HSV2) ||
		(name==HPV)){
		u = uniform01();
		if (u<_STI[i].get_proba_recurrence())
		{
			_STI[i].set_is_recurrent(true);
			// DEBUG
			//cout << _UID << "Newly aqcuired "<<STInameString(_STI[sti].get_name())<<" is recurrent."<<endl;
		}
	}
	
	
	// === OTHER STOCHASTIC PARAMETERS ===
	
	
	if (name == HPV)
	{
		// Variables _latentDuration and _infectiousDuration
		// are re-defined in the case of HPV because they are input
		// as maximum values these random variables can take (for each individual)
		
		// latent -> not infectious
		// For HPV, this value represents the MAXIMUM latent duration
		// Draw random variable for latent duration ~ Beta(3,3)
		
		double incub_max = _STI[i].get_latentDuration();
		_STI[i].set_latentDuration(incub_max * beta(3,3));

		// For HPV, "_infectiousDuration" represents the MAXIMUM infectious duration
		
		double infect_max = _STI[i].get_infectiousDuration();
		_STI[i].set_infectiousDuration(infect_max * beta(3,3));
		
//		cout <<"DEBUG -- HPV incub="<<_STI[sti].get_latentDuration();
//		cout << "   HPV infec="<< _STI[sti].get_infectiousDuration() <<endl;
	}
}



void Individual::add_STI_secondary_cases(STIname stiname,
										 unsigned long uid_infected){
	/// Add the UID of a secondary case for a given STI
	
	int i_sti = STI_find_index_position(stiname);
	_STI_secondary_cases[i_sti].push_back(uid_infected);
}


void Individual::clear_STI_secondary_cases(STIname stiname){
	/// Clear the record of all secondary cases
	/// for that STI (usually following cure of this STI)

	int i_sti = STI_find_index_position(stiname);
	_STI_secondary_cases[i_sti].clear();
}


bool Individual::STI_infected(STIname stiname){
	/// Check if this individual is being infected with that STI
	return (get_STIduration(stiname)>0?true:false);
}


bool Individual::STI_treated(STIname stiname){
	/// Check if this individual is being treated against that STI
	return (get_STItreatDuration(stiname)>0?true:false);
}



bool Individual::STI_anyInfection()
{
	/// Returns TRUE if this individual as at least one STI

	bool infected = false;	
	int i=0;
	
	while ( !infected && i<_STIduration.size() ) {
		infected = _STIduration[i]>0 ? true : infected ;
		i++;
	}	
	return infected;
}


int Individual::STI_nModelled(){
	// Returns the number of STIs modelled
	return _STI.size();
}


int Individual::STI_find_index_position(STIname s){
	/// FOR A GIVEN STI NAME
	/// RETURNS ITS POSITION
	/// (FROM THE VECTOR '_STI')
	
	int	res = -99;
	bool found = false;
	
	for (int i=0;i<_STI.size();i++){
		if (_STI[i].get_name()==s){
			res=i;
			found=true;
		}
	}
	string errmsg = " ERROR [STI_find_index_position]: STI name not found: " + STInameString(s);
	stopif(!found, errmsg);
	return res;
}


void Individual::STI_updateSusceptFactor()
{
	// TO DO
}


vector<string>	Individual::STI_listInfection()
{
	vector<string>	res(0);
	
	if ( STI_anyInfection() ){
		int n = _STIduration.size();
		for (int i=0; i<n; i++) {
			if (_STIduration[i]>0) {
				res.push_back(STInameString(_STI[i].get_name()));
			}
		}
	}
	return res;	
}

vector<double> Individual::STI_IC()
{
	/// Returns the infectivity to all STIs modeled.
	/// (assume durations are updated [in simulation])
	
	// Number of STis modelled
	int nSTI = _STIduration.size();
	
	// Infectivity for each STI
	// (at the current simulation time)
	vector<double> IC(nSTI,0.0);
	
	for (int i=0; i<nSTI; i++) {
		if (_STIduration[i]>0) 
		{
			// Retrieve infectivity curve
			IC[i] = _STI[i].infectivityCurve(_STIduration[i], _gender);
			
			// If this individual is under treatment for this STI
			// then adjust the infectivity level
			// (if no treatment, treat_reduction = 1)
			double treat_reduction = TRE(_STI[i].get_name());
			// If this individual has been vaccinated
			// with a vaccine that does not provide immunity
			// but reduces permanently infectiousness
			double vax_reduction = VRE(_STI[i].get_name());
			
			IC[i] = IC[i] * treat_reduction * vax_reduction;
			
			if (_STI[i].get_name() == HIV){
				// The increased HIV infectiousness
				// due to co-infection is implemented here.
				
				double incr_infect = 0.0;
				
				// DEBUG
//				cout << _UID << " has HIV since "<< get_STIduration(HIV)<<"; base IC="<<IC[i];
//				cout << " ; treatment: "<< _STItreatDuration[i]<<" ; treat_reduction=" << treat_reduction <<endl;
				
				for (int j=0; j<nSTI; j++){
					// Scan all other STIs because
					// they will increase HIV infectivity
					double Rebound_ji=0;
					double IC_j=0;
                    // test if co-infection with another STI
					if (i!=j && _STIduration[j]>0){
						Rebound_ji = _RebHIV[j];
						IC_j = _STI[j].infectivityCurve(_STIduration[j],_gender);
						
						incr_infect += Rebound_ji*IC_j;
						
						//DEBUG
//						cout <<"co-infected with ";
//						cout<< STInameString(_STI[j].get_name());
//						cout << " since :"<< _STIduration[j]<<" IC2="<<IC_j <<endl;
					}
				}
				// Integrity check:
				stopif(incr_infect<0, "Infectiousness increase is supposed to be positive.");
				
				// Final infectivity curve for STI "i"
				// taking into account all other co-infections:
				IC[i] = min(1.0, IC[i]*(1+incr_infect));
				
				// DEBUG
//				cout << endl << _UID << " increased HIV IC="<<IC[i]<<endl<<endl;
			}
		}
	}
	return IC;
}



void Individual::display_STI_IC()
{
	int n = _STI.size();
	cout << endl << " === INFECTIVITY CURVE FOR UID " <<_UID <<" ===" <<endl;
	for (int i=0; i<n; i++){
		cout << STInameString(_STI[i].get_name()) << " --> " << STI_IC()[i] <<endl;
	}
}



/* **************************************** */
/* ************* STI TREATMENT **************/
/* **************************************** */


void Individual::set_STItreatDuration(double x, STIname sti)
{
	stopif(x<0,"treatment duration cannot be negative!");
	
	int i_sti = positionSTIinVector(sti, _STI);
	_STItreatDuration[i_sti] = x;
}



void Individual::increase_STItreatDuration(double prd)
{
	for (int i=0; i<_STItreatDuration.size(); i++){
		if (_STItreatDuration[i]>0)
			_STItreatDuration[i]+= prd;
	}
}



void Individual::set_STItreatTMS(int x, STIname sti)
{
	int i_sti = positionSTIinVector(sti, _STI);
	_STItreatTMS[i_sti] = x;
}



void Individual::set_STItreatAdherence(double x, STIname sti)
{
	int i_sti = positionSTIinVector(sti, _STI);
	_STItreatAdherence[i_sti] = x;
}



double Individual::TRE(STIname sti)
{
	/// TREATMENT REDUCTION EFFECT ("TRE")
	/// CALCULATES MULTIPLICATIVE FACTOR
	/// THAT WILL REDUCE STI's INFECTIVITY CURVE
	/// THANKS TO TREATMENT
	
	double res=-9e9;
	int i_sti = positionSTIinVector(sti, _STI);
	double t = _STItreatDuration[i_sti];
	
	// If no treatment, no infectiousness reduction
	if (t==0) res = 1.0;
	
	// If treatment started
	if (t>0)
	{
		// retrieve treatment microbiological success flag
		// (assuming perfect adherence)
		int TMS = _STItreatTMS[i_sti];
		
		double adherence = _STItreatAdherence[i_sti];
		double tmp1 = TMS*(adherence*_STI[i_sti].TREstar(t) + 1-adherence);
		double tmp2 = 1-TMS;
		res = tmp1+tmp2;
		
		// DEBUG -------------------------------
		if(0){
			string fname = _DIR_OUT + "treat_attempts.out";
			ofstream f(fname,ios::app);
			
			f <<_UID<<",";
			f <<STInameString(sti)<<",";
			f <<TMS<<",";
			f <<adherence<<",";
			f <<_STI[i_sti].TREstar(t)<<",";
		}
		// -------------------------------------
	}
	return res;
}



double Individual::VRE(STIname stiname){
	/// VACCINATION REDUCTION EFFECT ("VRE")
	/// CALCULATES MULTIPLICATIVE FACTOR
	/// THAT WILL REDUCE STI's INFECTIVITY CURVE
	/// THANKS TO VACCINATION (in case vax does not provide immunity)
	
	double res	= -9e9;
	int i_sti	= positionSTIinVector(stiname, _STI);
	bool isVax	= _STI_vacc[i_sti];
	
	// If no vaccine, no infectiousness reduction
	if (!isVax) res = 1.0;
	// Else, retrieve the reduction value
	if (isVax)
		res = _STI[i_sti].get_VRE();
	
	return res;
}



void Individual::treat(STIname sti)
{
    /// ATTEMTPS TO TREAT STI
    /// FOR THIS INDIVIDUAL
    /// (STOCHASTIC)
	
	/// DETERMINES (FOR THIS INDIVIDUAL):
	/// - microbiological success of treatment
	/// - adherence
    
    int i_sti = positionSTIinVector(sti, _STI);
    
    // Check first if STI indeed infecting
	string errmsg ="Trying to treat STI "+STInameString(sti)+" that is not infecting UID "+ to_string(_UID);
	stopif(get_STIduration(sti)==0, errmsg);
	
	// Initialize treatment duration to
	// a small value (unit in years)
	// (being>0 is a test that treatment has started)
	set_STItreatDuration(0.001,sti);
	
    // --- Microbiological success of treatment ---
	// (success is stochastic and depends on each individual)
    double	p_fail	= _STI[i_sti].get_proba_treatmentFailure();
    int		TMS		= binom(1-p_fail, 1);

	set_STItreatTMS(TMS, sti);
    
    // --- Adherence ---
	// adh_prm[0] = adherence_max
	// adh_prm[1] = decay risk group
	// adh_prm[0] = asymptomatic factor
	vector<double> adh_prm = _STI[i_sti].get_adherence_param();
	
	// if symptomatic
    double A = adh_prm[0]*exp(-adh_prm[1]*(double)(_riskGroup));
	
	// if asymptomatic
	if (_STI[i_sti].get_is_symptomatic()) A = A*adh_prm[2];
	
	set_STItreatAdherence(A, sti);
}



/* ************* MISCELLEANOUS **************/







/* ************* HELPERS ************* */

void Individual::displayInfo()
{
	cout << endl << "==== Info on Individual ===="<<endl;
	
	cout << "UID: " <<_UID <<endl;
	cout << "Alive: " << _isAlive <<endl;
	cout << "Gender: " << _gender <<endl;
	cout << "Age: " << _age <<endl;
	
	cout << "-- Sexual Activity --" << endl;
	cout << "Risk group: " << _riskGroup << endl;
	cout << "nSexPartners: " << _nCurrSexPartner <<endl;
	cout << "nSpouses: "<< _nCurrSpouse <<endl;
	cout << "maxSexPartners: " << _nMaxCurrSexPartner <<endl;
	cout << "nSexActLastPrd: " << _nSexActs_period <<endl;
	
	if (_nCurrSexPartner>0)
	{
		cout << "UID of partners:"; 
		displayVector(_partnerUID);
	}
	
	cout << " -- STI --" <<endl;
	cout << "STIs modeled: ";
	for (int i=0; i<_STI.size(); i++) {
		cout<<STInameString(_STI[i].get_name())<<"; ";
	}
	cout<<endl<<endl;

	if (STI_anyInfection()) 
	{
		cout << "Infected with the following STIs:";
		vector<string> stiInfected = STI_listInfection();
		displayVector(stiInfected);
	 	cout << "STIs durations:"<< endl;
		displayVector(_STIduration);
	}
	else { cout<<"No STI infections" << endl; }

	cout <<endl<< "STI treatment durations in years (0=none)"<<endl;
	displayVector(_STItreatDuration);
	
	coutline(40);

}
