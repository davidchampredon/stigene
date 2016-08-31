//
//  STI.cpp
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 13-09-15.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#include "STI.h"


STI::STI(STIname name, string filename)
{
	/// ============================
	/// CONSTRUCT THE STI
	/// FROM A FILE
	/// ============================
	
	// TO DO: finish exhaustive implementation
	
	
	_name = name;
	
	_probaSexTransm_SAT.resize(3,0);
	_probaSexTransm_SAT[0] = 0.01; // Reduction in transmission when condom is used
	_probaSexTransm_SAT[2] = 1.00; // No reduction in transmission when riskiest sex act is performed
	// ** WARNING ** _probaSexTransm_SAT[1] is defined specifically for each STI
	
	
	switch (name)
	{
		case HIV:
			
			_type = virus;
			
			load_common_param(name, filename);
			
			_proba_recurrence =  0.0;
			
			_HIVparam.resize(7);
			_HIVparam[0] = getParameterFromFile("HIV_peakTime_weeks",filename)/52.0;	// Time after initial infection when viral load peaks
			_HIVparam[1] = getParameterFromFile("HIV_shape_acute",filename) ;			// Shape parameter of acute infection
			_HIVparam[2] = getParameterFromFile("HIV_var_acute_weeks",filename)/52.0;	// Dispersion of the duration of acute phase
			_HIVparam[3] = getParameterFromFile("HIV_frac_chronic_acute",filename) ;	// Fraction of peak viral load when chronic stage starts
			_HIVparam[4] = getParameterFromFile("HIV_chronic_duration_yrs",filename) ;	// Duration (in years) of the chronic infectious stage
			_HIVparam[5] = getParameterFromFile("HIV_chronic_growth",filename) ;		// Rate of viral load progression during the chronic stage
			_HIVparam[6] = getParameterFromFile("HIV_aids_duration_yrs",filename) ;		// Duration (in years) of AIDS (death as end-point)
			
			break;
			
			
		case HSV2:
			
			_type = virus;
			
			load_common_param(name, filename);
			
			_proba_recurrence		= getParameterFromFile("HSV2_proba_recurrence",filename);
			_recurrence_freq		= getParameterFromFile("HSV2_recurrence_freq",filename);
			_recurrence_duration	= getParameterFromFile("HSV2_recurrence_duration",filename);
			
			break;
		
			
			
		case HPV:
			
			_type = virus;
			
			load_common_param(name, filename);
			
			_proba_recurrence		= getParameterFromFile("HPV_proba_recurrence",filename);
			
			break;

			
			
		case Ct:
			_type = bacteria;
			load_common_param(name, filename);
			
			_proba_recurrence =  0.0;
			
			_infectiousDuration_2	=  getParameterFromFile("Ct_infectious_weeks_symptom", filename)/52.0;
			
			_shape_param.resize(4);
			
			_shape_param[0]			= getParameterFromFile("Ct_shape_a", filename);
			_shape_param[1]			= getParameterFromFile("Ct_shape_b", filename);
			
			_shape_param[2]			= getParameterFromFile("Ct_shape_a_symptom", filename);
			_shape_param[3]			= getParameterFromFile("Ct_shape_b_symptom", filename);
			
			break;
			
		
			
		case Ng:
			_type = bacteria;
			load_common_param(name, filename);
			_proba_recurrence =  0.0;
			
			_shape_param.resize(2);
			
			_shape_param[0]			= getParameterFromFile("Ng_shape_a", filename);
			_shape_param[1]			= getParameterFromFile("Ng_shape_b", filename);
			
			break;
		
	
			
		case Tp:
			_type = bacteria;
			load_common_param(name, filename);
			
			_proba_recurrence		=  getParameterFromFile("Tp_proba_recurrence", filename);
			_proba_recurrence_2		=  getParameterFromFile("Tp_proba_recurrence_2", filename);
			
			_latentDuration_2	=  getParameterFromFile("Tp_latent_weeks_2", filename)/52.0;
			_latentDuration_3	=  getParameterFromFile("Tp_latent_weeks_3", filename)/52.0;
			
			_infectiousDuration_2	=  getParameterFromFile("Tp_infectious_weeks_2", filename)/52.0;
			_infectiousDuration_2c	=  getParameterFromFile("Tp_infectious_weeks_2c", filename)/52.0;
			_infectiousDuration_3	=  getParameterFromFile("Tp_infectious_weeks_3", filename)/52.0;
			
			_maxInfectiousness_Tp1	= getParameterFromFile("Tp_maxInfectiousness_Tp1", filename);
			_maxInfectiousness_Tp2	= getParameterFromFile("Tp_maxInfectiousness_Tp2", filename);
			_maxInfectiousness_TpEL	= getParameterFromFile("Tp_maxInfectiousness_TpEL", filename);
			
			
			_shape_param.resize(8);
			
			_shape_param[0]			= getParameterFromFile("Tp_shape_a1", filename);
			_shape_param[1]			= getParameterFromFile("Tp_shape_b1", filename);
			
			_shape_param[2]			= getParameterFromFile("Tp_shape_a2", filename);
			_shape_param[3]			= getParameterFromFile("Tp_shape_b2", filename);
			
			_shape_param[4]			= getParameterFromFile("Tp_shape_a2c", filename);
			_shape_param[5]			= getParameterFromFile("Tp_shape_b2c", filename);
			
			_shape_param[6]			= getParameterFromFile("Tp_shape_aEL", filename);
			_shape_param[7]			= getParameterFromFile("Tp_shape_bEL", filename);
			
			break;
		
	
		
		
		case Hd:
			_type = bacteria;
			load_common_param(name, filename);
			_proba_recurrence =  0.0;
			
			_infectiousDuration_2	=  getParameterFromFile("Hd_infectious_weeks_symptom", filename)/52.0;
			
			_shape_param.resize(4);
			
			_shape_param[0]			= getParameterFromFile("Hd_shape_a", filename);
			_shape_param[1]			= getParameterFromFile("Hd_shape_b", filename);
			
			_shape_param[2]			= getParameterFromFile("Hd_shape_a_symptom", filename);
			_shape_param[3]			= getParameterFromFile("Hd_shape_b_symptom", filename);
			
			break;
			
			
			
		case Bv:
			_type = bacteria;
			load_common_param(name, filename);
			_proba_recurrence =  0.0;
			break;
		
		case Tv:
			_type = protozoa;
			load_common_param(name, filename);
			_proba_recurrence =  0.0;
			
			_infectiousDuration_f_a = getParameterFromFile("Tv_infectiousDuration_f_a_months", filename)/12.0;
			_infectiousDuration_f_s = getParameterFromFile("Tv_infectiousDuration_f_s_months", filename)/12.0;
			_infectiousDuration_m_a = getParameterFromFile("Tv_infectiousDuration_m_a_months", filename)/12.0;
			_infectiousDuration_m_s = getParameterFromFile("Tv_infectiousDuration_m_s_months", filename)/12.0;
			
			_asymptom_IC_reduc		= getParameterFromFile("Tv_asymptom_IC_reduc", filename);
			
			_shape_param.resize(8);
			
			_shape_param[0]			= getParameterFromFile("Tv_shape1_fem_asympt", filename);
			_shape_param[1]			= getParameterFromFile("Tv_shape2_fem_asympt", filename);
			
			_shape_param[2]			= getParameterFromFile("Tv_shape1_fem_sympt", filename);
			_shape_param[3]			= getParameterFromFile("Tv_shape2_fem_sympt", filename);
			
			_shape_param[4]			= getParameterFromFile("Tv_shape1_male_asympt", filename);
			_shape_param[5]			= getParameterFromFile("Tv_shape2_male_asympt", filename);
			
			_shape_param[6]			= getParameterFromFile("Tv_shape1_male_sympt", filename);
			_shape_param[7]			= getParameterFromFile("Tv_shape2_male_sympt", filename);
			break;
		
			
		case Cv:
			_type = fungus;
			// TO FINISH...
			break;
			
		default:
			string errmsg = "STI name unknown!";
			stopif(true, errmsg);
			break;
	}
}


STI::STI(STIname name,
		 string STIfeatures_filename,
		 string treatment_filename)
{
	STI(name,STIfeatures_filename);
	load_treatment_param(name, treatment_filename);
}





// ==============================================
// ==============================================
// ==============================================
// ==============================================




void STI::load_common_param(STIname sti_name, string filename)
{
	/// LOAD PARAMETERS THAT ARE DEFINED
	/// IN ALL STIs
	
	string sti_string = STInameString(sti_name);
	
	_proba_symptomatic_female	= getParameterFromFile(sti_string + "_proba_symptomatic_female",filename);
	_proba_symptomatic_male		= getParameterFromFile(sti_string + "_proba_symptomatic_male",filename);
	
	if (sti_name != HIV){
		_latentDuration		= getParameterFromFile(sti_string + "_latent_weeks",filename)/52.0;
		_infectiousDuration	= getParameterFromFile(sti_string + "_infectious_weeks",filename)/52.0; //not used for Tv
		_minInfectiousness	= getParameterFromFile(sti_string + "_minInfectiousness", filename);
	}
	_probaMaxSexTransm		= getParameterFromFile(sti_string + "_probaMaxSexTransm",filename); 
	_probaSexTransm_SAT[1]	= getParameterFromFile(sti_string + "_reduc_probaTrans_loRisk",filename); // Reduction factor in transmission when no condom, low risk sex
	
	_proba_MTCT				= getParameterFromFile(sti_string + "_proba_MTCT",filename);
	_circum_SF_reduction	= getParameterFromFile(sti_string + "_circum_SF_reduction",filename) ;
	_naturalClearanceDuration = getParameterFromFile(sti_string + "_clearance_duration",filename);
}



void STI::load_treatment_param(STIname sti_name, string filename)
{
	/// LOAD TREATMENT RELATED PARAMETERS
	/// FOR THIS STI
	/// - adhrence parameters
	/// - proba treatment biological failure
	/// - optimal treatment duration
	
	bool debug = false;
	if(debug){
		cout << "Loading treatment parameters for "<<STInameString(sti_name);
		cout << " from file "<<filename<<endl;
	}
	
	// STI independent variables
	vector<double> adherence_param;
	adherence_param.push_back(getParameterFromFile("adherence_max",filename));
	adherence_param.push_back(getParameterFromFile("decay_riskgroup",filename));
	adherence_param.push_back(getParameterFromFile("asymptom_reduct_factor",filename));
	set_adherence_param(adherence_param);
	
	// STI dependent variables
	string sti_string = STInameString(sti_name);
	double ptf	= getParameterFromFile(sti_string + "_treat_fail",filename);
	double otd	= getParameterFromFile(sti_string + "_treat_optDur",filename);

	set_proba_treatmentFailure(ptf);
	set_optimalTreatmentDuration(otd);
	
	if(debug){
		cout << "adherence:";displayVector(adherence_param);
		cout << "proba treatment failure: "<<ptf<<endl;
		cout << "optimal duration: "<<otd<<endl;
	}
}



void STI::load_vaccine_param(STIname sti_name, string filename)
{
	/// LOAD VACCINE RELATED PARAMETERS FOR THIS STI
	
	bool debug = false;
	if(debug){
		cout << "Loading vaccine parameters for " << STInameString(sti_name);
		cout << " from file " << filename << endl;
	}

	string sti_string	= STInameString(sti_name);
	
	// Vaccine failure probability
	double pvf			= getParameterFromFile(sti_string + "_vaccine_fail",filename);
	set_proba_vaccineFailure(pvf);
	
	// Vaccine Reduction Effect on infectivity curve
	// when no immunity is provided
	double VRE			= getParameterFromFile(sti_string+"_VRE", filename);
	set_VRE(VRE);

	// Vaccine Reduction Effect on susceptibility
	// when no immunity is provided
	double vacc_SF_reduc = getParameterFromFile(sti_string+"_vacc_waneRate", filename);
	
	// --- DEBUG
	//cout << " DEBUG-wane-" << sti_string	<< ": "<< vacc_SF_reduc ;
	//cout << " from file: "<< filename << endl;
	// ---------
	
	set_vacc_waneRate(vacc_SF_reduc);
}



double STI::infectivityCurve(double stiDuration, Gender g)
{
	/// DEFINES HERE THE INFECTIVITY CURVE
	/// OF EVERY STI
	
	double res = 0;
	
	// ========= HIV ===========
	
	if (_name == HIV)
	{
		// Retrieve parameters
		double Tvlmax = _HIVparam[0] ;		// Time after initial infection when viral load peaks
		double q = _HIVparam[1];
		double sigma = _HIVparam[2];		// Dispersion of the duration of acute phase
		double VLc = _HIVparam[3];			// Fraction of peak viral load when chronic stage starts
		double D = _HIVparam[4];			// Duration (in years) of the chronic infectious stage
		double Rchr = _HIVparam[5];			// Rate of viral load progression during the chronic stage
		double Daids = _HIVparam[6];		// Duration (in years) of AIDS (death as end-point)
		//double VLcoSTI = _HIVparam[7];	// Fraction of peak viral load when STI co-infected
		
		// Calculate time when chronic stage starts
		double Tc = sigma*pow(-log(VLc),1/2/q)+Tvlmax;
		
		if (stiDuration<Tc){
			// ACUTE STAGE
			res = exp(-pow((stiDuration-Tvlmax)/sigma,2*q));
		}
		
		double Taids = Tc+D;
		
		if (Tc <= stiDuration && stiDuration<Taids){
			// CHRONIC STAGE
			res = VLc*exp((stiDuration - Tc)*Rchr);
		}
		
		if (Taids <= stiDuration){
			// AIDS
			double tmp = VLc*exp((Taids - Tc)*Rchr);
			res = min(1.0,tmp * exp(-(stiDuration-Taids)*log(VLc)/Daids));
		}	
	}
	
	
	
	// ============== HSV2 ===============
	
	if (_name == HSV2)
	{		
		// latent -> not infectious
		if (stiDuration <= _latentDuration) res = 0.0;
		
		// First episode
		
		double first_recurrence = _latentDuration+_infectiousDuration;//+1.0/_recurrence_freq;
		
		if (stiDuration > _latentDuration && stiDuration < 2*first_recurrence){
			// first lesion last twice as long as recurrent lesion (hence _recurrence_duration/2)
			double tmp = pow(stiDuration-_latentDuration-_infectiousDuration/2.0,2.0)/_recurrence_duration/2.0;
			
			res = exp(-tmp)*(1.0-_minInfectiousness) + _minInfectiousness;
		}
		
		
		// Recurrent lesions
		
		if (stiDuration >= 2*first_recurrence)
		{
			double w;
			double sig;
			
			if (_is_symptomatic){
				w	= _recurrence_freq;
				sig	= _recurrence_duration;
			}
			
			if (!_is_symptomatic){
				w	= _recurrence_freq*0.66;
				sig	= _recurrence_duration/100;
			}
			
			double tmp2 =pow(sin(PI* stiDuration* w), 2.0)/sig;
			tmp2 = (exp(tmp2)-1)/(exp(1/sig)-1);
			tmp2 = tmp2*(1-_minInfectiousness)+_minInfectiousness;
			res = tmp2;
		}
	}
	
	
	
	// ============== HPV ===============
	
	if (_name == HPV)
	{
		// Shape parameter
		double k_hpv = 5.0;
		// normalizing constant
		double eta0 = exp(-1)/k_hpv;
		// Infectivity curve if not persistent HPV infection
		double t2 = max(stiDuration-_latentDuration,0.0)/_infectiousDuration;
		
		// WARNING: for HPV, recurrent means PERSISTENT
		if(!_is_recurrent){
			res = t2*exp(-k_hpv*t2)/eta0;
		}
			
		if(_is_recurrent){
			res = 1-exp(-t2/eta0);
		}
	}
	
	
	
	// ============== CHLAMYDIA ===============
	
	if (_name == Ct)
	{
		// Shape parameter depends on symptomatic status
		double a,b,t;
		
		if(!_is_symptomatic){
			t = max(stiDuration -_latentDuration,0.0)/_infectiousDuration;
			a = _shape_param[0];
			b = _shape_param[1];
		}

		if(_is_symptomatic){
			t = max(stiDuration -_latentDuration,0.0)/_infectiousDuration_2;
			a = _shape_param[2];
			b = _shape_param[3];
		}

		vector<double> prm;
		prm.push_back(a);
		prm.push_back(b);
				
		res = IC("pseudo_beta",t,prm); // pseudo_beta(t, a, b);
	}
	
	
	// ============== GONORROHEA ===============
	
	if (_name == Ng){
		double t = max(stiDuration -_latentDuration,0.0)/_infectiousDuration;
		res = IC("pseudo_beta",t,_shape_param); //pseudo_beta(t, a, b);
	}
	
	
	
	// ============== Chancroids ===============
	
	if (_name == Hd)
	{
		// Shape parameter depends on symptomatic status
		
		double a,b,t;
		
		if(!_is_symptomatic){
			t = max(stiDuration -_latentDuration,0.0)/_infectiousDuration;
			a = _shape_param[0];
			b = _shape_param[1];
		}
		
		if(_is_symptomatic){
			t = max(stiDuration -_latentDuration,0.0)/_infectiousDuration_2;
			a = _shape_param[2];
			b = _shape_param[3];
		}
		
		res = pseudo_beta(t, a, b);
	}
	

	
	
	// ============== SYPHILIS (Tp) ===============
	
	if (_name == Tp)
	{
		res = _minInfectiousness;
		
		// Latent initial period
		
		if (_latentDuration > stiDuration)
			res = 0.0;
		
		// Primary syphilis
		
		if (_latentDuration < stiDuration
			&& stiDuration < _latentDuration+_infectiousDuration)
		{
			double a1 = _shape_param[0];
			double b1 = _shape_param[1];
			double t = stiDuration - _latentDuration;
			double d1 = _infectiousDuration;
			double v1 = _maxInfectiousness_Tp1;
			
			res = v1 * pseudo_beta(t/d1, a1, b1);
			res = max(res,_minInfectiousness);
		}
		
		// Secondary syphilis
		
		// First, check if secondary syphilis develops
		if (_is_recurrent)
		{
			if (_latentDuration+_infectiousDuration < stiDuration
				&& stiDuration<_latentDuration_2){
				
				bool condylo = _is_symptomatic;
				
				// No condylomata
				if (!condylo){
					double a2 = _shape_param[2];
					double b2 = _shape_param[3];
					double t = stiDuration - (_latentDuration+_infectiousDuration);
					double d2 = _infectiousDuration_2;
					double v2 = _maxInfectiousness_Tp2;
					
					res = v2 * pseudo_beta(t/d2, a2, b2);
					res = max(res,_minInfectiousness);
				}
				
				// Condylomata develop
				if (condylo){
					double a2c = _shape_param[4];
					double b2c = _shape_param[5];
					double t = stiDuration - (_latentDuration+_infectiousDuration);
					double d2c = _infectiousDuration_2c;
					double v2c = 1.0;
					
					res = v2c * pseudo_beta(t/d2c, a2c, b2c);
					res = max(res,_minInfectiousness);
				}
				
			}
		}
		
		// Early latent syphilis.
		// First, check if early latent syphilis develop
		if (_is_recurrent_2){
			if (_latentDuration_3 < stiDuration){
				double a3 = _shape_param[6];
				double b3 = _shape_param[7];
				double t = stiDuration - _latentDuration_3;
				double d3 = _infectiousDuration_3;
				double v3 = _maxInfectiousness_TpEL;
				
				res = v3 * pseudo_beta(min(t/d3,1.0), a3, b3);
				// "res" is not floored because
				// not infectious after early latent psypihlis
			}
		}
	}
	
	
	// ============== Trichomoniasis ===============
	
	if (_name == Tv)
	{
		// Infectivity curve categorized
		// by gender and symptomatic status
	
		// _shape_param[0]: "a" for female asymptomatic
		// _shape_param[1]: "b" for female asymptomatic
		
		// _shape_param[2]: "a" for female symptomatic
		// _shape_param[3]: "b" for female symptomatic
		
		// _shape_param[4]: "a" for male asymptomatic
		// _shape_param[5]: "b" for male asymptomatic
		
		// _shape_param[6]: "a" for male symptomatic
		// _shape_param[7]: "b" for male symptomatic
		
		
		// "_infectiousDuration" not used, but other similar and more details parameters
		
		double L = _latentDuration;
		double D =-9E9,a=-9E9,b=-9E9,rho=-9E9;
		
		if (g==female){
			if (!_is_symptomatic){
				 D = _infectiousDuration_f_a;
				 a = _shape_param[0];
				 b = _shape_param[1];
				rho = _asymptom_IC_reduc;
			}

			if (_is_symptomatic){
				 D = _infectiousDuration_f_s;
				 a = _shape_param[2];
				 b = _shape_param[3];
				rho = 1.0;
			}
		}
		
		if (g==male){
			if (!_is_symptomatic){
				 D = _infectiousDuration_m_a;
				 a = _shape_param[4];
				 b = _shape_param[5];
				rho = _asymptom_IC_reduc;
			}
			
			if (_is_symptomatic){
				 D = _infectiousDuration_m_s;
				 a = _shape_param[6];
				 b = _shape_param[7];
				rho = 1.0;
			}
		}

		double t = max(stiDuration -L,0.0)/D;
				
		res = rho * pseudo_beta(t, a, b);
	}
	
	
	// ============== BV ===============
	
	if (_name == Bv){
		// Not implemented for now
		res = 0;
	}
	
	// =====
	
	return res;	
}




void STI::check_infectivityCurve(string nameFile)
{
	// writes "x" and "IC(x)" to a file
	// (can be plotted outside C++)
	
	ofstream f(nameFile.c_str());
	
	
	// Default time limit
	double lim = 3.0;

	double timestep = 1.0/365.0;
	
	if (_name == HIV)	lim = 10.0;
	if (_name == HSV2)	lim = 3.0;
	if (_name == HPV)
	{
		lim = 3.0;
		_latentDuration = 0.4;
		_infectiousDuration = 1.0;
	}
	
	_is_symptomatic = true;
	_is_recurrent = true;
	
	Gender g= female;
	
	for (double x=0; x<lim; x+=timestep)
	{
		f << x << "," << infectivityCurve(x,g) << endl;
	}
}


void STI::displayInfo()
{
	coutline(80);
	cout<<"\t STI INFORMATION";
	coutline(80);
	
	cout << "Name:  " << STInameString(_name) << endl;
	cout << "Type:  " << _type << endl;
	
	cout << "latent Duration:"		<< _latentDuration << " ("<< _latentDuration*365.0 <<" days)"<<endl;
	cout << "Infectious Duration:  "	<< _infectiousDuration << " ("<< _infectiousDuration*365.0 <<" days)" << endl;
	cout << "Natural Clearance Duration:  " << _naturalClearanceDuration << endl;
	
	cout << "Proba Recurrence:  "		<< _proba_recurrence << endl;
	cout << "Is Recurrent:  "			<< _is_recurrent << endl;
	cout << "Is Recurrent2:  "			<< _is_recurrent_2 << endl;
	cout << "Recurrence_duration:  "	<< _recurrence_duration << endl;
	cout << "Recurrence_freq:  "		<< _recurrence_freq << endl;
	cout << "Min Infectiousness:  "		<< _minInfectiousness << endl;
	
	cout << "Proba_symptomatic_female:  "	<< _proba_symptomatic_female << endl;
	cout << "Proba_symptomatic_male:  "		<< _proba_symptomatic_male << endl;
	cout << "Is symptomatic:  "			<< _is_symptomatic << endl;
	
	cout << "Proba treatment failure:  "	<< _proba_treatmentFailure<<endl;
	cout << "Optimal treatment duration:  "<< _optimalTreatmentDuration<<endl;
	
	cout << "Proba Reduction SAT:"; displayVector(_probaSexTransm_SAT);
	
	if(_name==HIV)
	{
		cout << "HIV specific parameters:";
		displayVector(_HIVparam);
	}
	
	coutline(80);
}


double STI::TREstar(double treatmentDuration)
{
	// TO DO: For now same shape for all STIs
	//  BUT MUST IMPLEMENT SPECIFIC CASES
	
	double od = _optimalTreatmentDuration;
	
	if (od<treatmentDuration)
		return 0;

	return (pow((od-treatmentDuration), 2)/od/od);
}




// ================================================================================
// ================================================================================
// ================================================================================
// ========      OUTSIDE CLASS      ===============================================
// ================================================================================
// ================================================================================
// ================================================================================


string STInameString(STIname sti)
{
	switch (sti) 
	{
		case HIV:
			return "HIV";
			break;
		case HSV2:
			return "HSV2";
			break;
		case HPV:
			return "HPV";
			break;
		case Ct:
			return "Ct";
			break;
		case Ng:
			return "Ng";
			break;
		case Tp:
			return "Tp";
			break;
		case Hd:
			return "Hd";
			break;
		case Bv:
			return "Bv";
			break;
		case Tv:
			return "Tv";
			break;	

		case Cv:
			return "Cv";
			break;

		default:
			cout << "ERROR [STInameString]: STI name "<< sti << " does not exist"<<endl;
			exit(1);
			break;
	}
}



STIname StringToSTIname(string name)
{
	/// CONVERT A 'string' TO A 'STIname'
	
	bool found = false;
	STIname s;
	
	if (name=="HIV") {s=HIV ; found=true;}
	if (name=="HSV2") {s=HSV2 ; found=true;}
	if (name=="HPV") {s=HPV ; found=true;}
	if (name=="Ct") {s=Ct ; found=true;}
	if (name=="Ng") {s=Ng ; found=true;}
	if (name=="Tp") {s=Tp ; found=true;}
	if (name=="Hd") {s=Hd ; found=true;}
	if (name=="Bv") {s=Bv ; found=true;}
	if (name=="Tv") {s=Tv ; found=true;}
	
	
	string errmsg = "ERROR [StringToSTIname]: STI name does not exist: '" + name+"'";
	stopif(!found,errmsg);

	return s;
}

vector<STIname> StringToSTIname_vector(vector<string> s)
{
	vector<STIname> res;
	
	for (int i=0; i<s.size(); i++) {
		res.push_back(StringToSTIname(s[i]));
	}
	
	return res;
}


void check_ALL_infectivityCurves()
{
	// TO DO: replace hard coded value cleanly
	int nSTI = 9;
	vector<STIname> stiName(nSTI);
	
	stiName[0]=HIV;
	stiName[1]=HSV2;
	stiName[2]=HPV;
	stiName[3]=Ct;
	stiName[4]=Ng;
	stiName[5]=Tp;
	stiName[6]=Hd;
	stiName[7]=Bv;
	stiName[8]=Tv;
	
	vector<STI> STItemplate(nSTI);
	
	for (int s=0; s<nSTI; s++) 
	{
		STI tmp(stiName[s], "in_STI.csv");
		STItemplate[s] = tmp;
	}
	
	STI sHIV = STItemplate[0];
	sHIV.check_infectivityCurve("checkIC_HIV.out");
	
	STI sHSV2 = STItemplate[1];
	sHSV2.check_infectivityCurve("checkIC_HSV2.out");
	
	
	STI sCt = STItemplate[3];
	sCt.check_infectivityCurve("checkIC_Ct.out");
	
}


int positionSTIinVector(STIname s, vector<STI> v)
{
	/// FOR A GIVEN VECTOR OF STIs
	/// RETURNS THE POSITION OF THE ELEMENT IN THIS VECTOR
	/// WHOSE NAME IS 's'
	
	int i=0;
	for (i=0; v[i].get_name()!=s && i<v.size(); i++) {
		// let i being incremented in the for loop
	}
	
	string errmsg = "STI name <" + STInameString(s) + "> does not exist!" ;
	stopif(i>=v.size(),errmsg);

	return i;
}



double IC(string fct_family, double t, vector<double> param)
{
	/// Return the value of the infectivity curve
	/// at time t, given a functional family and its parameters
	
	string errmsg = "Functional form name '"+ fct_family+"' unknown for infectivity curve!";
	stopif( ((fct_family!="pseudo_beta") &&
			(fct_family!="pseudo_gamma")),
		   errmsg);
	
	double res=-9999.99;
	
	if (fct_family=="pseudo_beta")
	{
		stopif(param.size()!=2, "Number of parameters for pseudo beta is not two");
		res = pseudo_beta(t, param[0], param[1]);
	}
	
	if (fct_family=="pseudo_gamma")
	{
		stopif(param.size()!=2, "Number of parameters for pseudo gamma is not two");
		res = pseudo_gamma(t, param[0], param[1]);
	}
	
	return res;
}

