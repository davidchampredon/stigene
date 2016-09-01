//
//  main.cpp
//  Epidemic_Models
//
//  Created by David Champredon
//  Copyright (c) 2012-14. All rights reserved.
//

//#include <sys/time.h>

#include "globalVar.h"

#include "dcTools.h"
#include "dcMatrix.h"
#include "RV.h"

#include "population.h"
#include "simulation.h"
#include "MCsimulation.h"
#include "compare_simulation.h"

#include <random>

#include "genetics.h"


int main(int argc, const char * argv[])
{
	
	bool doObj = true;

	system("date");
	system("pwd");
	
	
	// For performance monitoring
	// - do not delete -
	timeval tim;
	gettimeofday(&tim, NULL);
	double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	// ------------------------------------------
	
	string _DIR_IN = "../inputs/";
	
	// Switches to determine what will be done
	// for this program execution
	string main_switches_file = _DIR_IN + "main_switches.csv";
	
	bool doTest			= 0; //(bool)(getParameterFromFile("doTest", main_switches_file));
	bool doSingleRun	= 1; //(bool)(getParameterFromFile("doSingleRun", main_switches_file));
	
	
	string calib_files = "calibration_filename_wrapper_ZM4.csv";
	
	coutline(80);
	cout<<"doTest= "		<< doTest<<endl;
	cout<<"doSingleRun= "	<< doSingleRun<<endl;
	coutline(80);
	
	if (!doTest)
	{
		// ======================
		// === Initialization ===
		// ======================
		
		// Initialize empty Population object
		Population P(0);
		Population P2(0);
		
		bool debugInfo=true;
		
		unsigned long founder_size	= 250;
		double founder_femprop		= 0.5;
		double founder_cswprop		= 0.01;
		string folder_inputs		= "../inputs/";
		
		P.setup_for_simulation(founder_size,
							   founder_femprop,
							   founder_cswprop,
							   folder_inputs,
							   "in_STI.csv",
							   "in_STI_SFincrease.csv",
							   "in_HIVrebound.csv",
							   "in_STItreatment.csv",
							   "in_STI_vaccine.csv",
                               "in_STI_genomes.csv",
							   debugInfo);
		
		cout <<endl<< "  Populations setup done."<<endl;
//		P.displayInfo(true);
//		exit(0);
		
		// ======================
		// === Run simulation ===
		// ======================	
		
        double horizon	= 5.0; getParameterFromFile("horizon_years", _DIR_IN + "in_simulation.csv");
		double timeStep	= 2.0/365.0; //getParameterFromFile("timestep_days", _DIR_IN + "in_simulation.csv")/365.0;
		
		
		if (doSingleRun)
		{
			double			horizon_prtn		= 13.0 ; //50.0;
			double			timestep_prtn		= 20.0/365.0;
			bool			TraceNetwork		= false;
			unsigned int	iter_mc				= 1;
			int				displayProgress		= 11;
			
			string file_init_STI	= _DIR_IN + "in_STI_initial_prevalence.csv";
			
			vector<string> file_intervention;
			string file_interv_base =_DIR_IN + "in_scenario_VaxMass.csv";  //in_scenario_baseline.csv
			vectorFromCSVfile_string(file_intervention,file_interv_base.c_str(), 1);
			displayVector(file_intervention);
			file_intervention = trim(file_intervention);
			
			Simulation S;
			
			
			if(doObj){
				string folder_inputs = _DIR_IN;
				string folder_calib  = _DIR_CALIB;
				Simulation Sobj = runSimulation_one_obj(P,
														file_init_STI,
														file_intervention,
														horizon_prtn,
														timestep_prtn,
														horizon,
														timeStep,
														TraceNetwork,
														displayProgress,
														iter_mc,
														folder_inputs,
														folder_calib
														);
				
				dcDataFrame df = Sobj.get_df_sim();
				dcDataFrame export_pop = Sobj.get_population().export_to_dataframe();
				//df.display();
                
                Population POP = Sobj.get_population();
                STIname sn = POP.get_STI()[0].get_name();
                
                vector<string> g = POP.get_population_genomes(sn);
                vector<unsigned long> g_uid = POP.census_STIinfected_UID(sn);
                displayVector(g);
                displayVector(g_uid);
                cout << " dist = " << distance_genomes(g[0], g[1]);
                
                vector< vector<unsigned long> > z;
                STIname stiname = POP.get_STI()[0].get_name();
                for(unsigned long i=0; i < POP.get_size(); i++){
                    vector<unsigned long> seccases_uid = POP.get_individual(i).get_STI_secondary_cases(stiname);
                    unsigned long tmp = POP.get_individual(i).get_UID();
                    z.push_back(seccases_uid);
                }
                displayVector(z);
			}
		}
	} // ---- end of !doTest -----
	
	
	
	// --- TESTS ----------------------------
	
	if (doTest)
	{
        unsigned int	_RANDOM_SEED	= 1234;
        std::mt19937_64	_RANDOM_GENERATOR1(_RANDOM_SEED);
        
        unsigned long N = 3;
        double p = 0.02;
        unsigned long iter = 10000;
        
        // Distribution typed 'long'
        
        std::binomial_distribution<long> d(N, p);
        for(int i = 0; i < iter; i++) {
            unsigned long r = d(_RANDOM_GENERATOR1);
         cout << i << " " << r <<endl;
            if(r > N) {cout<<"PROBLEM"<<endl; exit(99);} // <-- never reached (which is good!)
        }
        
        // Distribution typed 'unsigned long'
        std::mt19937_64     _RANDOM_GENERATOR2(_RANDOM_SEED);
        
        std::binomial_distribution<int> d2(N, p);
        for(int i = 0; i < iter; i++) {
            unsigned long r = d2(_RANDOM_GENERATOR2);
            cout << i << " " << r <<endl;
            if(r > N) {cout<<"PROBLEM"<<endl; exit(99);} // <-- eventually reached, with r = huge number
        }

        
        
//		string genome = create_genome(1000);
//        string genome1 = genome;
//		cout << genome <<endl;
//		mutate_genes(genome, 12);
//		cout << genome <<endl;
//		
//        cout << "distance: " << distance_genomes(genome1, genome)<<endl;
//		
//        
//        string g = create_genome(3000);
//        cout<<g<<endl;
        
        
	}
	
	// ---- END TESTS -----------------------
	
	
	
	
	// --------------------------------------------------------------
	// COMPUTER TIME MONITORING - do not delete!
	
	gettimeofday(&tim, NULL);
	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	int minutes = (int)((t2-t1)/60.0);
	double sec = (t2-t1)-minutes*60.0;
	cout << endl << " - - - Computational time : ";
	cout << minutes<<" min "<<sec<<" sec" << endl;
	
	// --------------------------------------------------------------
	
	return 0;
}
