//
//  genetics.cpp
//  stiagent--tmp
//
//  Created by David CHAMPREDON on 2016-08-30.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include "genetics.h"
#include "globalVar.h"
#include "dcTools.h"


string create_genome(unsigned long length){
	/// CREATE A RANDOM GENOME
	
	string g(length,'!');
	
	std::uniform_int_distribution<uint> unif(1,4);
	
	for (unsigned long i=0; i<length; i++) {
		uint x = unif(_RANDOM_GENERATOR);
		if (x==1)      g[i] = 'A';
		else if (x==2) g[i] = 'C';
		else if (x==3) g[i] = 'G';
		else if (x==4) g[i] = 'T';
	}
	return g;
}


char random_mutation_gene(char old_gene){
	
	char mutated_gene = '!';
	
	std::uniform_int_distribution<uint> unif(1,3);
	uint x = unif(_RANDOM_GENERATOR);
	
	if (old_gene == 'A') {
		if (x==1)      mutated_gene = 'C';
		else if (x==2) mutated_gene = 'G';
		else if (x==3) mutated_gene = 'T';
	}
	
	else if (old_gene == 'C') {
		if (x==1)      mutated_gene = 'A';
		else if (x==2) mutated_gene = 'G';
		else if (x==3) mutated_gene = 'T';
	}
	
	else if (old_gene == 'G') {
		if (x==1)      mutated_gene = 'A';
		else if (x==2) mutated_gene = 'C';
		else if (x==3) mutated_gene = 'T';
	}
	
	else if (old_gene == 'T') {
		if (x==1)      mutated_gene = 'A';
		else if (x==2) mutated_gene = 'G';
		else if (x==3) mutated_gene = 'C';
	}
	return mutated_gene;
}



string	mutate_genome(string genome, unsigned long n_genes){
	/// RANDOMLY MUTATE GENES IN A GENOME
    
    string mutated_genome = genome;
	unsigned long gl = genome.length();
    
    stopif(gl==0, "Trying to mutate an empty genome.");
    
	std::uniform_int_distribution<unsigned long> unif(0, gl-1);
	
	for (uint i=0; i<n_genes; i++) {
		unsigned long idx = unif(_RANDOM_GENERATOR);
		char old_gene = genome.at(idx);
		char mutated_gene = random_mutation_gene(old_gene);
		mutated_genome[idx] = mutated_gene;
	}
    return mutated_genome;
}



unsigned long distance_genomes(string a, string b){
    
    unsigned long na = a.length();
    unsigned long nb = b.length();
    
    unsigned long n = (na<nb)?na:nb;
    
    unsigned long dist = 0;
    for (unsigned long i=0; i<n; i++) {
        if (a[i] != b[i]) dist++;
    }
    return dist;
}
