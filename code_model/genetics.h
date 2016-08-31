//
//  genetics.h
//  stiagent--tmp
//
//  Created by David CHAMPREDON on 2016-08-30.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef __stiagent__tmp__genetics__
#define __stiagent__tmp__genetics__

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <random>


using namespace std;

string  create_genome(unsigned long length);

unsigned long distance_genomes(string a, string b);

char	random_mutation_gene(char original_gene);

string	mutate_genome(string genome, unsigned long n_genes);

#endif /* defined(__stiagent__tmp__genetics__) */
