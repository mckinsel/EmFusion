/*
 * MarkovChain.cpp
 *
 *  Created on: Mar 4, 2010
 *      Author: mckinsel
 */

#include <assert.h>
#include "MarkovChain.h"
#include <iostream>

using namespace std;

MarkovChain::MarkovChain(unsigned int ord, int pseudo) {

	order = ord;
	pseudocount = pseudo;
	total_count = 0;

}

MarkovChain::~MarkovChain() {
	// TODO Auto-generated destructor stub
}

void MarkovChain::add_sequence(string seqtoadd) {
//	cout << seqtoadd.length() << " " << order << endl;
//	cout << seqtoadd.length() - order << " " << seqtoadd << endl;
	for(unsigned int i=0; i<seqtoadd.length() - order; i++){
		string_counts[seqtoadd.substr(i, order)]++;
		total_count++;
	}
}

long double MarkovChain::ordermer_prob(string ordermer) {
	assert(ordermer.length() == order);
	return (string_counts[ordermer] + pseudocount)/(long double)total_count;
}

long double MarkovChain::sequence_probability(string evalstring) {

	long double prob = 1;

	for(unsigned int i=0; i<evalstring.length()-order; i++){
		prob *= ordermer_prob(evalstring.substr(i, order));
	}
	return prob;
}
