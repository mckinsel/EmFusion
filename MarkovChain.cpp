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

	for(unsigned int i=0; i<seqtoadd.length() - order; i++){
		string_counts[seqtoadd.substr(i, order)]++;
		total_count++;
	}
}

long double MarkovChain::ordermer_prob(string ordermer) {
	if(ordermer.length() != order){
		cout << "ordermer " << ordermer << endl;
	}
	assert(ordermer.length() == order);
	return (string_counts[ordermer] + pseudocount)/((long double)total_count + pseudocount);
}

long double MarkovChain::sequence_probability(string evalstring) {

	long double prob = 1;
	if(evalstring.length() < 10){
	cout << "evalstring " << evalstring << " " << evalstring.length() << endl;
	}
	for(unsigned int i=0; i<evalstring.length() - order; i++){
		if(evalstring.length() < 10){
			cout << i << " " << evalstring.length() << evalstring.length() - order << endl;
		}
		prob *= ordermer_prob(evalstring.substr(i, order));
	}
	cout << "Evaluating MC probability of sequence " << evalstring << " and it's " << prob << endl;
	return prob;
}
