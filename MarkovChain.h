/*
 * MarkovChain.h
 *
 *  Created on: Mar 4, 2010
 *      Author: mckinsel
 */

#ifndef MARKOVCHAIN_H_
#define MARKOVCHAIN_H_

#include <string>
#include <tr1/unordered_map>



using namespace std;

typedef tr1::unordered_map<string, long int> umap;

class MarkovChain {
private:
	umap string_counts;
	int total_count;
	int order;
	int pseudocount;
	long double ordermer_prob(string);

public:
	MarkovChain(int, int);
	virtual ~MarkovChain();
	long double sequence_probability(string);
	void add_sequence(string);

};

#endif /* MARKOVCHAIN_H_ */
