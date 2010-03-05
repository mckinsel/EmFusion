/*
 * Random_EM_Map.h
 *
 *  Created on: Mar 4, 2010
 *      Author: mckinsel
 */

#ifndef RANDOM_EM_MAP_H_
#define RANDOM_EM_MAP_H_

#include "EM_Map.h"
#include "MarkovChain.h"
#include "BowtieEntry.h"

class Random_EM_Map: public EM_Map {
private:
	string seq1, seq2;
public:
	Random_EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, MarkovChain& mc);
	virtual ~Random_EM_Map();
	long double em_prob(long double);
};

#endif /* RANDOM_EM_MAP_H_ */
