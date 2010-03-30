/*
 * Random_EM_Map.cpp
 *
 *  Created on: Mar 4, 2010
 *      Author: mckinsel
 */

#include "Random_EM_Map.h"
#include <assert.h>

Random_EM_Map::Random_EM_Map () {
	d = 1;
	isoform_length = 1;
}
Random_EM_Map::Random_EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, MarkovChain& mc) {

	assert(bt1.base_read_id == bt2.base_read_id);

	base_read_id = bt1.base_read_id;
	isoform = "random";

	P = 1;

	P *= mc.sequence_probability(bt1.read->sequence);
	P *= mc.sequence_probability(bt2.read->sequence);

}

long double Random_EM_Map::em_prob(long double theta_i) {
	return theta_i*P;
}
Random_EM_Map::~Random_EM_Map() {
	// TODO Auto-generated destructor stub
}
