/*
 * EM_Map.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "EM_Map.h"
#include <assert.h>

EM_Map::EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, int i_length) {

	assert(bt1.base_read_id == bt2.base_read_id);
	assert(bt1.mapping == bt2.mapping);

	base_read_id = bt1.base_read_id;
	isoform = bt1.mapping;
	start = min(bt1.position, bt2.position);
	end = max(bt1.position, bt2.position);
	isoform_length = i_length;
	P = bt1.mapping_probability() * bt2.mapping_probability();


}

EM_Map::~EM_Map() {
	// TODO Auto-generated destructor stub
}

long double EM_Map::em_prob(long double theta_i) {
	return (1/(long double)isoform_length)*P*theta_i*d;
}
