/*
 * EM_Map.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "EM_Map.h"
#include <assert.h>

EM_Map::EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, int2doubleumap& d_prob, string2intumap& i_length) {

	assert(bt1.base_read_id == bt2.base_read_id);
	assert(bt1.mapping == bt2.mapping);

	base_read_id = bt1.base_read_id;
	isoform = bt1.mapping;
	int start = min(bt1.position, bt2.position);
	int end = max(bt1.position, bt2.position);
	isoform_length = i_length[isoform];
//	cout << "start " << start << " end " << end << " absdiff " << abs(start - end) << " d " << d_prob[abs(start - end)] <<endl;
	d = d_prob[abs(start - end)];
	P = bt1.mapping_probability() * bt2.mapping_probability();


}

EM_Map::EM_Map() {
//	start = -1;
//	end = -1;
	d = 1;
	isoform_length = 1;
}
EM_Map::~EM_Map() {
	// TODO Auto-generated destructor stub
}

long double EM_Map::em_prob(long double theta_i) {
//	cout << "theta_i " << theta_i << endl;
//	cout << "1/(long double)isoform_length " << (1/(long double)isoform_length) << endl;
//	cout << "P " << P << endl;
//	cout << "d " << d << endl;

	return (1/(long double)isoform_length)*P*d*theta_i;
}
