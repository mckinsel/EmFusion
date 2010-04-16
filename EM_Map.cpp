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

EM_Map::EM_Map () {
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

	if(isoform.substr(0, 2) == "F_"){
		r = 1;
	} else {
		r = 1;
	}
	return (1/(long double)isoform_length)*P*d*theta_i*r;
}

ostream& operator<<(ostream &stream, const EM_Map &emm) {
	stream << emm.base_read_id << "\t" << emm.isoform << "\t" << emm.isoform_length << "\t" << setprecision(25) << emm.d << "\t" << emm.P << endl;
	return stream;
}

istream& operator>>(istream &stream, EM_Map &emm) {
	stream >> emm.base_read_id;
	stream >> emm.isoform;
	stream >> emm.isoform_length;
	stream >> emm.d;
	stream >> emm.P;
//	cout << "read in a EM_Map " << emm.base_read_id << " " << emm.isoform << " " << emm.d << " " << emm.P << " " << emm.isoform_length << " " << emm.em_prob(1) << endl;
	return stream;
}
