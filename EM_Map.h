/*
 * EM_Map.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef EM_MAP_H_
#define EM_MAP_H_

#include <string>
#include "BowtieEntry.h"

using namespace std;

class EM_Map {
public:
	EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, int);
	virtual ~EM_Map();

	long double em_prob(long double);

	string base_read_id;
	string isoform;
	int start;
	int end;
	double d;
	int isoform_length;
	long double P;
};

#endif /* EM_MAP_H_ */
