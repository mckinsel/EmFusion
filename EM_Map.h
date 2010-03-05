/*
 * EM_Map.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef EM_MAP_H_
#define EM_MAP_H_

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include "BowtieEntry.h"

using namespace std;

typedef tr1::unordered_map<int, double> int2doubleumap;
typedef tr1::unordered_map<string, int> string2intumap;


class EM_Map {
public:
	EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, int2doubleumap&, string2intumap&);
	virtual ~EM_Map();

	long double em_prob(long double);

	string base_read_id;
	string isoform;
	int start;
	int end;
	double d;
	int isoform_length;
	long double P;

protected:
	EM_Map();
};

#endif /* EM_MAP_H_ */
