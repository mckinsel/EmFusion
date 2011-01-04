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
#include <iomanip>
#include <tr1/unordered_map>
#include "BowtieEntry.h"

using namespace std;

typedef tr1::unordered_map<int, double> int2doubleumap;
typedef tr1::unordered_map<string, int> string2intumap;


class EM_Map {
public:
	EM_Map();
	EM_Map(BowtieEntry& bt1, BowtieEntry& bt2, int2doubleumap&, string2intumap&);
	virtual ~EM_Map();

	long double em_prob(long double);

	string base_read_id;
	string isoform;
//	int start;
//	int end;
	double d;
	int isoform_length;
	long double P;
	double r;
	int start;
	int end;

	friend ostream &operator<<(ostream &stream, const EM_Map &emm);
	friend istream &operator>>(istream &stream, EM_Map &emm);

};

#endif /* EM_MAP_H_ */
