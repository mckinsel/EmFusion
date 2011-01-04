/*
 * BowtieEntry.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef BOWTIEENTRY_H_
#define BOWTIEENTRY_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Read.h"
#include "Utils.h"
using namespace std;

class BowtieEntry {
private:
	vector<int> get_indices(string mismatchstring);
	string _mapped_gene;
	string _mapped_transcript;
	string _full_mapped_transcript;
	void parse_mapping();

public:
	BowtieEntry (int);
	BowtieEntry (string read_id, string mapping, string strand, int position, string seq,
				string qual, string mismatches, int offset);
	virtual ~BowtieEntry();

	string read_id;
	string base_read_id;
	string mapping;
	string strand;
	int position;
	Read * read;
	vector<int> mismatch_indices;
	int offset;

	friend istream &operator>>(istream &stream, BowtieEntry &bt);

	long double mapping_probability();

	string mapped_gene();
	string mapped_transcript();
	string full_mapped_transcript();

};

#endif /* BOWTIEENTRY_H_ */
