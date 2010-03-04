/*
 * BowtieEntry.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "BowtieEntry.h"


BowtieEntry::BowtieEntry(int offs) {

	offset = offs;
}
BowtieEntry::BowtieEntry(string read_name, string map, string orient, int pos, string seq,
						 string qual, string mismatches, int off) {

	read_id = read_name;
	base_read_id = read_id.substr(0, read_id.size() - 2);
	mapping = map;
	strand = orient;
	position = pos;
	offset = off;
	read = new Read(read_id, seq, qual, offset);
	mismatch_indices = get_indices(mismatches);

}

BowtieEntry::~BowtieEntry() {
	// TODO Auto-generated destructor stub
}

istream& operator>>(istream &stream, BowtieEntry& bt) {
	char btline[1500];
	vector<string> lineparts;
	string mms;
//	cout << stream.eof() << endl;

	stream.getline(btline, (streamsize)1500);
//	cout << "btline " << btline << endl;
//	cout << stream.eof() << endl;
	myTokenize(string(btline), lineparts, "\t");

	if (lineparts.size() == 8){
		mms = lineparts.at(7);
	} else {
		mms = "";
	}
	if(lineparts.size()>6){
	bt.read_id = lineparts[0];
	bt.base_read_id = bt.read_id.substr(0, bt.read_id.size() - 2);
	bt.mapping = lineparts[2];
//	cout << "mapping " << bt.mapping <<endl;
	bt.strand = lineparts[1];
	bt.position = (int)atoi(lineparts.at(3).c_str());
	bt.read = new Read(lineparts[0], lineparts[4], lineparts[5], bt.offset);
	bt.mismatch_indices = bt.get_indices(mms);
	}
  return stream;
}

vector<int> BowtieEntry::get_indices(string mismatchstring) {
	vector<string> parts;
	vector<int> output;
	size_t coloni;

	myTokenize(mismatchstring, parts, ",");

	for(unsigned int i=0; i<parts.size(); i++){
		coloni = parts[i].find_first_of(":");
		parts[i].erase(coloni, parts[i].length()-coloni);
		output.push_back((int)atoi(parts.at(i).c_str()));
	}

	return output;
}

long double BowtieEntry::mapping_probability() {
	long double P;

	P = 1;

	for(unsigned int i=0; i<read->quality->quality_str.size(); i++){
		if(find(mismatch_indices.begin(), mismatch_indices.end(), i) == mismatch_indices.end()){ //if match
		//	cout << "error_prob at " << i << " " <<  read->quality->error_probabilities[i] <<endl;
			P *= (1 - read->quality->error_probabilities[i]);
		} else { //if mismatch
		//	cout << "error_prob at " << i << " " <<  read->quality->error_probabilities[i] <<endl;
			P *= read->quality->error_probabilities[i]/3;
		}

	}

	return P;

}
