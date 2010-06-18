/*
 * BowtieEntry.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "BowtieEntry.h"


BowtieEntry::BowtieEntry(int offs) {

	offset = offs;
	_mapped_transcript = "";
	_mapped_gene = "";
	read = NULL;
	position = -1;
	read_id = "";
	base_read_id = "";


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

	_mapped_transcript = "";
	_mapped_gene = "";

}

BowtieEntry::~BowtieEntry() {
	delete read;
}

istream& operator>>(istream &stream, BowtieEntry& bt) {
	char btline[1500];
	vector<string> lineparts;
	string mms;

	stream.getline(btline, (streamsize)1500);
	myTokenize(string(btline), lineparts, "\t");

	if (lineparts.size() == 8){
		mms = lineparts.at(7);
	} else {
		mms = "";
	}
	if(lineparts.size()>6 && !stream.eof()){
	bt.read_id = lineparts[0];
	bt.base_read_id = bt.read_id.substr(0, bt.read_id.size() - 2);
	bt.mapping = lineparts[2];
	bt.strand = lineparts[1];
	bt.position = (int)atoi(lineparts.at(3).c_str());
	delete bt.read;
	bt.read = new Read(lineparts[0], lineparts[4], lineparts[5], bt.offset);
	bt.mismatch_indices = bt.get_indices(mms);
	}
  return stream;
}

vector<int> BowtieEntry::get_indices(string mismatchstring) {
	vector<string> parts;
	vector<int> output;
	vector<int> corrected_output;
	size_t coloni;

	myTokenize(mismatchstring, parts, ",");

	for(unsigned int i=0; i<parts.size(); i++){
		coloni = parts[i].find_first_of(":");
		parts[i].erase(coloni, parts[i].length()-coloni);
		output.push_back((int)atoi(parts.at(i).c_str()));
	}

	if(strand == "+"){
		corrected_output = output;
	} else if(strand == "-"){
		int read_length = read->sequence.length();
		for(unsigned int i=0; i<output.size(); i++){
			corrected_output.push_back((-1*(output.at(i) - read_length + 1)));
		}
	}

	return corrected_output;
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

string BowtieEntry::mapped_gene() {
	if(_mapped_gene.length() > 0) {return _mapped_gene;}
	else {parse_mapping(); return _mapped_gene;}
}
string BowtieEntry::mapped_transcript() {
	if(_mapped_transcript.length() > 0) {return _mapped_transcript;}
	else {parse_mapping(); return _mapped_transcript;}
}

void BowtieEntry::parse_mapping() {

	char * maptok;
	char * cstr;

	cstr = new char [mapping.size() + 1];
	strcpy(cstr, mapping.c_str());

	maptok = strtok(cstr, "|");
	_mapped_transcript = string(maptok);

	maptok = strtok(NULL, "|");
	_mapped_gene = string(maptok);

	delete[] cstr;


}
