/*
 * FastaEntry.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: marcus
 */

#include "FastaEntry.h"

FastaEntry::FastaEntry(){
	id = "";
	sequence = "";
}
FastaEntry::FastaEntry(string name, string seq) {

	id = name;
	sequence = seq;

}

FastaEntry::~FastaEntry() {
	// TODO Auto-generated destructor stub
}

istream& operator>>(istream &stream, FastaEntry& fa) {
	char faline[1500];

	stream.getline(faline, (streamsize)1500);
	while(faline[0] != '>' && !stream.eof()){
		stream.getline(faline, (streamsize)1500);
	}

	if(stream.eof()) return stream;

	assert(faline[0] == '>');

	string id = string(faline).substr(1, string(faline).length() - 1);

	stream.getline(faline, (streamsize) 1500);

	string seq = string(faline);

	fa.id = id;
	fa.sequence = seq;

  return stream;
}
