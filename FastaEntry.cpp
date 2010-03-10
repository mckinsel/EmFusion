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
	assert(faline[0] == '>');

	string id = string(faline);

	stream.getline(faline, (streamsize) 1500);

	string seq = string(faline);

	fa.id = id;
	fa.sequence = seq;

  return stream;
}
