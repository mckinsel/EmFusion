/*
 * Read.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "Read.h"

Read::Read(int offs) {
	offset = offs;
}

Read::Read(string name, string seq, string qual, int offs) {
	id = name;
	base_id.substr(0, id.size() - 2);
	sequence = seq;
	quality = new Quality(qual, offs);
	offset = offs;
}

Read::~Read() {
	delete quality;
}

istream& operator>>(istream &stream, Read& rd) {

	char fastqline[1500];
	string id;
	string seq;
	string qual;

	stream.getline(fastqline, (streamsize)1500);
	assert(fastqline[0] == '@');
	id = string(fastqline);

	stream.getline(fastqline, (streamsize)1500);
	seq = string(fastqline);

	stream.getline(fastqline, (streamsize)2);
	stream.getline(fastqline, (streamsize)1500);
	qual = string(fastqline);

	delete rd.quality;

	rd.id = id;
	rd.sequence = seq;
	rd.quality = new Quality(qual, rd.offset);
	rd.base_id.substr(0, id.size() - 2);
	return stream;
}
