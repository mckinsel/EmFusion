/*
 * Read.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "Read.h"

Read::Read(int offs) {
	offset = offs;
	quality = NULL;
}

Read::Read(string name, string seq, string qual, int offs) {
	if(id[0] == '@') {
		id = name.substr(1, name.length() - 1);
	} else {
		id = name;
	}
	base_id = id.substr(0, id.size() - 2);
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
	while(fastqline[0] != '@'){
		stream.getline(fastqline, (streamsize)1500);
	}
	assert(fastqline[0] == '@');
	id = string(fastqline).substr(1, string(fastqline).length());

	stream.getline(fastqline, (streamsize)1500);
	seq = string(fastqline);

	stream.getline(fastqline, (streamsize)2);
	stream.getline(fastqline, (streamsize)1500);
	qual = string(fastqline);

	delete rd.quality;

	rd.id = id;
	rd.sequence = seq;
	rd.quality = new Quality(qual, rd.offset);
	rd.base_id = id.substr(0, id.size() - 2);
	return stream;
}

void Read::write_as_fasta(ostream & stream) {
	stream << ">" << id << endl;
	stream << sequence << endl;
}

void Read::write_as_fastq(ostream & stream) {
	stream << "@" << id << endl;
	stream << sequence << endl;
	stream << "+" << endl;
	stream << quality->quality_str << endl;
}
