/*
 * FastaEntry.h
 *
 *  Created on: Mar 10, 2010
 *      Author: marcus
 */

#ifndef FASTAENTRY_H_
#define FASTAENTRY_H_

#include <assert.h>
#include <string>
#include <fstream>

using namespace std;

class FastaEntry {
public:
	FastaEntry();
	FastaEntry(string name, string seq);
	virtual ~FastaEntry();

	string id;
	string sequence;

	friend istream &operator>>(istream &stream, FastaEntry &fa);
};

#endif /* FASTAENTRY_H_ */
