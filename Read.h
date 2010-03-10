/*
 * Read.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef READ_H_
#define READ_H_

#include <assert.h>
#include <fstream>
#include <string>
#include "Quality.h"
using namespace std;

class Read {
public:
	Read(int offset);
	Read(string name, string seq, string qual, int offset);
	virtual ~Read();

	string id;
	string sequence;
	string base_id;
	Quality * quality;
	int offset;

	friend istream &operator>>(istream &stream, Read &rd);

};

#endif /* READ_H_ */
