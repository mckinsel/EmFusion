/*
 * Read.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef READ_H_
#define READ_H_

#include <string>
#include "Quality.h"
using namespace std;

class Read {
public:
	Read(string name, string seq, string qual, int offset);
	virtual ~Read();

	string id;
	string sequence;
	Quality * quality;

};

#endif /* READ_H_ */
