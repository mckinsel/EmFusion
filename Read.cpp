/*
 * Read.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "Read.h"

Read::Read(string name, string seq, string qual, int offset) {
	id = name;
	sequence = seq;
	quality = new Quality(qual, offset);
}

Read::~Read() {
	delete quality;
}
