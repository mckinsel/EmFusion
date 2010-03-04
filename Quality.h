/*
 * Quality.h
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#ifndef QUALITY_H_
#define QUALITY_H_
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;

class Quality {
private:
	void calculate_error_probabilites();
	double phred_probability(int phredscore);

public:
	Quality(string , int);
	virtual ~Quality();

	string quality_str;
	int offset;
	vector<double> error_probabilities;

};

#endif /* QUALITY_H_ */
