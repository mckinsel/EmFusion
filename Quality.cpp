/*
 * Quality.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include "Quality.h"

Quality::Quality(string q_str, int offs) {
	quality_str = q_str;
	offset = offs;
	calculate_error_probabilites();
}

Quality::~Quality() {
	// TODO Auto-generated destructor stub
}

void Quality::calculate_error_probabilites(){

	for(unsigned int i=0; i<quality_str.length(); i++){
//		cout << "conv to int " << (int)quality_str.at(i) << " " << quality_str.at(i) << endl;
		int qualint = (int)quality_str.at(i);
//		cout << "qualint " << qualint << " minus offset " << qualint - offset << endl;
		error_probabilities.push_back(phred_probability(qualint - offset));
	}

}

double Quality::phred_probability(int phredscore){
//	cout << "phredscore " << phredscore << " " << phredscore/(long double)-10 << endl;
	return pow(10, phredscore/(long double)-10);
}
