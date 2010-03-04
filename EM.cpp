/*
 * EM.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include <fstream>
#include <tr1/unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include "EM_Map.h"
#include "BowtieEntry.h"
#include "Utils.h"

using namespace std;
using namespace __gnu_cxx;
using namespace tr1;



typedef tr1::unordered_map<string, vector<EM_Map*> > vectorumap;
typedef tr1::unordered_map<string, long double > longdoubleumap;


void EM_Update( vectorumap & read2emmap, vectorumap & isoform2emmap, longdoubleumap & th, longdoubleumap & newth, int N){

	longdoubleumap read_sums;

	vectorumap::iterator read_iterator;
	vectorumap::iterator isoform_iterator;

	int counter = 0;
	for(read_iterator=read2emmap.begin(); read_iterator != read2emmap.end(); read_iterator++){

		long double sumread = 0;

		string read_id = read_iterator->first;
		vector<EM_Map*> emmaps = read_iterator->second;
		for(unsigned int i=0; i < emmaps.size(); i++){
			sumread += emmaps.at(i)->em_prob(th[emmaps.at(i)->isoform]);
		}
		read_sums[read_id] += sumread;
		if(sumread == 0){
			cout << "Read " << read_id << " sumread is " << sumread << endl;
		}
		counter++;
		if(counter % 500000 == 0) cout << counter << " reads processed." << endl;

	}

	counter = 0;
	for(isoform_iterator=isoform2emmap.begin(); isoform_iterator != isoform2emmap.end(); isoform_iterator++){

		long double sumterm = 0;
		string isoform_id = isoform_iterator->first;
		vector<EM_Map*> isoemmaps = isoform_iterator->second;
		for(unsigned int i=0; i < isoemmaps.size(); i++){
			sumterm += isoemmaps.at(i)->em_prob(th[isoform_id])/read_sums[isoemmaps.at(i)->base_read_id];
		}
		newth[isoform_id] = sumterm/N;
		counter++;
		if(counter % 50000 == 0) cout << counter << " isoforms processed." << endl;

	}

}



int main(int argc, char * argv[]){

	char * btfilename = argv[1];
	int offset = atoi(argv[2]);
	ifstream btstream(btfilename);
	BowtieEntry bt1(offset);
	BowtieEntry bt2(offset); //Initialize with offset

	vectorumap read_to_emmaps;
	vectorumap isoform_to_emmaps;

	vectorumap::iterator isoform_iterator;
	vectorumap::iterator read_iterator;

	longdoubleumap theta;
	longdoubleumap newtheta;
	int counter = 0;



	while(btstream >> bt1 && !btstream.eof()){
		btstream >> bt2;

		EM_Map * pemmap;
		pemmap = new EM_Map(bt1, bt2, 10000); ///BADDDDDDDD, replace 10000

		read_to_emmaps[pemmap->base_read_id].push_back(pemmap);
		isoform_to_emmaps[pemmap->isoform].push_back(pemmap);

		delete bt1.read;
		delete bt2.read;
		counter++;
		if(counter%1000000==0) cout << counter << endl;
	}

	cout << "Bowtie reads read in." << endl;

	int N = (int)read_to_emmaps.size();
	int M = (int)isoform_to_emmaps.size();
	cout << "N " << N << " M " << M << endl;

	isoform_iterator = isoform_to_emmaps.begin();
	while(isoform_iterator != isoform_to_emmaps.end()){
		theta[isoform_iterator->first] = 1/(long double)M;
		++isoform_iterator;
	}

	cout << "Made first theta guess." << endl;
	int count = 0;

	while(count < 10){
		cout << "Starting iteration " << count << endl;
		EM_Update(read_to_emmaps, isoform_to_emmaps, theta, newtheta, N);
		theta = newtheta;
		count++;
	}

	ofstream outstream("em_output");

	for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
		outstream << thiter->first << "\t" << thiter->second << endl;
	}
	return 0;

}
