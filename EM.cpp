/*
 * EM.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include <fstream>
#include <stdlib.h>
#include <tr1/unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include "EM_Map.h"
#include "BowtieEntry.h"
#include "Utils.h"

using namespace std;
using namespace __gnu_cxx;

void EM_Update( tr1::unordered_map<const char *, vector<EM_Map*> > & read2emmap,
				tr1::unordered_map<const char *, vector<EM_Map*> > & isoform2emmap,
				tr1::unordered_map<const char *, long double> & th,
				tr1::unordered_map<const char *, long double> & newth, int N){

	tr1::unordered_map<const char *, long double> read_sums;

	tr1::unordered_map<const char *, vector<EM_Map*> >::iterator read_iterator;
	tr1::unordered_map<const char *, vector<EM_Map*> >::iterator isoform_iterator;

	int counter = 0;
	for(read_iterator=read2emmap.begin(); read_iterator != read2emmap.end(); read_iterator++){

		double sumread = 0;

		const char * read_id = read_iterator->first;
		vector<EM_Map*> emmaps = read_iterator->second;
		for(unsigned int i=0; i < emmaps.size(); i++){
			sumread += emmaps.at(i)->em_prob(th[emmaps.at(i)->isoform.c_str()]);
		}
//		cout << "sumread=" << sumread << endl;
		read_sums[read_id] = sumread;
		counter++;
//		if(counter %100000 == 0){
//		cout << counter << " reads processed." << endl;
//		}
	}
	
	counter = 0;
	for(isoform_iterator=isoform2emmap.begin(); isoform_iterator != isoform2emmap.end(); isoform_iterator++){

		double sumterm = 0;
		const char * isoform_id = isoform_iterator->first;
		vector<EM_Map*> isoemmaps = isoform_iterator->second;
		for(unsigned int i=0; i < isoemmaps.size(); i++){
			sumterm += isoemmaps.at(i)->em_prob(th[isoform_id])/read_sums[isoemmaps.at(i)->base_read_id.c_str()];
		}
//		cout << "sumterm=" << sumterm << " N=" << N << endl;
		newth[isoform_id] = sumterm/N;
		counter++;
//		if(counter % 100000 == 0){
//		cout << counter << " isoforms processed." << endl;
//		}

	}



}



int main(int argc, char * argv[]){

	char * btfilename = argv[1];
//	cout << argv[2] << " " << atoi(argv[2]) << endl;
	int offset = atoi(argv[2]);
	ifstream btstream(btfilename);
	BowtieEntry bt1(offset);
	BowtieEntry bt2(offset); //Initialize with offset

	tr1::unordered_map<const char *, vector<EM_Map*> > read_to_emmaps;
	tr1::unordered_map<const char *, vector<EM_Map*> > isoform_to_emmaps;
	tr1::unordered_map<const char *, vector<EM_Map*> >::iterator isoform_iterator;
	tr1::unordered_map<const char *, long double> theta;
	tr1::unordered_map<const char *, long double> newtheta;
	int counter = 0;
	while(btstream >> bt1 && !btstream.eof()){
//		cout << btstream.eof() << endl;
//		btstream >> bt1;
		btstream >> bt2;

		EM_Map * pemmap;
		pemmap = new EM_Map(bt1, bt2, 10000); ///BADDDDDDDD, replace 10000
		read_to_emmaps[pemmap->base_read_id.c_str()].push_back(pemmap);
//		cout << "isoform " << pemmap->isoform.c_str() << endl;
		isoform_to_emmaps[pemmap->isoform.c_str()].push_back(pemmap);
		
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
	//	cout << 1/(long double)M << endl;
		++isoform_iterator;
	}
	
	cout << "Made first theta guess." << endl;
	int count = 0;

	while(count < 1000){
		cout << "Starting iteration " << count << endl;
		EM_Update(read_to_emmaps, isoform_to_emmaps, theta, newtheta, N);
		theta = newtheta;
		count++;
	}

	ofstream outstream("em_output");

	for(tr1::unordered_map<const char *, long double>::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
		outstream << thiter->first << "\t" << thiter->second << endl;
	}
	return 0;

}
