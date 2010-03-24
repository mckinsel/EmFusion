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
#include <iomanip>
#include <cmath>
#include "EM_Map.h"
#include "Random_EM_Map.h"
#include "BowtieEntry.h"
#include "Utils.h"
#include "MarkovChain.h"
#include "main.h"

using namespace std;
using namespace tr1;



typedef tr1::unordered_map<string, vector<EM_Map*> > vectorumap;
typedef tr1::unordered_map<string, long double > longdoubleumap;
typedef tr1::unordered_map<int, double> int2doubleumap;
typedef tr1::unordered_map<string, int> string2intumap;

double getMarkovChain(string filename, MarkovChain& mc) {

//	First iterate through the unmapped file to build the MM
//	Then, iterate through it again to calculate the probabilities
//	of each of the unmapped sequences.

	double log_unmapped_prob = 0;
	string nextline;
	string fullseq = "";
	ifstream umfa(filename.c_str());
	while(umfa >> nextline){
		if(nextline.at(0) == '>'){
			if(fullseq.length() != 0){
				mc.add_sequence(fullseq);
				fullseq.clear();
			}

		} else {
			fullseq.append(nextline);
		}
	}
	mc.add_sequence(fullseq);

	umfa.close();

	ifstream umfa2(filename.c_str());

	while(umfa2 >> nextline){
			if(nextline.at(0) == '>' ){
				if(fullseq.length() != 0){
					log_unmapped_prob += log(mc.sequence_probability(fullseq));
					fullseq.clear();
				}
			} else {
				fullseq.append(nextline);
			}
		}
	log_unmapped_prob += log(mc.sequence_probability(fullseq));
	umfa2.close();
	return log_unmapped_prob;
}

long double log_likelihood(vectorumap & read2emmap, longdoubleumap & th){

	long double log_score = 0;
	vectorumap::iterator read_iterator;

	for(read_iterator=read2emmap.begin(); read_iterator != read2emmap.end(); read_iterator++){
		vector<EM_Map*> read_emmaps = read_iterator->second;

		long double read_score = 0;

		for(unsigned int i=0; i < read_emmaps.size(); i++){
			read_score += read_emmaps.at(i)->em_prob(th[read_emmaps.at(i)->isoform]);
		}
		if(read_score == 0){
			cout << read_iterator->first << endl;
		}
		log_score += log(read_score);
	}
	return log_score;
}
void EM_Update( vectorumap & read2emmap, vectorumap & isoform2emmap, longdoubleumap & th, longdoubleumap & newth, int N){

	longdoubleumap read_sums;

	vectorumap::iterator read_iterator;
	vectorumap::iterator isoform_iterator;

	//Calculate the denominators for the update rule.
	int counter = 0;
	for(read_iterator=read2emmap.begin(); read_iterator != read2emmap.end(); read_iterator++){

//		Iterate through every read and find the sum of all of its mappings.

		long double sumread = 0;

		string read_id = read_iterator->first;
		vector<EM_Map*> emmaps = read_iterator->second;
		for(unsigned int i=0; i < emmaps.size(); i++){
			sumread += emmaps.at(i)->em_prob(th[emmaps.at(i)->isoform]);
		}
		read_sums[read_id] = sumread;
//		if(sumread == 0 || true){
//			cout << "Read " << read_id << " sumread is " << sumread << endl;
//		}
		counter++;
	}

//	Now actually calculate a new theta.
	counter = 0;
	for(isoform_iterator=isoform2emmap.begin(); isoform_iterator != isoform2emmap.end(); isoform_iterator++){

//		Iterate through every isoform, including the random isoform.
		long double sumterm = 0;
		string isoform_id = isoform_iterator->first;
		vector<EM_Map*> isoemmaps = isoform_iterator->second;
		for(unsigned int i=0; i < isoemmaps.size(); i++){
			sumterm += isoemmaps.at(i)->em_prob(th[isoform_id])/read_sums[isoemmaps.at(i)->base_read_id];
		}
//		cout << "Isoform " << isoform_id << " sumterm is " << sumterm << endl;
		newth[isoform_id] = sumterm/N;
		counter++;
//		if(counter % 50000 == 0) cout << counter << " isoforms processed." << endl;

	}

}



int EM_main(int argc, char * argv[]){

	char * btfilename = argv[1];
	char * dist_prob_file = argv[2];
	char * reference_fasta = argv[3];
	char * unmapped_fasta = argv[4];
	int offset = atoi(argv[5]);

	cout << "Processed arguments." << endl;
//	Get Markov Chain for unmapped

	MarkovChain mc(3, 5);
	double log_unmapped_prob;
	log_unmapped_prob = getMarkovChain(unmapped_fasta, mc);

	cout << "Built MarkovChain." << endl;

//	Get mapping distance distribution
	int2doubleumap dist_prob;

	ifstream dpstream(dist_prob_file);
	string nextint, nextprob;

	while(dpstream >> nextint){
		dpstream >> nextprob;
		dist_prob[atoi(nextint.c_str())] = atof(nextprob.c_str());
	}

	cout << "Got distance probabilites." << endl;

//	Get gene lengths
	string2intumap isoform_lengths;
	string nextline;
	int lengthcounter = 0;
	string current_isoform = "nothing";
	ifstream refstream(reference_fasta);
	while(refstream >> nextline){
		if(nextline.at(0) == '>'){
			isoform_lengths[current_isoform] = lengthcounter;
			current_isoform = nextline.substr(1, nextline.length()-1);
			lengthcounter = 0;
		} else {
			lengthcounter += nextline.length();
		}
	}
	isoform_lengths[current_isoform] = lengthcounter;

//	Now iterate through bowtie file
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



	while(btstream >> bt1 && btstream >> bt2){

		EM_Map * pemmap;
		pemmap = new EM_Map(bt1, bt2, dist_prob, isoform_lengths);


//		cout << "Looking at mapping of " << bt1.base_read_id << " to " << bt1.mapping << endl;
//		cout << "It has a mapping probability of " << pemmap->P << endl;
		read_iterator = read_to_emmaps.find(pemmap->base_read_id);
		if(read_iterator == read_to_emmaps.end()){ //If it's the first time seeing the read.
			Random_EM_Map * prem;
			prem = new Random_EM_Map(bt1, bt2, mc);
			read_to_emmaps[prem->base_read_id].push_back(prem);
			isoform_to_emmaps[prem->isoform].push_back(prem);
		}
		read_to_emmaps[pemmap->base_read_id].push_back(pemmap);
		isoform_to_emmaps[pemmap->isoform].push_back(pemmap);

		delete bt1.read;
		delete bt2.read;
		counter++;
		if(counter%1000000==0) cout << "Reading Bowtie pair " << counter << endl;
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
	long double log_diff = 1;
	long double oldll, newll;

	while(log_diff > 1e-10){
		oldll = log_likelihood(read_to_emmaps, theta);
		EM_Update(read_to_emmaps, isoform_to_emmaps, theta, newtheta, N);
		newll = log_likelihood(read_to_emmaps, newtheta);
		if(count % 10 == 0){
			cout << "Starting iteration " << count << endl;
			cout << " Old log likelihood is " << oldll << endl;
			cout << " New log likelihood is " << newll << endl;
		}
		log_diff = abs(oldll - newll)/abs(oldll);
		theta = newtheta;
		count++;
	}

	ofstream outstream("em_output");

	for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
		outstream << thiter->first << "\t" << thiter->second << endl;
	}

	cout << "Final mapped log likelihood : " << setprecision(15) << newll << endl;
	cout << "Unmapped log likelihood : " << setprecision(15) << log_unmapped_prob << endl;
	cout << "Number of reads : " << N << endl;
	cout << "Number of isoforms : " << M << endl;
	return 0;

}
