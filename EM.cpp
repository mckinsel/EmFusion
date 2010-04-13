/*
 * EM.cpp
 *
 *  Created on: Mar 3, 2010
 *      Author: marcus
 */

#include <fstream>
#include <tr1/unordered_map>
#include <string>
#include <set>
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


void clear_EM_Map_vector(vector<EM_Map*> & emm) {
	for(unsigned int i = 0; i < emm.size(); i++) {
		delete emm.at(i);
	}
	emm.clear();
}
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

long double log_likelihood(string emfilename, string randomfilename, longdoubleumap & th){

	long double log_score = 0;
	longdoubleumap::iterator read_iterator;
	longdoubleumap readid_to_score;

	EM_Map * emmap;
	emmap = new EM_Map();

	ifstream emstream(emfilename.c_str());

	while (emstream >> *emmap && !emstream.eof()) {
		read_iterator = readid_to_score.find(emmap->base_read_id);
		if(read_iterator == readid_to_score.end()) {
			readid_to_score[emmap->base_read_id] = 0;
		}
		readid_to_score[emmap->base_read_id] += emmap->em_prob(th[emmap->isoform]);
	}

	emstream.close();
	delete emmap;

	Random_EM_Map * remmap;
	remmap = new Random_EM_Map();

	ifstream randomstream(randomfilename.c_str());

	while (randomstream >> *remmap && !randomstream.eof()) {
		readid_to_score[remmap->base_read_id] += remmap->em_prob(th["random"]);
	}

	randomstream.close();
	delete remmap;

	for(read_iterator = readid_to_score.begin(); read_iterator != readid_to_score.end(); read_iterator++) {
		log_score += log(read_iterator->second);
	}
	return log_score;
}
void EM_Update( string EMfilename, string randomfilename, longdoubleumap & th, longdoubleumap & newth, int N){

	longdoubleumap read_sums;

	longdoubleumap::iterator read_iterator;
	vectorumap::iterator isoform_iterator;

	ifstream EMstream;
	ifstream randomstream;

	EMstream.open(EMfilename.c_str());
	randomstream.open(randomfilename.c_str());

	EM_Map * emmap;
	emmap = new EM_Map();

	int counter = 0;

	while(EMstream >> *emmap && !EMstream.eof()) {
		counter++;
		if(counter%5000000 == 0) cout << "On BT entry " << counter << endl;
		read_iterator = read_sums.find(emmap->base_read_id);
		if (read_iterator == read_sums.end()) {
			read_sums[emmap->base_read_id] = 0;
//			cout << "Setting " << emmap->base_read_id << " to zero." << endl;
		}
		read_sums[emmap->base_read_id] += emmap->em_prob(th[emmap->isoform]);
//		cout << "read_sums at " << emmap->base_read_id << " " << read_sums[emmap->base_read_id] << endl;
	}

//	for(longdoubleumap::iterator it = read_sums.begin(); it != read_sums.end(); it++) {
//		cout << "read_sums at " << it->first << " is " << it->second << endl;
//	}

	Random_EM_Map * remmap;
	remmap = new Random_EM_Map();

	while (randomstream >> *remmap && !randomstream.eof()) {
		read_sums[remmap->base_read_id] += remmap->em_prob(th[remmap->isoform]);
	}

//	cout << "Done reading BT file the first time." << endl;
	counter = 0;
	EMstream.close();
	randomstream.close();

	EMstream.open(EMfilename.c_str());
	randomstream.open(randomfilename.c_str());

	string current_isoform_id = "";
	vector<EM_Map*> emmap_vect;


	while (EMstream >> *emmap && !EMstream.eof()) {
//		cout << "Working on read " << emmap->base_read_id << " mapped to " << emmap->isoform << endl;
//		cout << "The 'current_isoform_id' is " << current_isoform_id << endl;
		if(emmap->isoform != current_isoform_id && current_isoform_id != "") {
//			cout << "A new isoform!" << endl;
			long double sumterm = 0;
			for(unsigned int i=0; i < emmap_vect.size(); i++) {
				sumterm += emmap_vect.at(i)->em_prob(th[current_isoform_id])/read_sums[emmap_vect.at(i)->base_read_id];
			}
			newth[current_isoform_id] = sumterm/N;
			clear_EM_Map_vector(emmap_vect);
		}

		counter++;
		current_isoform_id = emmap->isoform;
		if(counter % 5000000 == 0) cout << "On BT entry " << counter << endl;
		emmap_vect.push_back(emmap);
		emmap = new EM_Map();
	}

	long double sumterm = 0;
	for(unsigned int i=0; i < emmap_vect.size(); i++) {
		sumterm += emmap_vect.at(i)->em_prob(th[current_isoform_id])/read_sums[emmap_vect.at(i)->base_read_id];
	}
	newth[current_isoform_id] = sumterm/N;
	clear_EM_Map_vector(emmap_vect);

	while (randomstream >> *remmap && !randomstream.eof()) {
		current_isoform_id = remmap->isoform;
		emmap_vect.push_back(remmap);
		remmap = new Random_EM_Map();
	}

	sumterm = 0;
	for(unsigned int i=0; i < emmap_vect.size(); i++) {
		sumterm += emmap_vect.at(i)->em_prob(th[current_isoform_id])/read_sums[emmap_vect.at(i)->base_read_id];
	}
	newth[current_isoform_id] = sumterm/N;
	clear_EM_Map_vector(emmap_vect);

	delete emmap;
	delete remmap;
}



int EM_main(int argc, char * argv[]){

	char * btfilename_mappingsorted = argv[1];
	char * dist_prob_file = argv[2];
	char * reference_fasta = argv[3];
	char * unmapped_fasta = argv[4];
	int offset = atoi(argv[5]);

	char EMfilename[150];
	strcpy(EMfilename, btfilename_mappingsorted);
	strcat(EMfilename, ".emint");

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

	cout << "Got distance probabilities." << endl;

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
	ifstream btstream(btfilename_mappingsorted);
	ofstream emintstream(EMfilename);

	char randomfilename[150];
	strcpy(randomfilename, btfilename_mappingsorted);
	strcat(randomfilename, ".random");
	ofstream randomstream(randomfilename);

	BowtieEntry bt1(offset);
	BowtieEntry bt2(offset); //Initialize with offset

//	vectorumap read_to_emmaps;
//	vectorumap isoform_to_emmaps;

	set<string> seen_readids;
	set<string>::iterator read_iterator;

	set<string> seen_isoforms;
	set<string>::iterator isoform_iterator;


	longdoubleumap theta;
	longdoubleumap newtheta;
	int counter = 0;



	while(btstream >> bt1 && btstream >> bt2){

		EM_Map * pemmap;
		pemmap = new EM_Map(bt1, bt2, dist_prob, isoform_lengths);

		read_iterator = seen_readids.find(pemmap->base_read_id);
		if(read_iterator == seen_readids.end()) {
			Random_EM_Map * prem;
			prem = new Random_EM_Map(bt1, bt2, mc);
			randomstream << *prem;
			delete prem;
			seen_readids.insert(pemmap->base_read_id);
		}
		seen_isoforms.insert(pemmap->isoform);
		emintstream << *pemmap;
		delete pemmap;
		counter++;
		if(counter%1000000==0) cout << "Reading Bowtie pair " << counter << endl;
	}
	emintstream.close();
	randomstream.close();

	cout << "Bowtie reads read in." << endl;

	int N = (int)seen_readids.size();
	int M = (int)seen_isoforms.size() + 1; //Don't forget random
	cout << "N " << N << " M " << M << endl;

	isoform_iterator = seen_isoforms.begin();
	while(isoform_iterator != seen_isoforms.end()){
		theta[*isoform_iterator] = 1/(long double)M;
		++isoform_iterator;
	}
	theta["random"] = 1/(long double)M;

	cout << "Made first theta guess." << endl;
	int count = 0;
	long double log_diff = 1;
	long double oldll = 0;
	long double newll = 0;

	ofstream intoutstream;

	while(log_diff > 1e-8){
		cout << "Starting iteration " << count << endl;
		if(count % 10 == 0) oldll = log_likelihood(EMfilename, randomfilename, theta);
		EM_Update(EMfilename, randomfilename, theta, newtheta, N);
		if(count % 10 == 0) newll = log_likelihood(EMfilename, randomfilename, newtheta);
		if(count % 10 == 0){
			cout << " Old log likelihood is " << oldll << endl;
			cout << " New log likelihood is " << newll << endl;
			cout << " Difference is " << abs(oldll - newll)/abs(oldll) << endl;
		}
		

		if(count % 10 == 0) {
			intoutstream.open("em_output.int" );
			for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
         		       intoutstream << thiter->first << "\t" << thiter->second << endl;
		        }
			intoutstream.close();
		}
		if(count % 10 == 0) log_diff = abs(oldll - newll)/abs(oldll);
		theta = newtheta;
		count++;
	}

	ofstream outstream("em_output");

	for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
		outstream << thiter->first << "\t" << thiter->second << endl;
	}

//	cout << "Final mapped log likelihood : " << setprecision(15) << newll << endl;
	cout << "Unmapped log likelihood : " << setprecision(15) << log_unmapped_prob << endl;
	cout << "Number of reads : " << N << endl;
	cout << "Number of isoforms : " << M << endl;
	return 0;

}
