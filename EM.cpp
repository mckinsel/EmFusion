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
#include <float.h>

using namespace std;
//using namespace tr1;

typedef tr1::unordered_map<string, vector<EM_Map*> > vectorumap;
typedef tr1::unordered_map<string, long double > longdoubleumap;
typedef tr1::unordered_map<int, double> int2doubleumap;
typedef tr1::unordered_map<string, int> string2intumap;

bool covers_fusion_point(EM_Map & em){
	char * maptok;
	char * cstr;
	string transcript;
	vector<char*> transcript_parts;

	int FUZZ = 8;

	cstr = new char [em.isoform.size() + 1];
	strcpy(cstr, em.isoform.c_str());
	maptok = strtok(cstr, "|");
	transcript = string(maptok);
	delete[] cstr;

	cstr = new char [transcript.size() + 1];
	strcpy(cstr, transcript.c_str());
	maptok = strtok(cstr, "_");

	while(maptok != NULL) {
		transcript_parts.push_back(maptok);
		maptok = strtok(NULL, "_");
	}

	int fusion_point = atoi( transcript_parts.at(transcript_parts.size() - 2) );
//    if(em.base_read_id == "MORF4_PA2G4_683_889_117"){
//        cout << fusion_point << " " << em.start << " " << em.end << endl;
//    }
	if(fusion_point > em.start + FUZZ && fusion_point < em.end - FUZZ) {
		return true;
	} else {
		return false;
	}
}

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

//	Random_EM_Map * remmap;
//	remmap = new Random_EM_Map();

//	ifstream randomstream(randomfilename.c_str());

//	while (randomstream >> *remmap && !randomstream.eof()) {
//		readid_to_score[remmap->base_read_id] += remmap->em_prob(th["random"]);
//	}

//	randomstream.close();
//	delete remmap;

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
//	ifstream randomstream;
	EMstream.open(EMfilename.c_str());
//	cout << "Opened " << EMfilename << endl;
//	randomstream.open(randomfilename.c_str());

	EM_Map * emmap;
	emmap = new EM_Map();

	int counter = 0;

	while(EMstream >> *emmap && !EMstream.eof()) {
		counter++;
//		if(counter%5000000 == 0) cout << "On BT entry " << counter << endl;
		read_iterator = read_sums.find(emmap->base_read_id);
		if (read_iterator == read_sums.end()) {
			read_sums[emmap->base_read_id] = 0;
//			cout << "Setting " << emmap->base_read_id << " to zero." << endl;
		}
		read_sums[emmap->base_read_id] += emmap->em_prob(th[emmap->isoform]);
//		cout << "read_sums at " << emmap->base_read_id << " " << read_sums[emmap->base_read_id] << endl;
	}

	for(longdoubleumap::iterator it = read_sums.begin(); it != read_sums.end(); it++) {
        if(it->second == 0)
		cout << "read_sums at " << it->first << " is " << it->second << endl;
	}

//	Random_EM_Map * remmap;
//	remmap = new Random_EM_Map();

//	while (randomstream >> *remmap && !randomstream.eof()) {
//		read_sums[remmap->base_read_id] += remmap->em_prob(th[remmap->isoform]);
//	}

//	for(longdoubleumap::iterator it = read_sums.begin(); it != read_sums.end(); it++) {
//		cout << "after random, read_sums at " << it->first << " is " << it->second << endl;
//	}

//	cout << "Done reading BT file the first time." << endl;
	counter = 0;
	EMstream.close();
//	randomstream.close();

	EMstream.open(EMfilename.c_str());
//	randomstream.open(randomfilename.c_str());

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

//            if(sumterm != sumterm) sumterm = DBL_MIN;
			newth[current_isoform_id] = sumterm/N;
//            cout << "Setting th for " << current_isoform_id << " to " << sumterm << "/" << N << "=" << sumterm/N <<endl;
//            cout << "It used to be " << th[current_isoform_id] << endl;
			clear_EM_Map_vector(emmap_vect);
		}

		counter++;
		current_isoform_id = emmap->isoform;
//		if(counter % 5000000 == 0) cout << "On BT entry " << counter << endl;
		emmap_vect.push_back(emmap);
		emmap = new EM_Map();
	}

	long double sumterm = 0;
	for(unsigned int i=0; i < emmap_vect.size(); i++) {
		sumterm += emmap_vect.at(i)->em_prob(th[current_isoform_id])/read_sums[emmap_vect.at(i)->base_read_id];
	}

//    if(sumterm != sumterm) sumterm = DBL_MIN;
	newth[current_isoform_id] = sumterm/N;
//    cout << "Setting th for " << current_isoform_id << " to " << sumterm << "/" << N << "=" << sumterm/N <<endl;
//    cout << "It used to be " << th[current_isoform_id] << endl;
	clear_EM_Map_vector(emmap_vect);

//	while (randomstream >> *remmap && !randomstream.eof()) {
//		current_isoform_id = remmap->isoform;
//		emmap_vect.push_back(remmap);
//		remmap = new Random_EM_Map();
//	}

//	sumterm = 0;
//	for(unsigned int i=0; i < emmap_vect.size(); i++) {
//		sumterm += emmap_vect.at(i)->em_prob(th[current_isoform_id])/read_sums[emmap_vect.at(i)->base_read_id];
//	}

//	newth[current_isoform_id] = sumterm/N;
//	clear_EM_Map_vector(emmap_vect);

	delete emmap;
//	delete remmap;
}

void get_initialization_values( char * i_file, longdoubleumap & th){

    ifstream init_stream(i_file);
    char nextline[1500];
    long double thval;
    string trans_id;
    char * tok;
    char * cstr;

    while(!init_stream.eof()){

        init_stream.getline(nextline, (streamsize)1500);
        if(init_stream.eof()) continue;
        cstr = new char [1500];
        strcpy(cstr, nextline);

        tok = strtok(cstr, "\t");
        trans_id = string(tok);

        tok = strtok(NULL, "\t");
        thval = atof(tok);

        th[trans_id] = thval;

        delete[] cstr;
        }

    }

int EM_main(int argc, char * argv[]){

	char * btfilename_mappingsorted = argv[1];
	char * dist_prob_file = argv[2];
	char * reference_fasta = argv[3];
//	char * unmapped_fasta = argv[4];
	int offset = atoi(argv[5]);

//    cout << "Set required parameters. " << argc << endl;
    char * initialization_file = NULL;
    if(argc > 6){
//        cout << "Setting initialization file." << endl;
        initialization_file = argv[6];
    }

    string already = "";
    bool already_run;
    if(argc > 7) already=string(argv[7]);

    if(already == "true"){
        already_run = true;
    } else {
        already_run = false;
    }


	char EMfilename[150];
	strcpy(EMfilename, btfilename_mappingsorted);
	strcat(EMfilename, ".emint");

    char randomfilename[150];
//	strcpy(randomfilename, btfilename_mappingsorted);
//	strcat(randomfilename, ".random");

    int N;
    int M;

	longdoubleumap theta;
	longdoubleumap newtheta;

	set<string> seen_readids;
	set<string>::iterator read_iterator;

	set<string> seen_isoforms;
	set<string>::iterator isoform_iterator;

    if(!already_run){
//	cout << "Processed arguments." << endl;
//	Get Markov Chain for unmapped

//	MarkovChain mc(12, 1);
//	double log_unmapped_prob;
//	log_unmapped_prob = getMarkovChain(unmapped_fasta, mc);

//	cout << "Built MarkovChain." << endl;

//	Get mapping distance distribution
	int2doubleumap dist_prob;

	ifstream dpstream(dist_prob_file);
	string nextint, nextprob;

	while(dpstream >> nextint){
		dpstream >> nextprob;
		assert(atof(nextprob.c_str()) <= 1);
//		cout << nextint << "\t" << nextprob << endl;
		dist_prob[atoi(nextint.c_str())] = atof(nextprob.c_str());
	}

//	cout << "Got distance probabilities." << endl;

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

//	ofstream randomstream(randomfilename);

	BowtieEntry bt1(offset);
	BowtieEntry bt2(offset); //Initialize with offset

//	vectorumap read_to_emmaps;
//	vectorumap isoform_to_emmaps;

	int counter = 0;


	while(btstream >> bt1 && btstream >> bt2){

		EM_Map * pemmap;
		pemmap = new EM_Map(bt1, bt2, dist_prob, isoform_lengths);
		read_iterator = seen_readids.find(pemmap->base_read_id);
		if(read_iterator == seen_readids.end()) {
//          Random_EM_Map * prem;
//			prem = new Random_EM_Map(bt1, bt2, mc);
//			randomstream << *prem;
//			delete prem;
			seen_readids.insert(pemmap->base_read_id);
		}
		seen_isoforms.insert(pemmap->isoform);
		emintstream << *pemmap;
		delete pemmap;
		counter++;
//		if(counter%1000000==0) cout << "Reading Bowtie pair " << counter << endl;
	}
	emintstream.close();
//	randomstream.close();

//	cout << "Bowtie reads read in." << endl;

	N = (int)seen_readids.size();
	M = (int)seen_isoforms.size()/* + 1*/; //Don't forget random
//	cout << "N " << N << " M " << M << endl;

    } else { //End of if !already_run block
        ifstream emintstream(EMfilename);

        char nextline[5000];
        char * cstr;
        char * tok;

        int rcounter = 0;
        while(!emintstream.eof()){
            emintstream.getline(nextline, (streamsize)5000);
            if(emintstream.eof()) continue;
            cstr = new char[5000];
            strcpy(cstr, nextline);

            tok = strtok(cstr, "\t");

            seen_readids.insert(string(tok));
//            cout << "adding " << string(tok) << " to seen_readids" << endl;

            tok = strtok(NULL, "\t");

            seen_isoforms.insert(string(tok));
//            cout << "adding " << string(tok) << " to seen_isoforms" << endl;
            delete cstr;
            rcounter ++;
            if(rcounter % 1000000 == 0) cout << rcounter << endl;
        }
        N = (int)seen_readids.size();
        M = (int)seen_isoforms.size()/* + 1*/;
//    	cout << "N " << N << " M " << M << endl;
    }
    if(argc == 6){
	    isoform_iterator = seen_isoforms.begin();
	    while(isoform_iterator != seen_isoforms.end()){
    		theta[*isoform_iterator] = 1/(long double)M;
		    ++isoform_iterator;
	    }
//	    theta["random"] = 1/(long double)M;

//	    cout << "Made first theta guess." << endl;
    } else if (argc > 6){
       get_initialization_values(initialization_file, theta);
    }


	int count = 0;
	long double log_diff = 1;
	long double oldll = 0;
	long double newll = 0;

	ofstream intoutstream;

	while(log_diff > 1e-5){
		cout << "Starting iteration " << count << endl;
		if(count % 1 == 0) oldll = log_likelihood(EMfilename, randomfilename, theta);
		EM_Update(EMfilename, randomfilename, theta, newtheta, N);
		if(count % 1 == 0) newll = log_likelihood(EMfilename, randomfilename, newtheta);
		if(count % 1 == 0){
			cout << " Old log likelihood is " << oldll << endl;
			cout << " New log likelihood is " << newll << endl;
			cout << " Difference is " << abs(oldll - newll)/abs(oldll) << endl;
		}


		if(count % 1 == 0) {
			intoutstream.open("em_output.int" );
			for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
         		       intoutstream << thiter->first << "\t" << setprecision(35) << thiter->second << endl;
		        }
			intoutstream.close();
		}
		if(count % 1 == 0) log_diff = abs(oldll - newll)/abs(oldll);
		theta = newtheta;
		count++;
	}

	ofstream outstream("em_output");

	for(longdoubleumap::iterator thiter=newtheta.begin(); thiter != newtheta.end(); thiter++){
		outstream << thiter->first << "\t" << thiter->second << endl;
	}

//	cout << "Final mapped log likelihood : " << setprecision(15) << newll << endl;
//	cout << "Unmapped log likelihood : " << setprecision(15) << log_unmapped_prob << endl;
//	cout << "Number of reads : " << N << endl;
//	cout << "Number of isoforms : " << M << endl;


	ifstream EMstream;
	EMstream.open(EMfilename);
	EM_Map * emmap;
	emmap = new EM_Map();

	bool is_fusion;

	tr1::unordered_map<string, bool> is_read_fusion;
	tr1::unordered_map<string, bool>::iterator is_read_fusion_iterator;

    set<string> concordant_reads;
    set<string>::iterator criterator;
    ifstream concordantstream("concordant.reads");
    string nextread;

    while(concordantstream >> nextread){
        concordant_reads.insert(nextread);
    }


	while(EMstream >> *emmap && !EMstream.eof()) {

		if(emmap->isoform.substr(0,2) == "F_"){
			is_fusion = true;
		} else {
			is_fusion = false;
		}

		is_read_fusion_iterator = is_read_fusion.find(emmap->base_read_id);
		if(is_read_fusion_iterator == is_read_fusion.end()){
			is_read_fusion[emmap->base_read_id] = is_fusion;
		} else {
			if(is_read_fusion[emmap->base_read_id] == true){ //If it's true, set it to what we find. If it's false, it stays false
				is_read_fusion[emmap->base_read_id] = is_fusion;
			}
		}
	}

    for(is_read_fusion_iterator = is_read_fusion.begin(); is_read_fusion_iterator != is_read_fusion.end(); is_read_fusion_iterator++){
        criterator = concordant_reads.find(is_read_fusion_iterator->first);
        if(criterator != concordant_reads.end()){
//            cout << "Setting read " << is_read_fusion_iterator->first << " to false because of concordant reads." << endl;
            is_read_fusion[is_read_fusion_iterator->first] = false;
        }
    }
//	cout << "Done setting is_fusion." << endl;

    ofstream fusreadsstream("fusion_reads.txt");
    for(is_read_fusion_iterator = is_read_fusion.begin(); is_read_fusion_iterator != is_read_fusion.end(); is_read_fusion_iterator++){
        if(is_read_fusion_iterator->second == true){
            fusreadsstream << is_read_fusion_iterator->first << endl;
        }
    }
    fusreadsstream.close();
//    cout << "Is MORF4_PA2G4_683_889_117? " << boolalpha << is_read_fusion["MORF4_PA2G4_683_889_117"] << endl;

	EMstream.close();
	EMstream.open(EMfilename);

	longdoubleumap f_read_sums;
	longdoubleumap::iterator f_read_iterator;

	while(EMstream >> *emmap && !EMstream.eof()) {
		if(!is_read_fusion[emmap->base_read_id]) continue;
		f_read_iterator = f_read_sums.find(emmap->base_read_id);
		if (f_read_iterator == f_read_sums.end()) {
			f_read_sums[emmap->base_read_id] = 0;
		}
		f_read_sums[emmap->base_read_id] += emmap->em_prob(theta[emmap->isoform]);
	}
	EMstream.close();


//	ifstream f_randomstream;
//	f_randomstream.open(randomfilename);
//	Random_EM_Map * remmap;
//	remmap = new Random_EM_Map();

//	while (f_randomstream >> *remmap && !f_randomstream.eof()) {
//		f_read_sums[remmap->base_read_id] += remmap->em_prob(theta[remmap->isoform]);
//	}

//	f_randomstream.close();
//	delete remmap;

	cout << "Done adding read_sums." << endl;

	EMstream.open(EMfilename);
	long double read_isoform_prob;
	longdoubleumap final_fusion_counts;
	longdoubleumap final_fusion_probs;
	longdoubleumap::iterator ffc_iterator;
	string f_current_isoform_id = "";

	while (EMstream >> *emmap && !EMstream.eof()) {
		if(!is_read_fusion[emmap->base_read_id]) continue;
		if(!covers_fusion_point(*emmap)) continue;
		read_isoform_prob = emmap->em_prob(theta[emmap->isoform]) / f_read_sums[emmap->base_read_id];
		ffc_iterator = final_fusion_counts.find(emmap->isoform);
		if(ffc_iterator == final_fusion_counts.end()){
			final_fusion_counts[emmap->isoform] = 0;
			final_fusion_probs[emmap->isoform] = 1;
		}
		final_fusion_counts[emmap->isoform] += read_isoform_prob;
		final_fusion_probs[emmap->isoform] *= (1 - read_isoform_prob);
/*        if(read_isoform_prob > .01){
            cout << emmap->base_read_id << " " << emmap->isoform << " " << read_isoform_prob << endl;
        } else if(emmap->base_read_id == "MORF4_PA2G4_683_889_117") {
            cout << emmap->base_read_id << " " << emmap->isoform << " " << read_isoform_prob << endl;
        }*/

	}

	delete emmap;

	ofstream finaloutstream;
	ofstream finalprobstream;

	finaloutstream.open("fusion_counts");
	for(ffc_iterator=final_fusion_counts.begin(); ffc_iterator != final_fusion_counts.end(); ffc_iterator++){
		finaloutstream << ffc_iterator->first << "\t" << ffc_iterator->second << endl;
	}
	finaloutstream.close();

	finalprobstream.open("fusion_probs");
	for(ffc_iterator=final_fusion_probs.begin(); ffc_iterator != final_fusion_probs.end(); ffc_iterator++){
		finalprobstream << ffc_iterator->first << "\t" << 1 - ffc_iterator->second << endl;
	}
	finalprobstream.close();
	return 0;

}



