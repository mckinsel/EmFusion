#include <fstream>
#include <algorithm>
#include "main.h"
#include "BowtieEntry.h"
#include "FastaEntry.h"
#include "Read.h"

void clear_bowtie_vector(vector<BowtieEntry*> & btv) {
	for(unsigned int i = 0; i < btv.size(); i++) {
		delete btv[i];
	}
	btv.clear();
}

int Sift_main(int argc, char * argv[]) {

	char * bowtiefile1 = argv[1];
	char * bowtiefile2 = argv[2];
	char * fastqfile1 = argv[3];
	char * fastqfile2 = argv[4];
	int offset = atoi(argv[5]);

	bool bustfiles;
	char * bustfile1 = NULL;
	char * bustfile2 = NULL;
	if(argc == 8){
		bustfiles = true;
		bustfile1 = argv[6];
		bustfile2 = argv[7];
	} else {
		bustfiles = false;
	}
	cout << "Bowtiefile1 " << bowtiefile1 << endl;
	cout << "Bowtiefile2 " << bowtiefile2 << endl;
	cout << "fastqfile1 " << fastqfile1 << endl;
	cout << "fastqfile2 " << fastqfile2 << endl;
	cout << "bustfile1 " << bustfile1 << endl;
	cout << "bustfile2 " << bustfile2 << endl;
	cout << "offset " << offset << endl;
//	For now, we'll just assume that these are all sorted.
//	In the future, we should implement an external sort here.
//	Or not, it probably doesn't matter.

	ifstream bowtiestream1(bowtiefile1);
	ifstream bowtiestream2(bowtiefile2);

	ifstream fastqstream1(fastqfile1);
	ifstream fastqstream2(fastqfile2);

	ifstream buststream1;
	ifstream buststream2;

	if(bustfiles) {
		buststream1.open(bustfile1);
		buststream2.open(bustfile2);
	}

	char outfile[250];
	strcpy(outfile, fastqfile1);
	ofstream mapping_distance_stream(strcat(outfile, ".mapdist"));
	strcpy(outfile, fastqfile1);
	ofstream discordant_mapping_stream(strcat(outfile, ".discord"));
	strcpy(outfile, fastqfile1);
	ofstream unmapped_fasta_stream(strcat(outfile, ".unmapped.fasta"));

	Read read1(offset);
	Read read2(offset);

	FastaEntry bust1;
	FastaEntry bust2;

	string cur_bust1_id = "";
	string cur_bust2_id = "";
	string cur_bowtie1_id = "";
	string cur_bowtie2_id = "";

	vector<BowtieEntry*> bowtieentries1;
	vector<BowtieEntry*> bowtieentries2;

	BowtieEntry * bte1 = new BowtieEntry(offset);
	BowtieEntry * bte2 = new BowtieEntry(offset);

	int iter_counter = 0;

	while(fastqstream1 >> read1 && fastqstream2 >> read2) {
		++iter_counter;
		if(iter_counter % 10000 == 0) {
			cout << "Working on read " << iter_counter << endl;
		}
//		cout << "working on read " << read1.base_id << " " << read2.base_id << endl;

		if(bustfiles){
//			cout << "Now looking in repeat files" << endl;
			while(cur_bust1_id.compare(read1.id) < 0 && buststream1 >> bust1){
//				cout << "Getting next record from repeat file 1" << endl;
				cur_bust1_id = bust1.id;
//				cout << "reading from repeat file 1 " << bust1.id << endl;
			}

			while(cur_bust2_id.compare(read2.id) < 0 && buststream2 >> bust2){
				cur_bust2_id = bust2.id;
//				cout << "reading from repeat file 2 " << bust2.id << endl;
			}
		}

		if(cur_bust1_id == read1.id || cur_bust2_id == read2.id) {
//			cout << "Match in repeat found, so continuing" << endl;
			continue;
		}

		clear_bowtie_vector(bowtieentries1);
		clear_bowtie_vector(bowtieentries2);

//		cout << "Bowtie entry vectors cleared" << endl;

		if(bte1->read_id == read1.id){ //See if there's a good match from the last iteration
			bowtieentries1.push_back(bte1);
			bte1 = new BowtieEntry(offset);
		}

		while(cur_bowtie1_id.compare(read1.id) <= 0 && bowtiestream1 >> *bte1) {
			if(bte1->read_id == read1.id){
				bowtieentries1.push_back(bte1);
				bte1 = new BowtieEntry(offset);
			}
			cur_bowtie1_id = bte1->read_id;

		}

		if(bte2->read_id == read2.id){ //See if there's a good match from the last iteration
			bowtieentries2.push_back(bte2);
			bte2 = new BowtieEntry(offset);
		}

		while(cur_bowtie2_id.compare(read2.id) <= 0 && bowtiestream2 >> *bte2) {
			if(bte2->read_id == read2.id){
				bowtieentries2.push_back(bte2);
				bte2 = new BowtieEntry(offset);
			}
			cur_bowtie2_id = bte2->read_id;

		}
//		cout << "Bowtie entries have been read in." << endl;
//		cout << "Current BT entries are: " << cur_bowtie1_id << " " << cur_bowtie2_id << endl;


		if(bowtieentries1.size() > 0 && bowtieentries2.size() > 0){ //If both ends of the mate pair have mappings
			vector<string> genes1;
			vector<string> genes2;

//			See if there is an overlap between the mapped genes for each side.
			for(unsigned int i = 0; i < bowtieentries1.size(); i++) {
				genes1.push_back(bowtieentries1.at(i)->mapped_gene());
			}
			for(unsigned int i = 0; i < bowtieentries2.size(); i++) {
				genes2.push_back(bowtieentries2.at(i)->mapped_gene());
			}
			sort(genes1.begin(), genes1.end());
			sort(genes2.begin(), genes2.end());

			for(unsigned int i = 0; i < genes1.size(); i++) {
//				cout << "Genes1 : " << genes1.at(i) << endl;
			}

			for(unsigned int i = 0; i < genes2.size(); i++) {
//				cout << "Genes2 : " << genes2.at(i) << endl;
			}
			vector<string> g_intersection(genes1.size() + genes2.size());
			vector<string>::iterator end_it;
			end_it = set_intersection(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), g_intersection.begin());

//			cout << "int(end_it - g_intersection.begin()) " << int(end_it - g_intersection.begin()) << endl;
			if (int(end_it - g_intersection.begin()) > 0){ //There is at least one concordant mapping

//				Write out distances between mappings on the same transcript, if there are any.
				for(unsigned int i = 0; i < bowtieentries1.size(); i++){
					for(unsigned int j = 0; j < bowtieentries2.size(); j++){
						if(bowtieentries1.at(i)->mapped_transcript() == bowtieentries2.at(j)->mapped_transcript()) {
//							cout << "Writing out pair " << bowtieentries1.at(i)->mapped_transcript() << " " << bowtieentries2.at(j)->mapped_transcript() << endl;
							mapping_distance_stream << abs(bowtieentries1.at(i)->position - bowtieentries2.at(j)->position) << endl;
						}
					}
				}
			} else { //Only discordant mappings

//				Write out discordant mappings
				for(unsigned int i = 0; i < bowtieentries1.size(); i++){
					for(unsigned int j = 0; j < bowtieentries2.size(); j++){
						vector<string> transcripts (2);
						vector<string> genes (2);
						vector<int> positions (2);

						if(bowtieentries1.at(i)->mapped_transcript().compare(bowtieentries2.at(j)->mapped_transcript()) < 0) {
							discordant_mapping_stream << bowtieentries1.at(i)->mapped_transcript() << "\t";
							discordant_mapping_stream << bowtieentries2.at(j)->mapped_transcript() << "\t";
							discordant_mapping_stream << bowtieentries1.at(i)->mapped_gene() << "\t";
							discordant_mapping_stream << bowtieentries2.at(j)->mapped_gene() << "\t";
							discordant_mapping_stream << bowtieentries1.at(i)->position << "\t";
							discordant_mapping_stream << bowtieentries2.at(j)->position << endl;
						} else if (bowtieentries1.at(i)->mapped_transcript().compare(bowtieentries2.at(j)->mapped_transcript()) > 0) {
							discordant_mapping_stream << bowtieentries2.at(j)->mapped_transcript() << "\t";
							discordant_mapping_stream << bowtieentries1.at(i)->mapped_transcript() << "\t";
							discordant_mapping_stream << bowtieentries2.at(j)->mapped_gene() << "\t";
							discordant_mapping_stream << bowtieentries1.at(i)->mapped_gene() << "\t";
							discordant_mapping_stream << bowtieentries2.at(j)->position << "\t";
							discordant_mapping_stream << bowtieentries1.at(i)->position << endl;
						} else {
							cout << "THE PROGRAM IS BROKEN!" << endl;
							exit(1);
						}
					}
				}
			}
		}
		else if (bowtieentries1.size() == 0 && bowtieentries2.size() == 0){
			read1.write_as_fasta(unmapped_fasta_stream);
			read2.write_as_fasta(unmapped_fasta_stream);
		}
	}
	delete bte1;
	delete bte2;
	clear_bowtie_vector(bowtieentries1);
	clear_bowtie_vector(bowtieentries2);
	return 0;
}
