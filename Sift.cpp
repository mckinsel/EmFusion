#include <fstream>
#include "main.h"
#include "BowtieEntry.h"
#include "FastaEntry.h"
#include "Read.h"

int Sift_main(int argc, char * argv[]) {

	char * bowtiefile1 = argv[1];
	char * bowtiefile2 = argv[2];
	char * fastqfile1 = argv[3];
	char * fastqfile2 = argv[4];
	int offset = atoi(argv[5]);

	bool bustfiles;
	char * bustfile1, * bustfile2;
	if(argc == 8){
		bustfiles = true;
		bustfile1 = argv[6];
		bustfile2 = argv[7];
	} else {
		bustfiles = false;
	}

	ifstream fastqstream1(fastqfile1);
	ifstream fastqstream2(fastqfile2);

	ifstream buststream1;
	ifstream buststream2;

	if(bustfiles) {
		buststream1.open(bustfile1);
		buststream2.open(bustfile2);
	}

	Read read1(offset);
	Read read2(offset);

	FastaEntry bust1;
	FastaEntry bust2;

	string cur_bust1_id = "";
	string cur_bust2_id = "";

	while(fastqstream1 >> read1) {
		fastqstream2 >> read2;

		if(bustfiles){
			while(cur_bust1_id.compare(read1.id) < 0 && !buststream1.eof()){
				buststream1 >> bust1;
				cur_bust1_id = bust1.id;
			}

			while(cur_bust2_id.compare(read2.id) < 0 && !buststream2.eof()){
				buststream2 >> bust2;
				cur_bust2_id = bust2.id;
			}
		}


	}


	return 0;

}
