#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include "../main.h"

using namespace std;


namespace {

TEST(EMTest, Resolves3MultimapCase) {

//	In this test, we have six reads. Two map to gene 1, two map to gene 2, and two map to a fusion of gene 1 and 2.
//	There are three reference sequences to disinguish between:
//         R1,R2
//						R3,R4
//				  R5,R6
//	G1 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//	G2 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//	FG GGGGGGGGGGGGGAAAAAAAAAAAAAAAAA

	char btfilename[50], dist_prob_file[50], reference_fasta[50], unmapped_fasta[50];
	strcpy(btfilename, "../Test/EM_MapTest.txt.sorted");
	strcpy(dist_prob_file, "../Test/EMTest_distprob.txt");
	strcpy(reference_fasta, "../Test/EMTest_reference.txt");
	strcpy(unmapped_fasta, "../Test/EMTest_unmapped.txt");
	char offs[] = "33";
	char self[] = "self";
	char * em_argv[] = {self, btfilename, dist_prob_file, reference_fasta, unmapped_fasta, offs};

	remove("em_output");
	EM_main(6, em_argv);


	ifstream emoutstream("em_output");
	EXPECT_TRUE(emoutstream);


	string nextline;
	while (emoutstream >> nextline) {
		string emo_transid = nextline;
		emoutstream >> nextline;
		double emo_theta = atof(nextline.c_str());

		if(emo_transid == "Trans1|Gene1|known|6|100|1000|1" || emo_transid == "Trans2|Gene2|known|6|100|1000|1") {
			EXPECT_LT(emo_theta, 1e-7);
		} else if(emo_transid == "TFus12|GFus12|known|2|100|1000|1") {
			EXPECT_GT(emo_theta, .9999);
		} else if(emo_transid == "random") {
			EXPECT_LT(emo_theta, 1e-10);
		} else {
			EXPECT_TRUE(false);
		}
	}

}

TEST(EMTest, Resolves4MultimapCase) {


//	In this test, we have eight reads. Two map to gene 1, four map to gene 2, and two map to a fusion of gene 1 and 2.
//	There are three reference sequences to disinguish between:
//         R1,R2
//	       GGGGG
//						R3,R4
//	                    AAAAA
//				  R5,R6
//				  GGAAA
//	       R7,R8
//	       AAAAA
//	G1 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
//	G2 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//	FG GGGGGGGGGGGGGAAAAAAAAAAAAAAAAA
	char btfilename[50], dist_prob_file[50], reference_fasta[50], unmapped_fasta[50];
	strcpy(btfilename, "../Test/EM_Test2.txt.sorted");
	strcpy(dist_prob_file, "../Test/EMTest_distprob.txt");
	strcpy(reference_fasta, "../Test/EMTest_reference.txt");
	strcpy(unmapped_fasta, "../Test/EMTest_unmapped.txt");
	char offs[] = "33";
	char self[] = "self";
	char * em_argv[] = {self, btfilename, dist_prob_file, reference_fasta, unmapped_fasta, offs};

	remove("em_output");
	EM_main(6, em_argv);


	ifstream emoutstream("em_output");
	EXPECT_TRUE(emoutstream);


	string nextline;
	while (emoutstream >> nextline) {
		string emo_transid = nextline;
		emoutstream >> nextline;
		double emo_theta = atof(nextline.c_str());

		if(emo_transid == "Trans1|Gene1|known|6|100|1000|1") {
			EXPECT_LT(emo_theta, 1e-7);
		} else if(emo_transid == "Trans2|Gene2|known|6|100|1000|1") {
			EXPECT_LT(abs(emo_theta - .333333333333), 1e-5);
		} else if(emo_transid == "TFus12|GFus12|known|2|100|1000|1") {
			EXPECT_LT(abs(emo_theta - .666666666666), 1e-5);
		} else if(emo_transid == "random") {
			EXPECT_LT(emo_theta, 1e-10);
		} else {
			EXPECT_TRUE(false);
		}
	}

}
}
