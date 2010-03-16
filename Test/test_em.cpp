#include <gtest/gtest.h>
#include <string>
#include "../main.h"

using namespace std;
//int EM_main(int argc, char * argv[]){
//
//	char * btfilename = argv[1];
//	char * dist_prob_file = argv[2];
//	char * reference_fasta = argv[3];
//	char * unmapped_fasta = argv[4];
//	int offset = atoi(argv[5]);

namespace {

class EMTest : public ::testing::Test {
protected:

	char btfilename[50], dist_prob_file[50], reference_fasta[50], unmapped_fasta[50];

	EMTest () {
		strcpy(btfilename, "../Test/EM_MapTest.txt");
		strcpy(dist_prob_file, "../Test/EM_Test_distprob.txt");
		strcpy(reference_fasta, "../Test/EMTest_reference.txt");
		strcpy(unmapped_fasta, "../Test/EMTest_unmapped.txt");
	}

virtual ~EMTest() {

 }
};

TEST_F(EMTest, Resolves3MultimapCase) {
	char offs[] = "33";
	char self[] = "self";
	char * em_argv[] = {self, btfilename, dist_prob_file, reference_fasta, unmapped_fasta, offs};
	EM_main(6, em_argv);
}

}
