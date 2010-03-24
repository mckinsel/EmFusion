#include <gtest/gtest.h>
#include <fstream>
#include "../main.h"

using namespace std;

namespace {

TEST(SiftTest, Works) {
	char btfilename1[50], btfilename2[50];
	char fastqfilename1[50], fastqfilename2[50];
	char offs[] = "33";
	char repeatfilename1[50], repeatfilename2[50];
	char self[] = "self";
	strcpy(btfilename1, "../Test/SiftTest_BT1.txt");
	strcpy(btfilename2, "../Test/SiftTest_BT2.txt");
	strcpy(fastqfilename1, "../Test/SiftTest_FQ1.txt");
	strcpy(fastqfilename2, "../Test/SiftTest_FQ2.txt");
	strcpy(repeatfilename1, "../Test/SiftTest_Repeat1.txt");
	strcpy(repeatfilename2, "../Test/SiftTest_Repeat2.txt");

	char * sift_argv[] = {self, btfilename1, btfilename2, fastqfilename1, fastqfilename2, offs,
							repeatfilename1,repeatfilename2};

	Sift_main(8, sift_argv);

	ifstream unmapped_stream("../Test/SiftTest_FQ1.txt.unmapped.fasta");
	string nextline;

	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, ">Read7/1");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "TCAGGCGGCCCAAGACACTGCGACTCCGGAGGCAGCCCAG");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, ">Read7/2");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "TCAGGCGGCCCAAGACACTGCGACTCCGGAGGCAGCCCAG");

	unmapped_stream.close();

	unmapped_stream.open("../Test/SiftTest_FQ1.txt.mapdist");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "80");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "360");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "100");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "110");
	unmapped_stream >> nextline;
	EXPECT_EQ(nextline, "120");

	unmapped_stream.close();

	char nl[201];
	unmapped_stream.open("../Test/SiftTest_FQ1.txt.discord");
	unmapped_stream.getline(nl, 200);
	EXPECT_EQ(string(nl), "Trans4\tTrans5\tGene3\tGene4\t400\t300");
	unmapped_stream.getline(nl, 200);
	EXPECT_EQ(string(nl), "Trans4\tTrans7\tGene3\tGene6\t400\t300");
	unmapped_stream.getline(nl, 200);
	EXPECT_EQ(string(nl), "Trans5\tTrans6\tGene4\tGene5\t300\t400");
	unmapped_stream.getline(nl, 200);
	EXPECT_EQ(string(nl), "Trans6\tTrans7\tGene5\tGene6\t400\t300");


	unmapped_stream.close();
}

}
