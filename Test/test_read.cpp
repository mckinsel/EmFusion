#include <gtest/gtest.h>
#include <fstream>
#include "../Read.h"

namespace {
class ReadTest : public ::testing::Test {
 protected:
	Read * r0;
	Read * r1;
	Read * r2;
	Read * r3;
	Read * r4;
	Read * r5;

	ifstream read_test_stream;

	ofstream read_output_stream;

	string r41id, r41seq, r41qual;
	string r42id, r42seq, r42qual;

	string fastaline1, fastaline2;
	string fastqline1, fastqline2, fastqline3, fastqline4;

  ReadTest() {
	  r0 = new Read(50);
	  r1 = new Read("Test Read 1", "ACTGACTGACGCGATGTCAGT", "TTSSZZSSttSSQQSSSAASS", 33);
	  r2 = new Read("TEST Read\\2", "ACTG", "SSSSSSSS", 64);
	  r3 = new Read(33);
	  r4 = new Read(64);
	  r5 = new Read(0);

	  read_test_stream.open("../Test/ReadTest.txt");
	  read_output_stream.open("../Test/read_test_output.txt");

	  read_test_stream >> *r3;
	  read_test_stream >> *r4;

	  r41id = r4->id;
	  r41seq = r4->sequence;
	  r41qual = r4->quality->quality_str;

	  read_test_stream >> *r4;

	  r42id = r4->id;
	  r42seq = r4->sequence;
	  r42qual = r4->quality->quality_str;

	  read_test_stream >> *r5;

	  r1->write_as_fasta(read_output_stream);
	  r5->write_as_fastq(read_output_stream);

	  read_output_stream.close();
	  read_test_stream.close();

	  char streamchar[150];
	  ifstream teststream("../Test/read_test_output.txt");
	  teststream.getline(streamchar, (streamsize)150);
	  fastaline1 = string(streamchar);
	  teststream.getline(streamchar, (streamsize)150);
	  fastaline2 = string(streamchar);
	  teststream.getline(streamchar, (streamsize)150);
	  fastqline1 = string(streamchar);
	  teststream.getline(streamchar, (streamsize)150);
	  fastqline2 = string(streamchar);
	  teststream.getline(streamchar, (streamsize)150);
	  fastqline3 = string(streamchar);
	  teststream.getline(streamchar, (streamsize)150);
	  fastqline4 = string(streamchar);

  }

  virtual ~ReadTest() {
	  delete r0;
	  delete r1;
	  delete r2;
	  delete r3;
	  delete r4;
	  delete r5;
  }

};

TEST_F(ReadTest, InitializesValues) {
	EXPECT_EQ(r0->id, "");
	EXPECT_EQ(r0->base_id, "");
	EXPECT_EQ(r0->sequence, "");
	EXPECT_EQ(NULL, r0->quality);
	EXPECT_EQ(r0->offset, 50);

	EXPECT_EQ(r1->id, "Test Read 1");
	EXPECT_EQ(r1->base_id, "Test Read");
	EXPECT_EQ(r1->sequence, "ACTGACTGACGCGATGTCAGT");
	EXPECT_EQ(r1->quality->quality_str, "TTSSZZSSttSSQQSSSAASS");
	EXPECT_EQ(r1->offset, 33);

	EXPECT_EQ(r2->id, "TEST Read\\2");
	EXPECT_EQ(r2->base_id, "TEST Read");
	EXPECT_EQ(r2->sequence, "ACTG");
	EXPECT_EQ(r2->quality->quality_str, "SSSSSSSS");
	EXPECT_EQ(r2->offset, 64);
}

TEST_F(ReadTest, ReadsFromFileStream) {
	EXPECT_EQ(r3->id, "R:e\\dtes1ing3\\1");
	EXPECT_EQ(r3->id.length(), 15);
	EXPECT_EQ(r3->base_id, "R:e\\dtes1ing3");
	EXPECT_EQ(r3->sequence, "ACTGACTGACTGACTGACTG");
	EXPECT_EQ(r3->quality->quality_str, "QQQQQQ@@@@@@aaaaaaaa");

	EXPECT_EQ(r4->id, "R:e\\dtes1ing5:1");
	EXPECT_EQ(r4->id.length(), 15);
	EXPECT_EQ(r4->base_id, "R:e\\dtes1ing5");
	EXPECT_EQ(r4->sequence, "ACTGACTGACTGACTGACTG");
	EXPECT_EQ(r4->quality->quality_str, "\\t\\t\\r\\n\\\"\\n\\r\\txxxx");
	EXPECT_EQ(r4->quality->quality_str.length(), 20);
}

TEST_F(ReadTest, ReturnsValues) {
	EXPECT_EQ(r41id, "R:e\\dte@1i ng4\\4");
	EXPECT_EQ(r41seq, "ACTGACTGGGGGACTGACTG");
	EXPECT_EQ(r41qual, "QQAAQQ@@@@@@aaazzaaa");
}

TEST_F(ReadTest, WritesFASTA) {
	EXPECT_EQ(fastaline1, ">Test Read 1");
	EXPECT_EQ(fastaline2, "ACTGACTGACGCGATGTCAGT");
}

TEST_F(ReadTest, WritesFASTQ) {
	EXPECT_EQ(fastqline1, "@R:e\\dtes1ing5:1");
	EXPECT_EQ(fastqline2, "ACTGACTGACTGACTGACTG");
	EXPECT_EQ(fastqline3, "+");
	EXPECT_EQ(fastqline4, "\\t\\t\\r\\n\\\"\\n\\r\\txxxx");
}
}
