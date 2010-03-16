#include "../BowtieEntry.h"
#include <gtest/gtest.h>
#include <fstream>
//#include <string>


namespace {

class BowtieEntryTest : public ::testing::Test {
protected:
	BowtieEntry * bt0;
	BowtieEntry * bt1;
	BowtieEntry * bt2;
	BowtieEntry * bt3;

	ifstream bowtie_stream;

BowtieEntryTest () {

	bt0 = new BowtieEntry(33);
	bt1 = new BowtieEntry("Test BT Entry 1", "Trans1|Gene1|real|X|12345|67890|-1", "+", 123,
						"ACTGACTG", "a~\\tRHYSavrq;slbnqePIgqmw", "9:G>C,19:G>A", 64);
	bt2 = new BowtieEntry(0);
	bt3 = new BowtieEntry(100);

	bowtie_stream.open("../Test/BowtieEntryTest.txt");

	bowtie_stream >> *bt2;
	bowtie_stream >> *bt2;
	bowtie_stream >> *bt3;
}

virtual ~BowtieEntryTest() {
	  delete bt0;
	  delete bt1;
	  delete bt2;
	  delete bt3;
 }
};



TEST_F(BowtieEntryTest, InitializesValues) {
	EXPECT_EQ(NULL, bt0->read);
	EXPECT_EQ(bt0->read_id, "");
	EXPECT_EQ(bt0->base_read_id, "");
	EXPECT_EQ(bt0->mapping, "");
	EXPECT_EQ(bt0->strand, "");
	EXPECT_EQ(bt0->position, -1);
	EXPECT_EQ(bt0->mismatch_indices.size(), 0);
	EXPECT_EQ(bt0->offset, 33);

	EXPECT_EQ(bt1->read_id, "Test BT Entry 1");
	EXPECT_EQ(bt1->base_read_id, "Test BT Entry");
	EXPECT_EQ(bt1->mapping, "Trans1|Gene1|real|X|12345|67890|-1");
	EXPECT_EQ(bt1->strand, "+");
	EXPECT_EQ(bt1->position, 123);
	EXPECT_EQ(bt1->mismatch_indices.size(), 2);
	EXPECT_EQ(bt1->mismatch_indices.at(0), 9);
	EXPECT_EQ(bt1->mismatch_indices.at(1), 19);
	EXPECT_EQ(bt1->offset, 64);
	EXPECT_EQ(bt1->read->sequence, "ACTGACTG");
	EXPECT_EQ(bt1->read->quality->quality_str, "a~\\tRHYSavrq;slbnqePIgqmw");
}

TEST_F(BowtieEntryTest, ReadsFromFileStream) {
	EXPECT_EQ(bt2->read_id, "HWI-EAS76_4:1:100:1000:1131/1");
	EXPECT_EQ(bt2->base_read_id, "HWI-EAS76_4:1:100:1000:1131");
	EXPECT_EQ(bt2->mapping, "ENST00000309052|ENSG00000172954|known|2|30670123|30867091|1");
	EXPECT_EQ(bt2->strand, "-");
	EXPECT_EQ(bt2->position, 4497);
	EXPECT_EQ(bt2->mismatch_indices.size(), 1);
	EXPECT_EQ(bt2->mismatch_indices.at(0), 19);
	EXPECT_EQ(bt2->offset, 0);
	EXPECT_EQ(bt2->read->sequence, "GATGCCCCTGTGTGCATTTAAAATAAAGACTCCAACGGAG");
	EXPECT_EQ(bt2->read->quality->quality_str, "GPDKCDDISSSJSRGSSSSSISSSSSSSSOSSSSSSSSSS");

	EXPECT_EQ(bt3->read_id, "HWI-EAS76_4:1:100:1000:1145/1");
	EXPECT_EQ(bt3->base_read_id, "HWI-EAS76_4:1:100:1000:1145");
	EXPECT_EQ(bt3->mapping, "ENST00000374677|ENSG00000112473|known|6|33168650|33172216|1");
	EXPECT_EQ(bt3->strand, "+");
	EXPECT_EQ(bt3->position, 2212);
	EXPECT_EQ(bt3->mismatch_indices.size(), 0);
	EXPECT_EQ(bt3->offset, 100);
	EXPECT_EQ(bt3->read->sequence, "GGAGGAGTCGGGGATAAACATCAAACATCAATCGTGTTTC");
	EXPECT_EQ(bt3->read->quality->quality_str, "SSSSSSSSSSSSSSSSSSQSSSSSSSSSSSMSSSSHPCOM");
}


TEST_F(BowtieEntryTest, GetsMappedGeneAndTransscript) {
	EXPECT_EQ(bt1->mapped_gene(), "Gene1");
	EXPECT_EQ(bt1->mapped_transcript(), "Trans1");
	EXPECT_EQ(bt2->mapped_gene(), "ENSG00000172954");
	EXPECT_EQ(bt2->mapped_transcript(), "ENST00000309052");
	EXPECT_EQ(bt3->mapped_gene(), "ENSG00000112473");
	EXPECT_EQ(bt3->mapped_transcript(), "ENST00000374677");
}

TEST_F(BowtieEntryTest, CalculatesMappingProbability) {
	EXPECT_DOUBLE_EQ(bt1->mapping_probability(), -0.0000000170595831734277407);
	EXPECT_DOUBLE_EQ(bt2->mapping_probability(), 0.0000000016706221662520599);
	EXPECT_FLOAT_EQ(bt3->mapping_probability(), 3.108508787959149133405226e+72);
}

}

