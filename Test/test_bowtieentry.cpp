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

BowtieEntryTest () {

	bt0 = new BowtieEntry(33);
	bt1 = new BowtieEntry("Test BT Entry 1", "Trans1|Gene1|real|X|12345|67890|-1", "+", 123,
						"ACTGACTG", "a~\tRHYS", "9:G>C,19:G>A", 64);
	bt2 = new BowtieEntry(0);
	bt3 = new BowtieEntry(100);


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
}

TEST_F(BowtieEntryTest, ReadsFromFileStream) {

}


TEST_F(BowtieEntryTest, GetsMappedGeneAndTransscript) {

}

TEST_F(BowtieEntryTest, CalculatesMappingProbability) {

}

}

