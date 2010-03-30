#include <gtest/gtest.h>
#include "../EM_Map.h"

namespace {

class EM_MapTest : public ::testing::Test {
protected:
	EM_Map * em0;
	EM_Map * em1;
	EM_Map * em2;
	EM_Map * em3;

	int2doubleumap d_prob;
	string2intumap i_length;

	ifstream btfilestream;
	BowtieEntry *bt0, *bt1;

	EM_MapTest () {
		btfilestream.open("../Test/EM_MapTest.txt");

		d_prob.insert(int2doubleumap::value_type(100, .1)); //Stick in dummy probability for distance
		i_length.insert(string2intumap::value_type("Trans1|Gene1|known|6|100|1000|1", 400)); //Same for isoform lengths
		i_length.insert(string2intumap::value_type("Trans2|Gene2|known|6|100|1000|1", 400));
		i_length.insert(string2intumap::value_type("TFus12|GFus12|known|2|100|1000|1", 400));

		bt0 = new BowtieEntry(33);
		bt1 = new BowtieEntry(33);

		btfilestream >> *bt0;
		btfilestream >> *bt1;
		em0 = new EM_Map(*bt0, *bt1, d_prob, i_length);

		btfilestream >> *bt0;
		btfilestream >> *bt1;
		em1 = new EM_Map(*bt0, *bt1, d_prob, i_length);

		btfilestream >> *bt0;
		btfilestream >> *bt1;
		em2 = new EM_Map(*bt0, *bt1, d_prob, i_length);

		btfilestream >> *bt0;
		btfilestream >> *bt1;
		em3 = new EM_Map(*bt0, *bt1, d_prob, i_length);

	}

virtual ~EM_MapTest() {
	btfilestream.close();
	delete bt0,
	delete bt1;
	delete em0;
	delete em1;
	delete em2;
 }
};

TEST_F(EM_MapTest, InitializesValues) {
	EXPECT_EQ(em0->base_read_id, "Read1");
//	EXPECT_EQ(em0->start, 0);
//	EXPECT_EQ(em0->end, 100);
	EXPECT_EQ(em0->d, .1);
	EXPECT_EQ(em0->isoform_length, 400);
	EXPECT_EQ(em0->isoform, "Trans1|Gene1|known|6|100|1000|1");

	EXPECT_EQ(em1->base_read_id, "Read1");
//	EXPECT_EQ(em1->start, 50);
//	EXPECT_EQ(em1->end, 150);
	EXPECT_EQ(em1->d, .1);
	EXPECT_EQ(em1->isoform_length, 400);
	EXPECT_EQ(em1->isoform, "TFus12|GFus12|known|2|100|1000|1");

	EXPECT_EQ(em2->base_read_id, "Read2");
//	EXPECT_EQ(em2->start, 0);
//	EXPECT_EQ(em2->end, 100);
	EXPECT_EQ(em2->d, .1);
	EXPECT_EQ(em2->isoform_length, 400);
	EXPECT_EQ(em2->isoform, "Trans1|Gene1|known|6|100|1000|1");

}

TEST_F(EM_MapTest, CalculatesP) {
	EXPECT_DOUBLE_EQ(em0->P, 0.9980127120145297903164305);
	EXPECT_DOUBLE_EQ(em1->P, 0.9980127120145297903164305);
	EXPECT_FLOAT_EQ(em2->P, 0.9953424204896054172664321);
	EXPECT_DOUBLE_EQ(em3->P, 0.9966766719733826107585628);

}

TEST_F(EM_MapTest, CalculatesEMProb) {
	EXPECT_FLOAT_EQ(em0->em_prob(.05), 0.0000124751589001816224928);
	EXPECT_FLOAT_EQ(em1->em_prob(.2), 0.0000499006356007264899712);
	EXPECT_FLOAT_EQ(em2->em_prob(.001), 0.0000002488356051224013319);
	EXPECT_FLOAT_EQ(em3->em_prob(.05), 0.0000124584583996625006714);
}
}
