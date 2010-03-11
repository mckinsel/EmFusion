#include "../main.h"
#include "../Quality.h"
#include <gtest/gtest.h>

namespace {

// The fixture for testing class Foo.
class QualityTest : public ::testing::Test {
 protected:
	Quality * q1;
	Quality * q2;
	Quality * q3;


  QualityTest() {
	  q1 = new Quality("ABCDEFGHIJKLMNOP", 0);
	  q2 = new Quality("ABCDEFGHIJKLMNOP", 33);
	  q3 = new Quality("\t \"~\r`", -10);
  }

  virtual ~QualityTest() {
	  delete q1;
	  delete q2;
	  delete q3;
  }

};
}

TEST_F(QualityTest, ConvertsCharToProb) {
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(0), 0.0000003162277660168379191);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(1), 0.0000002511886431509582273);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(2), 0.0000001995262314968878707);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(3), 0.0000001584893192461114080);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(4), 0.0000001258925411794166167);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(5), 0.0000000999999999999999955);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(6), 0.0000000794328234724282196);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(7), 0.0000000630957344480192964);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(8), 0.0000000501187233627272497);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(9), 0.0000000398107170553496896);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(10), 0.0000000316227766016837919);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(11), 0.0000000251188643150958207);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(12), 0.0000000199526231496887864);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(13), 0.0000000158489319246111428);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(14), 0.0000000125892541179416610);
	EXPECT_DOUBLE_EQ(q1->error_probabilities.at(15), 0.0000000100000000000000002);
}

TEST_F(QualityTest, HandlesOffset) {
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(0), 0.00063095734448019298);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(1), 0.00050118723362727253);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(2), 0.00039810717055349735);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(3), 0.00031622776601683794);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(4), 0.00025118864315095795);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(5), 0.00019952623149688788);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(6), 0.00015848931924611142);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(7), 0.00012589254117941674);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(8), 0.00010000000000000000);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(9), 0.00007943282347242822);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(10), 0.00006309573444801929);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(11), 0.00005011872336272725);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(12), 0.00003981071705534969);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(13), 0.00003162277660168380);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(14), 0.00002511886431509582);
	EXPECT_DOUBLE_EQ(q2->error_probabilities.at(15), 0.00001995262314968879);
}

TEST_F(QualityTest, IgnoresEscapes) {
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(0), 0.000000000063095734448019424597);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(1), 0.000000000000251188643150958201);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(2), 0.000063095734448019292995787910);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(3), 0.000039810717055349694576925346);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(4), 0.000000000000025118864315095823);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(5), 0.000000000063095734448019424597);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(6), 0.000000000000398107170553496920);
	EXPECT_DOUBLE_EQ(q3->error_probabilities.at(7), 0.000000000025118864315095823148);
}

int Test_main(int argc, char * argv[]) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
