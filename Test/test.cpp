//#include "test.h"
#include <gtest/gtest.h>
#include <fstream>

int Test_main(int argc, char * argv[]) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
