
#include <string>
#include "hapSamp.hpp"
#include "gtest/gtest.h"

using namespace std;

TEST(HapSampTest, small) {

  string hapFile = "../../samples/test10/test10.madeUpData1.hap";
  string sampFile = "../../samples/test10/test10.madeUpData1.sample";

  HapSamp chunk1(hapFile, sampFile);
}
