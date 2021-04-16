#include "../CppUnitTestFramework.hpp"
#include "../src/PixelDifferenceMethods.cpp"

namespace
{
struct PixelDifferenceMethodsTests
{
};
} // namespace

TEST_CASE(PixelDifferenceMethodsTests, SubtractZeroMatrix)
{
    Eigen::Matrix3d original;
    original << 0, 0, 0,
        0, 0, 0,
        0, 0, 0;
    Eigen::Matrix3d result;
    result << 0, 0, 0,
        0, 0, 0,
        0, 0, 0;
    Eigen::Matrix3d difference;
    difference = original - result;
    REQUIRE(difference.isZero());
}
