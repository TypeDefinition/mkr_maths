#include <gtest/gtest.h>
#include "maths/matrix.h"
#include "maths/matrix_util.h"

using namespace mkr;

TEST(matrix, zero) {
    EXPECT_TRUE(matrix1x1{} == matrix1x1::zero());
    EXPECT_TRUE(matrix1x2{} == matrix1x2::zero());
    EXPECT_TRUE(matrix1x3{} == matrix1x3::zero());
    EXPECT_TRUE(matrix1x4{} == matrix1x4::zero());

    EXPECT_TRUE(matrix2x1{} == matrix2x1::zero());
    EXPECT_TRUE(matrix2x2{} == matrix2x2::zero());
    EXPECT_TRUE(matrix2x3{} == matrix2x3::zero());
    EXPECT_TRUE(matrix2x4{} == matrix2x4::zero());

    EXPECT_TRUE(matrix3x1{} == matrix3x1::zero());
    EXPECT_TRUE(matrix3x2{} == matrix3x2::zero());
    EXPECT_TRUE(matrix3x3{} == matrix3x3::zero());
    EXPECT_TRUE(matrix3x4{} == matrix3x4::zero());

    EXPECT_TRUE(matrix4x1{} == matrix4x1::zero());
    EXPECT_TRUE(matrix4x2{} == matrix4x2::zero());
    EXPECT_TRUE(matrix4x3{} == matrix4x3::zero());
    EXPECT_TRUE(matrix4x4{} == matrix4x4::zero());
}

TEST(matrix, equality) {
    EXPECT_TRUE(matrix1x1{} != matrix1x1::identity());
    EXPECT_TRUE(matrix2x2{} != matrix2x2::identity());
    EXPECT_TRUE(matrix3x3{} != matrix3x3::identity());
    EXPECT_TRUE(matrix4x4{} != matrix4x4::identity());

    EXPECT_FALSE(matrix1x1{} == matrix1x1::identity());
    EXPECT_FALSE(matrix2x2{} == matrix2x2::identity());
    EXPECT_FALSE(matrix3x3{} == matrix3x3::identity());
    EXPECT_FALSE(matrix4x4{} == matrix4x4::identity());
}

TEST(matrix_test, add) {
    {
        matrix2x4 a{{2.0f, 4.0f, 7.0f, 8.0f,
                     2.0f, 1.0f, 0.0f, 15.0f}};
        matrix2x4 b{{15.0f, 2.0f, 1.0f, 0.0f,
                     0.0f, 1.0f, 77.0f, 1045.0f}};
        matrix2x4 c{{17.0f, 6.0f, 8.0f, 8.0f,
                     2.0f, 2.0f, 77.0f, 1060.0f}};
        EXPECT_TRUE(c == a + b);
    }
}

TEST(matrix_test, sub) {
    {
        matrix2x4 a{{2.0f, 4.0f, 7.0f, 8.0f,
                     2.0f, 1.0f, 0.0f, 15.0f}};
        matrix2x4 b{{15.0f, 2.0f, 1.0f, 0.0f,
                     0.0f, 1.0f, 77.0f, 1045.0f}};
        matrix2x4 c{{17.0f, 6.0f, 8.0f, 8.0f,
                     2.0f, 2.0f, 77.0f, 1060.0f}};
        EXPECT_TRUE(c - b == a);
    }
}

TEST(matrix_test, mult) {
    {
        matrix1x1 a{{4.0f}};
        matrix1x1 b{{16.0f}};
        EXPECT_TRUE(a * a == b);
        EXPECT_FALSE(a * a == a);
    }

    {
        matrix4x1 a{{1.0f, 2.0f, 3.0f, 4.0f}};
        matrix1x4 b{{1.0f, 2.0f, 3.0f, 4.0f}};
        matrix1x1 c{{30.0f}};
        EXPECT_TRUE(a * b == c);
    }

    {
        matrix1x3 a{{88.0f, 45.0f, -12.0f}};
        matrix3x1 b{{-15.0f, -103.0f, 9.0f}};
        matrix3x3 c{{-1320.0f, -675.0f, 180.0f,
                     -9064.0f, -4635.0f, 1236.0f,
                     792.0f, 405.0f, -108.0f}};
        EXPECT_TRUE(a * b == c);
    }

    {
        matrix4x4 a{{1.0f, 2.0f, 3.0f, 4.0f,
                     5.0f, 6.0f, 7.0f, 8.0f,
                     9.0f, 10.0f, 11.0f, 12.0f,
                     13.0f, 14.0f, 15.0f, 16.0f}};
        matrix4x4 b{{17.0f, 18.0f, 19.0f, 20.0f,
                     21.0f, 22.0f, 23.0f, 24.0f,
                     25.0f, 26.0f, 27.0f, 28.0f,
                     29.0f, 30.0f, 31.0f, 32.0f}};
        matrix4x4 c{{538.0f, 612.0f, 686.0f, 760.0f,
                     650.0f, 740.0f, 830.0f, 920.0f,
                     762.0f, 868.0f, 974.0f, 1080.0f,
                     874.0f, 996.0f, 1118.0f, 1240.0f}};
        EXPECT_TRUE(a * b == c);
        EXPECT_TRUE(c == a * b);
        EXPECT_TRUE(a * matrix4x4::identity() == a);
        EXPECT_TRUE(b * matrix4x4::identity() == b);
        EXPECT_TRUE(a * matrix4x4::identity() != b);
        EXPECT_FALSE(c * matrix4x4::identity() == b);
    }
}

TEST(matrix_test, transpose) {
    {
        matrix3x2 a{{2.0f, 4.0f, 6.0f,
                     8.0f, 10.0f, 12.0f}};
        matrix2x3 b{{2.0f, 6.0f, 10.0f,
                     4.0f, 8.0f, 12.0f}};
        EXPECT_TRUE(a.transposed() == b);
    }

    {
        EXPECT_TRUE((matrix<10, 10>::identity() == matrix<10, 10>::identity().transposed()));
    }
}

TEST(matrix_test, inverse) {
    {
        EXPECT_TRUE(matrix_util::invertible(matrix1x1::diagonal(4.0f)));
        EXPECT_TRUE(matrix_util::invertible(matrix2x2::diagonal(5.0f)));
        EXPECT_TRUE(matrix_util::invertible(matrix3x3::diagonal(6.0f)));
        EXPECT_TRUE(matrix_util::invertible(matrix4x4::identity()));

        matrix4x4 a{{1.0f, 1.0f, 1.0f, 1.0f,
                     2.0f, 2.0f, 2.0f, 2.0f,
                     3.0f, 3.0f, 3.0f, 3.0f,
                     4.0f, 4.0f, 4.0f, 4.0f}};
        EXPECT_FALSE(matrix_util::invertible(a));
    }

    {
        matrix1x1 a{{27.0f}};
        matrix1x1 b{{1.0f / 27.0f}};
        EXPECT_TRUE(matrix_util::inverse_matrix(a) == b);
        EXPECT_TRUE(matrix_util::inverse_matrix(matrix1x1::identity()) == matrix1x1::identity());
    }

    {
        matrix2x2 a{{2.0f, 7.0f,
                     4.0f, 3.0f}};
        matrix2x2 b{{-3.0f / 22.0f, 7.0f / 22.0f,
                     2.0f / 11.0f, -1.0f / 11.0f}};
        EXPECT_TRUE(matrix_util::inverse_matrix(a) == b);
        EXPECT_TRUE(matrix_util::inverse_matrix(matrix2x2::identity()) == matrix2x2::identity());
    }

    {
        matrix3x3 a{{1.0f, 4.0f, 7.0f,
                     2.0f, 5.0f, 2.0f,
                     3.0f, 6.0f, 9.0f}};
        matrix3x3 b{{-11.0f / 12.0f, -1.0f / 6.0f, 3.0f / 4.0f,
                     1.0f / 3.0f, 1.0f / 3.0f, -1.0f / 3.0f,
                     1.0f / 12.0f, -1.0f / 6.0f, 1 / 12.0f}};
        EXPECT_TRUE(matrix_util::inverse_matrix(a) == b);
        EXPECT_TRUE(matrix_util::inverse_matrix(matrix3x3::identity()) == matrix3x3::identity());
    }

    {
        matrix4x4 a{{2.0f, 8.0f, 3.0f, 0.0f,
                     3.0f, 12.0f, 3.0f, 5.0f,
                     7.0f, 23.0f, 12.0f, 11.0f,
                     11.0f, -5.0f, -4.0f, 25.0f}};
        matrix4x4 b{{1039.0f / 1501.0f, -431.0f / 1501.0f, -215.0f / 3002.0f, 267.0f / 3002.0f,
                     -35.0f / 1501.0f, 263.0f / 1501.0f, -65.0f / 1501.0f, -24.0f / 1501.0f,
                     -99.0f / 1501.0f, -414.0f / 1501.0f, 245.0f / 1501.0f, -25.0f / 1501.0f,
                     -480.0f / 1501.0f, 176.0f / 1501.0f, 147.0f / 3002.0f, -15.0f / 3002.0f}};
        EXPECT_TRUE(matrix_util::inverse_matrix(a).value() == b);
        EXPECT_TRUE(matrix_util::inverse_matrix(matrix4x4::identity()) == matrix4x4::identity());
    }
}