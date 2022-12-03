#include <gtest/gtest.h>
#include "maths/matrix.h"
#include "maths/quaternion.h"
#include "maths/matrix_util.h"

using namespace mkr;

TEST(quaternion_test, conjugate) {
    {
        quaternion a{4.0f, 23.0f, 6.0f, 7.0f};
        quaternion b{4.0f, -23.0f, -6.0f, -7.0f};
        EXPECT_FALSE(a == b);
        EXPECT_TRUE(a.conjugated() == b);
        EXPECT_FALSE(a.conjugated() != b);

        a.conjugate();
        EXPECT_FALSE(a.conjugated() == b);
        EXPECT_FALSE(a != b);
        EXPECT_TRUE(a == b);
    }
}

TEST(quaternion_test, inverse) {
    {
        quaternion a{4.0f, 23.0f, 6.0f, 8.0f};
        quaternion b{4.0f / 645.0f, -23.0f / 645.0f, -6.0f / 645.0f, -8.0f / 645.0f};
        EXPECT_FALSE(a == b);
        EXPECT_TRUE(a.inversed() == b);
        EXPECT_FALSE(a.inversed() != b);

        a.inverse();
        EXPECT_FALSE(a.inversed() == b);
        EXPECT_FALSE(a != b);
        EXPECT_TRUE(a == b);
    }
}

TEST(quaternion_test, rotate) {
    {
        quaternion q{vector3::x_axis, 45.0f * maths_util::deg2rad};
        matrix4x4 m = matrix_util::rotation_matrix({45.0f * maths_util::deg2rad, 0.0f, 0.0f});
        EXPECT_TRUE(q.to_rotation_matrix() == m);
    }

    {
        quaternion q{vector3::y_axis, 120.0f * maths_util::deg2rad};
        matrix4x4 m = matrix_util::rotation_matrix({0.0f, 120.0f * maths_util::deg2rad, 0.0f});
        EXPECT_TRUE(q.to_rotation_matrix() == m);
    }

    {
        quaternion q{vector3::z_axis, 943.0f * maths_util::deg2rad};
        matrix4x4 m = matrix_util::rotation_matrix({0.0f, 0.0f, 943.0f * maths_util::deg2rad});
        EXPECT_TRUE(q.to_rotation_matrix() == m);
    }

    {
        quaternion qx{vector3::x_axis, 123.0f * maths_util::deg2rad};
        quaternion qy{vector3::y_axis, 456.0f * maths_util::deg2rad};
        quaternion qz{vector3::z_axis, 789.0f * maths_util::deg2rad};
        quaternion q = qx * qy * qz;
        matrix4x4 m = matrix_util::rotation_matrix(vector3{123.0f, 456.0f, 789.0f} * maths_util::deg2rad);
        EXPECT_TRUE(q.to_rotation_matrix() == m);
    }

    {
        quaternion q{vector3::x_axis, 0.0f};
        vector3 axis;
        float angle;
        q.to_axis_angle(axis, angle);
        EXPECT_TRUE(0.0f == angle);
        EXPECT_TRUE(vector3::x_axis == axis);
    }

    {
        quaternion q{vector3::y_axis, 0.0f};
        vector3 axis;
        float angle;
        q.to_axis_angle(axis, angle);
        EXPECT_TRUE(0.0f == angle);
        EXPECT_TRUE(vector3::x_axis == axis);
    }

    {
        quaternion q{vector3::z_axis, 0.0f};
        vector3 axis;
        float angle;
        q.to_axis_angle(axis, angle);
        EXPECT_TRUE(0.0f == angle);
        EXPECT_TRUE(vector3::x_axis == axis);
    }

    {
        quaternion q{-vector3::x_axis * 12.0f, 6.0f * maths_util::pi};
        vector3 axis;
        float angle;
        q.normalised().to_axis_angle(axis, angle);
        EXPECT_TRUE(0.0f == angle);
        EXPECT_TRUE(vector3::x_axis == axis);
    }

    {
        quaternion q{-vector3::y_axis, 2.0f * maths_util::pi};
        vector3 axis;
        float angle;
        q.to_axis_angle(axis, angle);
        EXPECT_TRUE(0.0f == angle);
        EXPECT_TRUE(vector3::x_axis == axis);
    }

    {
        quaternion q{-vector3::z_axis * 37.0f, 7.0f * maths_util::pi};
        q.normalise();
        vector3 axis; float angle;
        q.to_axis_angle(axis, angle);
        EXPECT_TRUE((maths_util::pi - angle) < 0.0001f); // Hardcode a larger epsilon due to floating point error.
        EXPECT_TRUE(vector3::z_axis == axis);
    }
}