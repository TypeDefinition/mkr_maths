#pragma once

#include "maths/vector3.h"
#include "maths/matrix.h"

namespace mkr {
    /**
     * @brief
     * Quaternion
     *
     * [https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/index.htm]
     * [https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html]
     *
     * Multiplication of Basis Elements:
     * i^2 = j^2 = k^2 = ijk = -1
     *
     * Quaternion Multiplication Table:
     *   | 1   i   j   k
     * ----------------
     * 1 | 1   i   j   k
     * i | i  -1   k  -j
     * j | j  -k  -1   i
     * k | k   j  -i  -1
     */
    class quaternion {
    public:
        static const quaternion zero;
        static const quaternion identity;

        /**
         * The W component of the quaternion. It is the scalar component.
         * @warning Don't modify this directly unless you know quaternions inside out.
         */
        float w_;
        /**
         * The X component of the quaternion. It is part of the vector component.
         * @warning Don't modify this directly unless you know quaternions inside out.
         */
        float x_;
        /**
         * The Y component of the quaternion. It is part of the vector component.
         * @warning Don't modify this directly unless you know quaternions inside out.
         */
        float y_;
        /**
         * The Z component of the quaternion. It is part of the vector component.
         * @warning Don't modify this directly unless you know quaternions inside out.
         */
        float z_;

        /**
         * Rotate a point around an axis.
         * @param _point The point to rotate.
         * @param _angle The rotation angle.
         * @param _rotation_axis The unit rotation axis.
         * @return The rotated point.
         * @warning _rotation_axis must be a unit vector.
         */
        [[nodiscard]] static vector3 rotate(const vector3& _point, float _angle, const vector3& _rotation_axis);

        /**
         * Rotate a point.
         * @param _point The point to rotate.
         * @param _rotation The unit rotation.
         * @return The rotated point.
         * @warning _rotation must be a rotational quaternion.
         */
        [[nodiscard]] static vector3 rotate(const vector3& _point, const quaternion& _rotation);

        /**
         * Spherical Linear Interpolation between 2 rotational quaternions.
         * @param _start The start rotation.
         * @param _end The end rotation.
         * @param _ratio The ratio to interpolate between _start and _end. A value of 0 will return _start, and a value of 1 will return _end.
         * @param _clamp_ratio If set to true, the _ratio is clamped to between 0 and 1. The default value is false.
         * @return The interpolated rotation.
         * @warning _start and _end must be rotational quaternions.
         */
        [[nodiscard]] static quaternion slerp(const quaternion& _start, const quaternion& _end, float _ratio, bool _clamp_ratio = false);

        /**
         * Constructs a quaternion.
         * @param _w The W component of the quaternion. It is the scalar component.
         * @param _x The X component of the quaternion. It is part of the vector component.
         * @param _y The Y component of the quaternion. It is part of the vector component.
         * @param _z The Z component of the quaternion. It is part of the vector component.
         */
        explicit quaternion(float _w = 1.0f, float _x = 0.0f, float _y = 0.0f, float _z = 0.0f);

        /**
         * Constructs a quaternion.
         * @param _w The W component of the quaternion. It is the scalar component.
         * @param _xyz The XYZ component of the quaternion. It is the vector component.
         */
        explicit quaternion(float _w, const vector3& _xyz);

        /**
         * Constructs a rotational quaternion.
         * @param _rotation_axis The unit rotation axis.
         * @param _angle The rotation angle.
         * @warning _rotation_axis must be a unit vector.
         */
        explicit quaternion(const vector3& _rotation_axis, float _angle);

        bool operator==(const quaternion& _rhs) const;

        bool operator!=(const quaternion& _rhs) const;

        quaternion operator+(const quaternion& _rhs) const;

        quaternion& operator+=(const quaternion& _rhs);

        quaternion operator-(const quaternion& _rhs) const;

        quaternion& operator-=(const quaternion& _rhs);

        quaternion operator*(const quaternion& _rhs) const;

        quaternion& operator*=(const quaternion& _rhs);

        quaternion operator*(float _scalar) const;

        quaternion& operator*=(float _scalar);

        /**
         * Get the dot product of this quaternion and another quaternion.
         * @param _quaternion The other quaternion to calculate the dot product with.
         * @return The dot product of this quaternion and the other quaternion.
         */
        [[nodiscard]] float dot(const quaternion& _quaternion) const;

        /**
         * Normalise this quaternion.
         */
        void normalise();

        /**
         * Get a normalised copy of this quaternion.
         * @return A normalised copy of this quaternion.
         */
        [[nodiscard]] quaternion normalised() const;

        /**
         * Checks if this quaternion is a zero quaternion.
         * @return
         */
        [[nodiscard]] bool is_zero() const;

        /**
         * Checks if this quaternion is a unit quaternion.
         * @return
         */
        [[nodiscard]] bool is_unit() const;

        /**
         * Return the length of the quaternion.
         * @return The length of the quaternion.
         */
        [[nodiscard]] float length() const;

        /**
         * Return the squared length of the quaternion.
         * @return The squared length of the quaternion.
         */
        [[nodiscard]] float length_squared() const;

        /**
         * Set this quaternion to a rotation given an angle in radians and an axis.
         * @param _angle The angle rotation in radians.
         * @param _rotation_axis The axis of rotation.
         * @warning _rotation_axis must be a unit vector.
         */
        void set_rotation(const vector3& _rotation_axis, float _angle);

        /**
         * Conjugate this quaternion.
         * @note For rotational quaternions, the inverse is equal to the conjugate.
         */
        void conjugate();

        /**
         * Get a conjugated copy of this quaternion.
         * @return A conjugated copy of this quaternion.
         * @note For rotational quaternions, the inverse is equal to the conjugate.
         */
        [[nodiscard]] quaternion conjugated() const;

        /**
         * Inverse this quaternion.
         * @note For rotational quaternions, the inverse is equal to the conjugate.
         */
        void inverse();

        /**
         * Get a inversed copy of this quaternion.
         * @return A inversed copy of this quaternion.
         * @note For rotational quaternions, the inverse is equal to the conjugate.
         */
        [[nodiscard]] quaternion inversed() const;

        /**
         * Get this quaternion as a rotation in radians, around an axis.
         * @param _angle The angle of rotation result.
         * @param _rotation_axis The axis of rotation result.
         * @warning This quaternion must be a rotational (unit) quaternion.
         */
        void to_axis_angle(vector3& _rotation_axis, float& _angle) const;

        /**
         * Get this quaternion as a matrix4x4 rotation matrix.
         * @return This quaternion as a matrix4x4 rotation matrix.
         * @warning This function assumes that this quaternion is a rotational (unit) quaternion.
         */
        [[nodiscard]] matrix4x4 to_rotation_matrix() const;

        friend quaternion operator*(float _scalar, const quaternion& _quaternion);

        [[nodiscard]] std::string to_string(const int _precision = 4) const {
            std::ostringstream out;
            out.precision(_precision);
            out << std::fixed;
            out << w_ << "," << x_ << ", " << y_ << ", " << z_;
            return out.str();
        }

        friend std::ostream& operator<<(std::ostream& _stream, const quaternion& _quaternion) {
            return _stream << _quaternion.to_string();
        }
    };
}