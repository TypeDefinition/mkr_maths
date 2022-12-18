#include "quaternion.h"

namespace mkr {
    const quaternion quaternion::zero = quaternion{0.0f, 0.0f, 0.0f, 0.0f};
    const quaternion quaternion::identity = quaternion{1.0f, 0.0f, 0.0f, 0.0f};

    vector3 quaternion::rotate(const vector3& _point, float _angle, const vector3& _rotation_axis) {
        /**
         * Rotation Formula:
         * R * P * R.Inverse
         */
        quaternion rotation{_rotation_axis, _angle};
        quaternion point{0.0f, _point.x_, _point.y_, _point.z_};
        quaternion result = rotation * point * rotation.conjugated();
        return vector3{result.x_, result.y_, result.z_};
    }

    quaternion quaternion::slerp(const quaternion& _start, const quaternion& _end, float _ratio, bool _clamp_ratio) {
        const float ratio = _clamp_ratio ? maths_util::clamp(_ratio, 0.0f, 1.0f) : _ratio;

        /**
         * Let start quaternion be S.
         * Let end quaternion be E.
         * Let the quaternion needed to transform S to E be D (difference).
         * We need to find what is D.
         *
         * D * S = E
         * Therefore,
         * D * S * S.Inverse = E * s.Inverse
         * D = E * S.Inverse (S * S.Inverse cancels each other out)
         */
        quaternion d = _end * _start.conjugated();

        /**
          * Since S * D = E, if we want for example to only turn halfway, then we need to half the angle that d turns.
          * For example, D could represent a rotation of 70 degrees around the Vector::UP axis.
          * So to do a half rotation, we need to do a rotation of 35 degrees around the Vector::UP axis.
          * To do that, we must first figure out what the angle and rotational axis of D is.
          */
        vector3 d_axis;
        float d_angle;
        d.to_axis_angle(d_axis, d_angle);

        // Now that we have the angle and axis to rotate around, we have our result.
        return quaternion{d_axis, d_angle * ratio} * _start;
    }

    quaternion::quaternion(float _w, float _x, float _y, float _z)
            : w_(_w), x_(_x), y_(_y), z_(_z) {}

    quaternion::quaternion(float _w, const vector3& _xyz)
            : w_(_w), x_(_xyz.x_), y_(_xyz.y_), z_(_xyz.z_) {}

    quaternion::quaternion(const vector3& _rotation_axis, float _angle) {
        set_rotation(_rotation_axis, _angle);
    }

    bool quaternion::operator==(const quaternion& _rhs) const {
        return maths_util::approx_equal(w_, _rhs.w_) &&
               maths_util::approx_equal(x_, _rhs.x_) &&
               maths_util::approx_equal(y_, _rhs.y_) &&
               maths_util::approx_equal(z_, _rhs.z_);
    }

    bool quaternion::operator!=(const quaternion& _rhs) const {
        return !((*this) == _rhs);
    }

    quaternion quaternion::operator+(const quaternion& _rhs) const {
        return quaternion{w_ + _rhs.w_, x_ + _rhs.x_, y_ + _rhs.y_, z_ + _rhs.z_};
    }

    quaternion& quaternion::operator+=(const quaternion& _rhs) {
        *this = (*this) + _rhs;
        return *this;
    }

    quaternion quaternion::operator-(const quaternion& _rhs) const {
        return quaternion{w_ - _rhs.w_, x_ - _rhs.x_, y_ - _rhs.y_, z_ - _rhs.z_};
    }

    quaternion& quaternion::operator-=(const quaternion& _rhs) {
        *this = (*this) - _rhs;
        return *this;
    }

    quaternion quaternion::operator*(const quaternion& _rhs) const {
        const vector3 this_vector{x_, y_, z_};
        const vector3 rhs_vector{_rhs.x_, _rhs.y_, _rhs.z_};

        /**
         * Quaternion Multiplication Formula:
         * (sa,va) * (sb,vb) = (sa*sb-va•vb, va×vb + sa*vb + sb*va)
         */
        const float w = w_ * _rhs.w_ - this_vector.dot(rhs_vector);
        const vector3 xyz = w_ * rhs_vector + _rhs.w_ * this_vector + this_vector.cross(rhs_vector);

        return quaternion{w, xyz};
    }

    quaternion& quaternion::operator*=(const quaternion& _rhs) {
        *this = (*this) * _rhs;
        return *this;
    }

    quaternion quaternion::operator*(float _scalar) const {
        return quaternion{w_ * _scalar, x_ * _scalar, y_ * _scalar, z_ * _scalar};
    }

    quaternion& quaternion::operator*=(float _scalar) {
        *this = (*this) * _scalar;
        return *this;
    }

    float quaternion::dot(const quaternion& _quaternion) const {
        return (w_ * _quaternion.w_) + (x_ * _quaternion.x_) + (y_ * _quaternion.y_) + (z_ * _quaternion.z_);
    }

    void quaternion::normalise() {
        const float length = this->length();
        if (maths_util::approx_equal(length, 0.0f)) {
            w_ = x_ = y_ = z_ = 0.0f;
            return;
        }
        w_ /= length;
        x_ /= length;
        y_ /= length;
        z_ /= length;
    }

    quaternion quaternion::normalised() const {
        const float length = this->length();
        if (maths_util::approx_equal(length, 0.0f)) {
            return quaternion::zero;
        }
        return quaternion{w_ / length, x_ / length, y_ / length, z_ / length};
    }

    bool quaternion::is_zero() const {
        return maths_util::approx_equal(0.0f, length_squared());
    }

    bool quaternion::is_unit() const {
        return maths_util::approx_equal(1.0f, length_squared());
    }

    float quaternion::length() const {
        return std::sqrt(length_squared());
    }

    float quaternion::length_squared() const {
        return (w_ * w_) + (x_ * x_) + (y_ * y_) + (z_ * z_);
    }

    void quaternion::conjugate() {
        x_ = -x_;
        y_ = -y_;
        z_ = -z_;
    }

    quaternion quaternion::conjugated() const {
        return quaternion{w_, -x_, -y_, -z_};
    }

    void quaternion::inverse() {
        (*this) = inversed();
    }

    quaternion quaternion::inversed() const {
        return 1 / length_squared() * conjugated();
    }

    void quaternion::set_rotation(const vector3& _rotation_axis, float _angle) {
        /**
         * Axis-Angle Rotation To Quaternion:
         * Quaternion = cos(a/2) + [x*sin(a/2)]i + [y*sin(a/2)]j + [z*sin(a/2)]k
         */
        vector3 xyz = _rotation_axis * std::sin(_angle * 0.5f);
        x_ = xyz.x_;
        y_ = xyz.y_;
        z_ = xyz.z_;
        w_ = std::cos(_angle * 0.5f);
    }

    void quaternion::to_axis_angle(vector3& _rotation_axis, float& _angle) const {
        vector3 xyz{x_, y_, z_};
        if (xyz.is_zero()) {
            _angle = 0.0f;
            _rotation_axis = vector3::x_axis;
            return;
        }

        _angle = std::acos(w_) * 2.0f;
        _rotation_axis = xyz.normalised();
    }

    matrix4x4 quaternion::to_rotation_matrix() const {
        return matrix4x4{{1.0f - 2.0f * y_ * y_ - 2.0f * z_ * z_, 2.0f * x_ * y_ + 2.0f * w_ * z_, 2.0f * x_ * z_ - 2.0f * w_ * y_, 0.0f,
                          2.0f * x_ * y_ - 2.0f * w_ * z_, 1.0f - 2.0f * x_ * x_ - 2.0f * z_ * z_, 2.0f * y_ * z_ + 2.0f * w_ * x_, 0.0f,
                          2.0f * x_ * z_ + 2.0f * w_ * y_, 2.0f * y_ * z_ - 2.0f * w_ * x_, 1.0f - 2.0f * x_ * x_ - 2.0f * y_ * y_, 0.0f,
                          0.0f, 0.0f, 0.0f, 1.0f}};
    }

    quaternion operator*(float _scalar, const quaternion& _quaternion) {
        return _quaternion * _scalar;
    }
}