#include "maths/vector2.h"

namespace mkr {
    const vector2 vector2::zero{0.0f, 0.0f};
    const vector2 vector2::up{0.0f, 1.0f};
    const vector2 vector2::down{0.0f, -1.0f};
    const vector2 vector2::left{1.0f, 0.0f};
    const vector2 vector2::right{-1.0f, 0.0f};
    const vector2 vector2::x_direction{1.0f, 0.0f};
    const vector2 vector2::y_direction{0.0f, 1.0f};

    vector2::vector2(float _x, float _y)
            : x_(_x), y_(_y) {}

    vector2 vector2::operator-() const {
        return vector2(-x_, -y_);
    }

    bool vector2::operator==(const vector2 &_rhs) const {
        return maths_util::approx_equal(x_, _rhs.x_) &&
               maths_util::approx_equal(y_, _rhs.y_);
    }

    bool vector2::operator!=(const vector2 &_rhs) const {
        return !(*this == _rhs);
    }

    vector2 vector2::operator+(const vector2 &_rhs) const {
        return vector2(x_ + _rhs.x_, y_ + _rhs.y_);
    }

    vector2 &vector2::operator+=(const vector2 &_rhs) {
        x_ += _rhs.x_;
        y_ += _rhs.y_;
        return *this;
    }

    vector2 vector2::operator-(const vector2 &_rhs) const {
        return vector2(x_ - _rhs.x_, y_ - _rhs.y_);
    }

    vector2 &vector2::operator-=(const vector2 &_rhs) {
        x_ -= _rhs.x_;
        y_ -= _rhs.y_;
        return *this;
    }

    vector2 vector2::operator*(float _scalar) const {
        return vector2(_scalar * x_, _scalar * y_);
    }

    vector2 &vector2::operator*=(float _scalar) {
        x_ *= _scalar;
        y_ *= _scalar;
        return *this;
    }

    void vector2::normalise() {
        const float length = this->length();
        if (maths_util::approx_equal(length, 0.0f)) {
            zeroise();
            return;
        }
        x_ /= length;
        y_ /= length;
    }

    vector2 vector2::normalised() const {
        const float length = this->length();
        return maths_util::approx_equal(length, 0.0f) ? vector2::zero : vector2{x_ / length, y_ / length};
    }

    void vector2::zeroise() {
        x_ = y_ = 0.0f;
    }

    bool vector2::is_zero() const {
        return maths_util::approx_equal(0.0f, length_squared());
    }

    bool vector2::is_unit() const {
        return maths_util::approx_equal(1.0f, length_squared());
    }

    bool vector2::is_parallel(const vector2 &_vector) const {
        return !is_zero() &&
               !_vector.is_zero() &&
               maths_util::approx_equal(x_ * _vector.y_ - y_ * _vector.x_, 0.0f);
    }

    bool vector2::is_perpendicular(const vector2 &_vector) const {
        return !is_zero() &&
               !_vector.is_zero() &&
               maths_util::approx_equal(0.0f, dot(_vector));
    }

    float vector2::length() const {
        return std::sqrt(length_squared());
    }

    float vector2::length_squared() const {
        return x_ * x_ + y_ * y_;
    }

    float vector2::dot(const vector2 &_vector) const {
        return x_ * _vector.x_ + y_ * _vector.y_;
    }

    vector2 vector2::project(const vector2 &_vector) const {
        const float other_length_squared = _vector.length_squared();
        if (maths_util::approx_equal(other_length_squared, 0.0f) ||
            maths_util::approx_equal(length_squared(), 0.0f)) {
            return vector2::zero;
        }
        return (dot(_vector) * _vector) * (1.0f / other_length_squared);
    }

    float vector2::angle_between(const vector2 &_vector) const {
        return std::acos(dot(_vector) / (length() * _vector.length()));
    }

    vector2 operator*(float _scalar, const vector2 &_vector) {
        return _vector * _scalar;
    }
}