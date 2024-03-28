#include "maths/maths_util.h"
#include "maths/colour.h"

namespace mkr {
    bool colour::operator==(const colour& _rhs) const {
        return maths_util::approx_equal(r_, _rhs.r_) &&
               maths_util::approx_equal(g_, _rhs.g_) &&
               maths_util::approx_equal(b_, _rhs.b_) &&
               maths_util::approx_equal(a_, _rhs.a_);
    }

    bool colour::operator!=(const colour& _rhs) const {
        return !(*this == _rhs);
    }

    colour colour::operator+(const colour& _rhs) const {
        return colour(r_ + _rhs.r_, g_ + _rhs.g_, b_ + _rhs.b_, a_ + _rhs.a_);
    }

    colour& colour::operator+=(const colour& _rhs) {
        r_ += _rhs.r_;
        g_ += _rhs.g_;
        b_ += _rhs.b_;
        a_ += _rhs.a_;
        return *this;
    }

    colour colour::operator-(const colour& _rhs) const {
        return colour(r_ - _rhs.r_, g_ - _rhs.g_, b_ - _rhs.b_, a_ - _rhs.a_);
    }

    colour& colour::operator-=(const colour& _rhs) {
        r_ -= _rhs.r_;
        g_ -= _rhs.g_;
        b_ -= _rhs.b_;
        a_ -= _rhs.a_;
        return *this;
    }

    colour colour::operator*(const colour& _rhs) const {
        return colour(r_ * _rhs.r_, g_ * _rhs.g_, b_ * _rhs.b_, a_ * _rhs.a_);
    }

    colour& colour::operator*=(const colour& _rhs) {
        r_ *= _rhs.r_;
        g_ *= _rhs.g_;
        b_ *= _rhs.b_;
        a_ *= _rhs.a_;
        return *this;
    }

    colour colour::operator*(float _scalar) const {
        return colour(r_ * _scalar, g_ * _scalar, b_ * _scalar, a_ * _scalar);
    }

    colour& colour::operator*=(float _scalar) {
        r_ *= _scalar;
        g_ *= _scalar;
        b_ *= _scalar;
        a_ *= _scalar;
        return *this;
    }

    colour operator*(float _scalar, const colour& _colour) {
        return _colour * _scalar;
    }
}