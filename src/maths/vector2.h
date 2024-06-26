#pragma once

#include <sstream>
#include "maths/maths_util.h"

namespace mkr {
    /**
     * A 2D vector with x and y components.
     */
    class vector2 {
    public:
        static vector2 zero() { return vector2{0.0f, 0.0f}; }
        static vector2 up() { return vector2{0.0f, 1.0f}; }
        static vector2 down() { return vector2{0.0f, -1.0f}; }
        static vector2 left() { return vector2{1.0f, 0.0f}; }
        static vector2 right() { return vector2{-1.0f, 0.0f}; }
        static vector2 x_axis() { return vector2{1.0f, 0.0f}; }
        static vector2 y_axis() { return vector2{0.0f, 1.0f}; }

        /// The x component.
        float x_;
        /// The y component.
        float y_;

        /**
         * Constructs the vector.
         * @param _x The x component.
         * @param _y The y component.
         */
        vector2(float _x = 0.0f, float _y = 0.0f);

        vector2 operator-() const;

        bool operator==(const vector2& _rhs) const;

        bool operator!=(const vector2& _rhs) const;

        vector2 operator+(const vector2& _rhs) const;

        vector2& operator+=(const vector2& _rhs);

        vector2 operator-(const vector2& _rhs) const;

        vector2& operator-=(const vector2& _rhs);

        vector2 operator*(float _rhs) const;

        vector2& operator*=(float _rhs);

        vector2 operator*(const vector2& _rhs) const;

        vector2& operator*=(const vector2& _rhs);

        /**
         * Normalise this vector.
         */
        void normalise();

        /**
         * Returns a normalised copy of this vector.
         * @return A normalised copy of this vector.
         */
        [[nodiscard]] vector2 normalised() const;

        /**
         * Checks if this vector is a zero vector.
         * @return Returns true if the vector is a zero vector, else return false.
         */
        [[nodiscard]] bool is_zero() const;

        /**
         * Checks if this vector is a unit vector.
         */
        [[nodiscard]] bool is_unit() const;

        /**
         * Checks if 2 vectors are parallel.
         * @param _vector The vector to compare to.
         * @return Returns true if the 2 vectors are parallel, else returns false.
         */
        [[nodiscard]] bool is_parallel(const vector2& _vector) const;

        /**
         * Checks if 2 vectors are perpendicular.
         * @param _vector The vector to compare to.
         * @return Returns true if the 2 vectors are perpendicular, else returns false.
         */
        [[nodiscard]] bool is_perpendicular(const vector2& _vector) const;

        /**
         * Returns the length of this vector.
         * @return The length of this vector.
         */
        [[nodiscard]] float length() const;

        /**
         * Returns the squared length of this vector.
         * @return The squared length of this vector.
         */
        [[nodiscard]] float length_squared() const;

        /**
         * Returns the dot product of 2 vectors.
         * @param _vector　The vector to dot with.
         * @return The dot product of 2 vectors.
         */
        [[nodiscard]] float dot(const vector2& _vector) const;

        /**
         * Returns the projection of this vector onto another vector.
         * @param _vector The vector to project this vector onto.
         * @return The projection of this vector onto another vector.
         */
        [[nodiscard]] vector2 project(const vector2& _vector) const;

        /**
         * Returns the angle between 2 vectors.
         * @param _vector The other vector to find the angle with.
         * @return The angle between 2 vectors.
         */
        [[nodiscard]] float angle_between(const vector2& _vector) const;

        friend vector2 operator*(float _lhs, const vector2& _vector);

        [[nodiscard]] std::string to_string(const int _precision = 4) const {
            std::ostringstream out;
            out.precision(_precision);
            out << std::fixed;
            out << x_ << ", " << y_;
            return out.str();
        }

        friend std::ostream& operator<<(std::ostream& _stream, const vector2& _vector) {
            return _stream << _vector.to_string();
        }
    };
}