#pragma once

#include <array>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>
#include "maths/maths_util.h"
#include "maths/vector3.h"

namespace mkr {
    template<size_t Columns, size_t Rows>
    class matrix {
    private:
        std::array<float, Columns * Rows> values_ = {};

    public:
        constexpr matrix() = default;

        constexpr matrix(std::array<float, Columns * Rows> _values) : values_{_values} {}

        static constexpr bool is_square_matrix = (Columns == Rows);

        /**
         * Returns Columns * Rows.
         * @return Columns * Rows.
         */
        static constexpr size_t size() { return Columns * Rows; }

        /**
         * Returns a zero matrix.
         * @return A zero matrix.
         */
        static constexpr matrix zero() { return matrix{}; }

        /**
         * Returns a diagonal matrix.
         * @return A diagonal matrix.
         * @warning This function is only defined for square matrices.
         */
        static constexpr matrix diagonal(float _value) requires is_square_matrix {
            matrix result;
            for (size_t i = 0; i < Columns; ++i) { result[i][i] = _value; }
            return result;
        }

        /**
         * Returns an identity matrix.
         * @return An identity matrix.
         * @warning This function is only defined for square matrices.
         */
        static constexpr matrix identity() requires is_square_matrix { return diagonal(1.0f); }

        /**
         * @brief Returns the transpose of a matrix.
         * @return The transpose of a matrix.
         */
        matrix<Rows, Columns> transposed() const {
            matrix<Rows, Columns> mat;
            for (size_t i = 0; i < Columns; ++i) {
                for (size_t j = 0; j < Rows; ++j) {
                    mat[j][i] = (*this)[i][j];
                }
            }
            return mat;
        }

        inline const float* operator[](size_t _column) const {
            return &values_[_column * Rows];
        }

        inline float* operator[](size_t _column) {
            return &values_[_column * Rows];
        }

        bool operator==(const matrix& _rhs) const {
            bool equal = true;
            for (size_t i = 0; i < size(); i++) {
                equal &= maths_util::approx_equal(values_[i], _rhs.values_[i]);
            }
            return equal;
        }

        bool operator!=(const matrix& _rhs) const {
            return !(*this == _rhs);
        }

        matrix operator+(const matrix& _rhs) const {
            matrix result;
            for (size_t i = 0; i < size(); ++i) {
                result.values_[i] = values_[i] + _rhs.values_[i];
            }
            return result;
        }

        matrix& operator+=(const matrix& _rhs) {
            *this = (*this) + _rhs;
            return *this;
        }

        matrix operator-(const matrix& _rhs) const {
            matrix result;
            for (size_t i = 0; i < size(); ++i) {
                result.values_[i] = values_[i] - _rhs.values_[i];
            }
            return result;
        }

        matrix& operator-=(const matrix& _rhs) {
            *this = (*this) - _rhs;
            return *this;
        }

        template<size_t RHSColumns>
        matrix<RHSColumns, Rows> operator*(const matrix<RHSColumns, Columns>& _rhs) const {
            matrix<RHSColumns, Rows> result;
            for (size_t i = 0; i < RHSColumns; ++i) {
                for (size_t j = 0; j < Rows; ++j) {
                    for (size_t k = 0; k < Columns; ++k) {
                        result[i][j] += (*this)[k][j] * _rhs[i][k];
                    }
                }
            }
            return result;
        }

        matrix& operator*=(const matrix& _rhs) requires is_square_matrix {
            *this = (*this) * _rhs;
            return *this;
        }

        matrix operator*(float _scalar) const {
            matrix result;
            for (size_t i = 0; i < size(); ++i) {
                result.values_[i] = values_[i] * _scalar;
            }
            return result;
        }

        matrix& operator*=(float _scalar) {
            *this = (*this) * _scalar;
            return *this;
        }

        vector3 operator*(const vector3 _rhs) const requires (Columns == 4 && Rows == 4) {
            auto point = (*this) * matrix<1, 4>{{_rhs.x_, _rhs.y_, _rhs.z_, 1.0f}};
            return vector3{point[0][0] / point[0][3], point[0][1] / point[0][3], point[0][2] / point[0][3]};
        }

        friend matrix operator*(float _scalar, const matrix& _matrix) {
            return _matrix * _scalar;
        }

        /**
         * @brief Returns a string representation of this matrix.
         * 
         * @param _precision The number of decimals to show for each element
         * @param _padding The padding for each element so that it is printed aligned
         * @return A string representation of this matrix.
         */
        [[nodiscard]] std::string to_string(const int _precision = 4) const {
            std::ostringstream out;
            out.precision(_precision);
            out << std::fixed;
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Columns; ++j) {
                    out << (*this)[j][i];
                    if (j < Columns - 1) { out << ", "; }
                }
                out << '\n';
            }
            return out.str();
        }

        friend std::ostream& operator<<(std::ostream& _stream, const matrix<Columns, Rows>& _matrix) {
            return _stream << _matrix.to_string();
        }
    };

    // Do not allow matrices with 0 rows or columns.
    template<size_t Rows>
    class matrix<0, Rows>;

    template<size_t Columns>
    class matrix<Columns, 0>;

    typedef matrix<1, 1> matrix1x1;
    typedef matrix<1, 2> matrix1x2;
    typedef matrix<1, 3> matrix1x3;
    typedef matrix<1, 4> matrix1x4;
    typedef matrix<1, 5> matrix1x5;

    typedef matrix<2, 1> matrix2x1;
    typedef matrix<2, 2> matrix2x2;
    typedef matrix<2, 3> matrix2x3;
    typedef matrix<2, 4> matrix2x4;
    typedef matrix<2, 5> matrix2x5;

    typedef matrix<3, 1> matrix3x1;
    typedef matrix<3, 2> matrix3x2;
    typedef matrix<3, 3> matrix3x3;
    typedef matrix<3, 4> matrix3x4;
    typedef matrix<3, 5> matrix3x5;

    typedef matrix<4, 1> matrix4x1;
    typedef matrix<4, 2> matrix4x2;
    typedef matrix<4, 3> matrix4x3;
    typedef matrix<4, 4> matrix4x4;
    typedef matrix<4, 5> matrix4x5;

    typedef matrix<5, 1> matrix5x1;
    typedef matrix<5, 2> matrix5x2;
    typedef matrix<5, 3> matrix5x3;
    typedef matrix<5, 4> matrix5x4;
    typedef matrix<5, 5> matrix5x5;
}