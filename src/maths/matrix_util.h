#pragma once

#include <optional>
#include "maths/matrix.h"
#include "maths/maths_util.h"
#include "maths/vector3.h"

namespace mkr {
    class matrix_util {
    public:
        matrix_util() = delete;

        /**
         * @brief Returns the transpose of a matrix.
         * @return The transpose of a matrix.
         */
        template<size_t Columns, size_t Rows>
        static matrix<Rows, Columns> transpose_matrix(const matrix<Columns, Rows> &_matrix) {
            matrix<Rows, Columns> mat;
            for (size_t i = 0; i < Columns; ++i) {
                for (size_t j = 0; j < Rows; ++j) {
                    mat[j][i] = _matrix[i][j];
                }
            }
            return mat;
        }

        /**
         * @brief Get the determinant of a 1 by 1 matrix
         *
         * @param _matrix the matrix that will be used
         * @return float the determinant of the matrix
         */
        static inline float determinant(const matrix1x1 &_matrix) {
            return _matrix[0][0];
        };

        /**
         * @brief Get the determinant of a 2 by 2 matrix
         *
         * @param _matrix the matrix that will be used
         * @return float the determinant of the matrix
         */
        static inline float determinant(const matrix2x2 &_matrix) {
            return (_matrix[0][0] * _matrix[1][1]) -
                   (_matrix[1][0] * _matrix[0][1]);
        };

        /**
         * @brief Get the determinant of a 3 by 3 matrix
         *
         * @param _matrix the matrix that will be used
         * @return float the determinant of the matrix
         */
        static inline float determinant(const matrix3x3 &_matrix) {
            return (_matrix[0][0] * _matrix[1][1] * _matrix[2][2]) +
                   (_matrix[1][0] * _matrix[2][1] * _matrix[0][2]) +
                   (_matrix[2][0] * _matrix[0][1] * _matrix[1][2]) -
                   (_matrix[0][2] * _matrix[1][1] * _matrix[2][0]) -
                   (_matrix[1][2] * _matrix[2][1] * _matrix[0][0]) -
                   (_matrix[2][2] * _matrix[0][1] * _matrix[1][0]);
        };

        /**
         * @brief Get the determinant of a 4 by 4 matrix
         *
         * @param _matrix the matrix that will be used
         * @return float the determinant of the matrix
         */
        static inline float determinant(const matrix4x4 &_matrix) {
            const float m00 = _matrix[0][5] * _matrix[0][10] * _matrix[0][15] -
                              _matrix[0][5] * _matrix[0][11] * _matrix[0][14] -
                              _matrix[0][9] * _matrix[0][6] * _matrix[0][15] +
                              _matrix[0][9] * _matrix[0][7] * _matrix[0][14] +
                              _matrix[0][13] * _matrix[0][6] * _matrix[0][11] -
                              _matrix[0][13] * _matrix[0][7] * _matrix[0][10];
            const float m01 = -_matrix[0][4] * _matrix[0][10] * _matrix[0][15] +
                              _matrix[0][4] * _matrix[0][11] * _matrix[0][14] +
                              _matrix[0][8] * _matrix[0][6] * _matrix[0][15] -
                              _matrix[0][8] * _matrix[0][7] * _matrix[0][14] -
                              _matrix[0][12] * _matrix[0][6] * _matrix[0][11] +
                              _matrix[0][12] * _matrix[0][7] * _matrix[0][10];
            const float m02 = _matrix[0][4] * _matrix[0][9] * _matrix[0][15] -
                              _matrix[0][4] * _matrix[0][11] * _matrix[0][13] -
                              _matrix[0][8] * _matrix[0][5] * _matrix[0][15] +
                              _matrix[0][8] * _matrix[0][7] * _matrix[0][13] +
                              _matrix[0][12] * _matrix[0][5] * _matrix[0][11] -
                              _matrix[0][12] * _matrix[0][7] * _matrix[0][9];
            const float m03 = -_matrix[0][4] * _matrix[0][9] * _matrix[0][14] +
                              _matrix[0][4] * _matrix[0][10] * _matrix[0][13] +
                              _matrix[0][8] * _matrix[0][5] * _matrix[0][14] -
                              _matrix[0][8] * _matrix[0][6] * _matrix[0][13] -
                              _matrix[0][12] * _matrix[0][5] * _matrix[0][10] +
                              _matrix[0][12] * _matrix[0][6] * _matrix[0][9];

            return _matrix[0][0] * m00 +
                   _matrix[0][1] * m01 +
                   _matrix[0][2] * m02 +
                   _matrix[0][3] * m03;
        };

        /**
         * @brief Get the determinant of given square matrix
         *
         * @tparam Columns the columns/rows of the matrix
         * @param _matrix the matrix that will be used
         * @return float the determinant of the matrix
         */
        template<size_t Columns>
        static float determinant(const matrix<Columns, Columns> &_matrix) {
            float det = 0.0f;
            for (size_t col = 0; col < Columns; ++col) {
                int sign = 1 + ((int) col & 1) * -2; // Determine sign-ness of that iteration.
                float cofactor = _matrix[col][0];
                float minor_det = determinant(minor_matrix(_matrix, col, 0));
                det += (float) sign * cofactor * minor_det;
            }
            return det;
        }

        /**
         * @brief Get the minor square matrix given a matrix of a bigger size
         *
         * @tparam Columns the size of the original matrix
         * @param _matrix the original matrix
         * @param _cofactor_col the column that will not be used the minor matrix
         * @param _cofactor_row the row that will not be used in the minor matrix
         * @return matrix<Columns - 1, Columns - 1> the minor square matrix of size - 1
         */
        template<size_t Columns>
        static matrix<Columns - 1, Columns - 1>
        minor_matrix(const matrix<Columns, Columns> &_matrix, size_t _cofactor_col, size_t _cofactor_row) {
            matrix<Columns - 1, Columns - 1> mat;
            for (size_t major_col = 0, minor_col = 0; major_col < Columns; ++major_col) {
                if (major_col == _cofactor_col) continue;
                for (size_t major_row = 0, minor_row = 0; major_row < Columns; ++major_row) {
                    if (major_row == _cofactor_row) continue;
                    mat[minor_col][minor_row] = _matrix[major_col][major_row];
                    ++minor_row;
                }
                ++minor_col;
            }
            return mat;
        }

        /**
         * @brief Get the cofactor matrix of a 1 by 1 matrix
         *
         * @param _matrix the original matrix
         * @return matrix<1, 1> the cofactor matrix generated from the original matrix
         */
        static inline matrix<1, 1> cofactor_matrix(const matrix<1, 1> &_matrix) {
            return _matrix;
        }

        /**
         * @brief Get the cofactor matrix of a given matrix
         *
         * @tparam Columns the size of the matrix
         * @param _matrix the original matrix
         * @return matrix<Columns, Columns> the cofactor matrix generated from the original matrix
         */
        template<size_t Columns>
        static matrix<Columns, Columns> cofactor_matrix(const matrix<Columns, Columns> &_matrix) {
            matrix<Columns, Columns> mat;
            for (size_t col = 0; col < Columns; ++col) {
                for (size_t row = 0; row < Columns; ++row) {
                    int sign = 1 + ((int) (col + row) & 1) * -2; // Determine sign-ness of that iteration.
                    float det = determinant(minor_matrix(_matrix, col, row));
                    mat[col][row] = (float) sign * det;
                }
            }
            return mat;
        }

        /**
         * @brief Get the adjugate matrix of a given matrix
         *
         * @tparam Size the size of the original matrix
         * @param _matrix the original matrix
         * @return matrix<Size, Size> the adjugate matrix gernerated from the original matrix
         */
        template<size_t Size>
        static matrix<Size, Size> adjugate_matrix(const matrix<Size, Size> &_matrix) {
            return transpose_matrix(cofactor_matrix(_matrix));
        }

        /**
         * @brief boolean to check if an inverse matrix exist
         *
         * @param _matrix the matrix to check
         * @return bool true if the matrix has an inverse, else false
         */
        template<size_t Size>
        static bool invertible(const matrix<Size, Size> &_matrix) {
            return !maths_util::approx_equal(0.0f, determinant(_matrix));
        }

        /**
         * @brief Get the inverse matrix of 1 by 1 matrix
         *
         * @param _matrix the matrix to find the inverse of
         * @return std::optional<matrix1x1> the inverse matrix if it exists, else std::nullopt
         */
        static std::optional<matrix1x1> inverse_matrix(const matrix1x1 &_matrix) {
            const float det = determinant(_matrix);
            if (maths_util::approx_equal(0.0f, det)) return std::nullopt;

            matrix1x1 mat;
            mat[0][0] = 1.0f / _matrix[0][0];
            return mat;
        }

        /**
         * @brief Get the inverse matrix of 2 by 2 matrix
         *
         * @param _matrix the matrix to find the inverse of
         * @return std::optional<matrix2x2> the inverse matrix if it exists, else std::nullopt
         */
        static std::optional<matrix2x2> inverse_matrix(const matrix2x2 &_matrix) {
            const float det = determinant(_matrix);
            if (maths_util::approx_equal(0.0f, det)) return std::nullopt;

            matrix2x2 mat = _matrix;
            mat[0][1] = -mat[0][1];
            mat[1][0] = -mat[1][0];
            std::swap(mat[0][0], mat[1][1]);
            return 1.0f / det * mat;
        }

        /**
         * @brief Get the inverse matrix of 3 by 3 matrix
         *
         * @param _matrix the matrix to find the inverse of
         * @return std::optional<matrix3x3> the inverse matrix if it exists, else std::nullopt
         */
        static std::optional<matrix3x3> inverse_matrix(const matrix3x3 &_matrix) {
            matrix3x3 mat;
            mat[0][0] = _matrix[1][1] * _matrix[2][2] - _matrix[1][2] * _matrix[2][1];
            mat[1][0] = -(_matrix[1][0] * _matrix[2][2] - _matrix[1][2] * _matrix[2][0]);
            mat[2][0] = _matrix[1][0] * _matrix[2][1] - _matrix[1][1] * _matrix[2][0];

            const float det = mat[0][0] * _matrix[0][0] +
                              mat[1][0] * _matrix[0][1] +
                              mat[2][0] * _matrix[0][2];
            if (maths_util::approx_equal(0.0f, det)) return std::nullopt;

            mat[0][1] = -(_matrix[0][1] * _matrix[2][2] - _matrix[0][2] * _matrix[2][1]);
            mat[1][1] = _matrix[0][0] * _matrix[2][2] - _matrix[0][2] * _matrix[2][0];
            mat[2][1] = -(_matrix[0][0] * _matrix[2][1] - _matrix[0][1] * _matrix[2][0]);
            mat[0][2] = _matrix[0][1] * _matrix[1][2] - _matrix[0][2] * _matrix[1][1];
            mat[1][2] = -(_matrix[0][0] * _matrix[1][0] - _matrix[0][2] * _matrix[1][2]);
            mat[2][2] = _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
            return 1.0f / det * mat;
        }

        /**
         * @brief Get the inverse matrix of 4 by 4 matrix
         *
         * @param _matrix the matrix to find the inverse of
         * @return std::optional<matrix4x4> the inverse matrix if it exists, else std::nullopt
         */
        static std::optional<matrix4x4> inverse_matrix(const matrix4x4 &_matrix) {
            matrix4x4 mat;
            mat[0][0] = _matrix[0][5] * _matrix[0][10] * _matrix[0][15] -
                        _matrix[0][5] * _matrix[0][11] * _matrix[0][14] -
                        _matrix[0][9] * _matrix[0][6] * _matrix[0][15] +
                        _matrix[0][9] * _matrix[0][7] * _matrix[0][14] +
                        _matrix[0][13] * _matrix[0][6] * _matrix[0][11] -
                        _matrix[0][13] * _matrix[0][7] * _matrix[0][10];
            mat[0][4] = -_matrix[0][4] * _matrix[0][10] * _matrix[0][15] +
                        _matrix[0][4] * _matrix[0][11] * _matrix[0][14] +
                        _matrix[0][8] * _matrix[0][6] * _matrix[0][15] -
                        _matrix[0][8] * _matrix[0][7] * _matrix[0][14] -
                        _matrix[0][12] * _matrix[0][6] * _matrix[0][11] +
                        _matrix[0][12] * _matrix[0][7] * _matrix[0][10];
            mat[0][8] = _matrix[0][4] * _matrix[0][9] * _matrix[0][15] -
                        _matrix[0][4] * _matrix[0][11] * _matrix[0][13] -
                        _matrix[0][8] * _matrix[0][5] * _matrix[0][15] +
                        _matrix[0][8] * _matrix[0][7] * _matrix[0][13] +
                        _matrix[0][12] * _matrix[0][5] * _matrix[0][11] -
                        _matrix[0][12] * _matrix[0][7] * _matrix[0][9];
            mat[0][12] = -_matrix[0][4] * _matrix[0][9] * _matrix[0][14] +
                         _matrix[0][4] * _matrix[0][10] * _matrix[0][13] +
                         _matrix[0][8] * _matrix[0][5] * _matrix[0][14] -
                         _matrix[0][8] * _matrix[0][6] * _matrix[0][13] -
                         _matrix[0][12] * _matrix[0][5] * _matrix[0][10] +
                         _matrix[0][12] * _matrix[0][6] * _matrix[0][9];

            const float det = _matrix[0][0] * mat[0][0] +
                              _matrix[0][1] * mat[0][4] +
                              _matrix[0][2] * mat[0][8] +
                              _matrix[0][3] * mat[0][12];
            if (maths_util::approx_equal(0.0f, det)) return std::nullopt;

            mat[0][1] = -_matrix[0][1] * _matrix[0][10] * _matrix[0][15] +
                        _matrix[0][1] * _matrix[0][11] * _matrix[0][14] +
                        _matrix[0][9] * _matrix[0][2] * _matrix[0][15] -
                        _matrix[0][9] * _matrix[0][3] * _matrix[0][14] -
                        _matrix[0][13] * _matrix[0][2] * _matrix[0][11] +
                        _matrix[0][13] * _matrix[0][3] * _matrix[0][10];
            mat[0][5] = _matrix[0][0] * _matrix[0][10] * _matrix[0][15] -
                        _matrix[0][0] * _matrix[0][11] * _matrix[0][14] -
                        _matrix[0][8] * _matrix[0][2] * _matrix[0][15] +
                        _matrix[0][8] * _matrix[0][3] * _matrix[0][14] +
                        _matrix[0][12] * _matrix[0][2] * _matrix[0][11] -
                        _matrix[0][12] * _matrix[0][3] * _matrix[0][10];
            mat[0][9] = -_matrix[0][0] * _matrix[0][9] * _matrix[0][15] +
                        _matrix[0][0] * _matrix[0][11] * _matrix[0][13] +
                        _matrix[0][8] * _matrix[0][1] * _matrix[0][15] -
                        _matrix[0][8] * _matrix[0][3] * _matrix[0][13] -
                        _matrix[0][12] * _matrix[0][1] * _matrix[0][11] +
                        _matrix[0][12] * _matrix[0][3] * _matrix[0][9];
            mat[0][13] = _matrix[0][0] * _matrix[0][9] * _matrix[0][14] -
                         _matrix[0][0] * _matrix[0][10] * _matrix[0][13] -
                         _matrix[0][8] * _matrix[0][1] * _matrix[0][14] +
                         _matrix[0][8] * _matrix[0][2] * _matrix[0][13] +
                         _matrix[0][12] * _matrix[0][1] * _matrix[0][10] -
                         _matrix[0][12] * _matrix[0][2] * _matrix[0][9];
            mat[0][2] = _matrix[0][1] * _matrix[0][6] * _matrix[0][15] -
                        _matrix[0][1] * _matrix[0][7] * _matrix[0][14] -
                        _matrix[0][5] * _matrix[0][2] * _matrix[0][15] +
                        _matrix[0][5] * _matrix[0][3] * _matrix[0][14] +
                        _matrix[0][13] * _matrix[0][2] * _matrix[0][7] -
                        _matrix[0][13] * _matrix[0][3] * _matrix[0][6];
            mat[0][6] = -_matrix[0][0] * _matrix[0][6] * _matrix[0][15] +
                        _matrix[0][0] * _matrix[0][7] * _matrix[0][14] +
                        _matrix[0][4] * _matrix[0][2] * _matrix[0][15] -
                        _matrix[0][4] * _matrix[0][3] * _matrix[0][14] -
                        _matrix[0][12] * _matrix[0][2] * _matrix[0][7] +
                        _matrix[0][12] * _matrix[0][3] * _matrix[0][6];
            mat[0][10] = _matrix[0][0] * _matrix[0][5] * _matrix[0][15] -
                         _matrix[0][0] * _matrix[0][7] * _matrix[0][13] -
                         _matrix[0][4] * _matrix[0][1] * _matrix[0][15] +
                         _matrix[0][4] * _matrix[0][3] * _matrix[0][13] +
                         _matrix[0][12] * _matrix[0][1] * _matrix[0][7] -
                         _matrix[0][12] * _matrix[0][3] * _matrix[0][5];
            mat[0][14] = -_matrix[0][0] * _matrix[0][5] * _matrix[0][14] +
                         _matrix[0][0] * _matrix[0][6] * _matrix[0][13] +
                         _matrix[0][4] * _matrix[0][1] * _matrix[0][14] -
                         _matrix[0][4] * _matrix[0][2] * _matrix[0][13] -
                         _matrix[0][12] * _matrix[0][1] * _matrix[0][6] +
                         _matrix[0][12] * _matrix[0][2] * _matrix[0][5];
            mat[0][3] = -_matrix[0][1] * _matrix[0][6] * _matrix[0][11] +
                        _matrix[0][1] * _matrix[0][7] * _matrix[0][10] +
                        _matrix[0][5] * _matrix[0][2] * _matrix[0][11] -
                        _matrix[0][5] * _matrix[0][3] * _matrix[0][10] -
                        _matrix[0][9] * _matrix[0][2] * _matrix[0][7] +
                        _matrix[0][9] * _matrix[0][3] * _matrix[0][6];
            mat[0][7] = _matrix[0][0] * _matrix[0][6] * _matrix[0][11] -
                        _matrix[0][0] * _matrix[0][7] * _matrix[0][10] -
                        _matrix[0][4] * _matrix[0][2] * _matrix[0][11] +
                        _matrix[0][4] * _matrix[0][3] * _matrix[0][10] +
                        _matrix[0][8] * _matrix[0][2] * _matrix[0][7] -
                        _matrix[0][8] * _matrix[0][3] * _matrix[0][6];
            mat[0][11] = -_matrix[0][0] * _matrix[0][5] * _matrix[0][11] +
                         _matrix[0][0] * _matrix[0][7] * _matrix[0][9] +
                         _matrix[0][4] * _matrix[0][1] * _matrix[0][11] -
                         _matrix[0][4] * _matrix[0][3] * _matrix[0][9] -
                         _matrix[0][8] * _matrix[0][1] * _matrix[0][7] +
                         _matrix[0][8] * _matrix[0][3] * _matrix[0][5];
            mat[0][15] = _matrix[0][0] * _matrix[0][5] * _matrix[0][10] -
                         _matrix[0][0] * _matrix[0][6] * _matrix[0][9] -
                         _matrix[0][4] * _matrix[0][1] * _matrix[0][10] +
                         _matrix[0][4] * _matrix[0][2] * _matrix[0][9] +
                         _matrix[0][8] * _matrix[0][1] * _matrix[0][6] -
                         _matrix[0][8] * _matrix[0][2] * _matrix[0][5];
            return 1.0f / det * mat;
        }

        /**
         * @brief Get the inverse matrix of square matrix
         *
         * @param _matrix the matrix to find the inverse of
         * @return std::optional<matrix<Columns, Columns>> the inverse matrix if it exists, else std::nullopt
         */
        template<size_t Columns>
        static std::optional<matrix<Columns, Columns>> inverse_matrix(const matrix<Columns, Columns> &_matrix) {
            const float det = determinant(_matrix);
            if (maths_util::approx_equal(0.0f, det)) return std::nullopt;
            return 1.0f / det * adjugate_matrix(_matrix);
        }

        /**
         * @brief Generate the homogeneous translation matrix given a translation vector
         * 
         * @param _translation the translation vectors data that will be used
         * @return matrix4x4 the homogeneous translation matrix
         */
        static matrix4x4 translation_matrix(const vector3 &_translation) {
            /**
             * Translation Matrix
             * |   1   0   0   x   |
             * |   0   1   1   y   |
             * |   0   0   1   z   |
             * |   0   0   0   1   |
             */
            matrix4x4 mat = matrix4x4::identity();
            mat[3][0] = _translation.x_;
            mat[3][1] = _translation.y_;
            mat[3][2] = _translation.z_;

            return mat;
        }

        /**
         * @brief Generate the homogeneous rotation matrix about the X-axis given angle in radians.
         * 
         * @param _angle angle to rotate about in radians
         * @return matrix4x4 homogeneous rotation matrix about the X-axis given angle in radians.
         */
        static matrix4x4 rotation_matrix_x(float _angle) {
            /**
             * Rotation Matrix on X-axis
             * |   1   0   0   0   |
             * |   0  cos -sin 0   |
             * |   0  sin  cos 0   |
             * |   0   0   0   1   |
             */
            matrix4x4 mat = matrix4x4::identity();
            mat[1][1] = mat[2][2] = std::cos(_angle);
            mat[1][2] = std::sin(_angle);
            mat[2][1] = -mat[1][2];

            return mat;
        }

        /**
         * @brief Generate the homogeneous rotation matrix about the Y-axis given angle in radians.
         * 
         * @param _angle angle to rotate about in radians
         * @return matrix4x4 homogeneous rotation matrix about the Y-axis given angle in radians.
         */
        static matrix4x4 rotation_matrix_y(float _angle) {
            /**
             * Rotation Matrix on Y-axis
             * |  cos  0  sin  0   |
             * |   0   1   0   0   |
             * | -sin  0  cos  0   |
             * |   0   0   0   1   |
             */
            matrix4x4 mat = matrix4x4::identity();
            mat[0][0] = mat[2][2] = std::cos(_angle);
            mat[2][0] = std::sin(_angle);
            mat[0][2] = -mat[2][0];

            return mat;
        }

        /**
         * @brief Generate the homogeneous rotation matrix about the Z-axis given angle in radians.
         * 
         * @param _angle angle to rotate about in radians
         * @return matrix4x4 homogeneous rotation matrix about the Z-axis given angle in radians.
         */
        static matrix4x4 rotation_matrix_z(float _angle) {
            /**
             * Rotation Matrix on Y-axis
             * |  cos -sin 0   0   |
             * |  sin  cos 0   0   |
             * |   0   0   1   0   |
             * |   0   0   0   1   |
             */
            matrix4x4 mat = matrix4x4::identity();
            mat[0][0] = mat[1][1] = std::cos(_angle);
            mat[0][1] = std::sin(_angle);
            mat[1][0] = -mat[0][1];

            return mat;
        }

        /**
         * @brief Generate the homogeneous rotation matrix about all 3 XYZ-axis 
         * given a vector of euler angles rotations.
         * 
         * @param _euler_angles the euler angles to the XYZ axis
         * @return matrix4x4 homogeneous rotation matrix about the 3 XYZ-axis.
         */
        static matrix4x4 rotation_matrix(const vector3 &_euler_angles) {
            return rotation_matrix_x(_euler_angles.x_) *
                   rotation_matrix_y(_euler_angles.y_) *
                   rotation_matrix_z(_euler_angles.z_);
        }

        /**
         * @brief Generate the homogeneous scale matrix given scale vector
         * 
         * @param _scale the vector representating the scale of the object in XYZ
         * @return matrix4x4 
         */
        static matrix4x4 scale_matrix(const vector3 &_scale) {
            /**
             * Scale Matrix
             * |   x   0   0   0   |
             * |   0   y   1   0   |
             * |   0   0   z   0   |
             * |   0   0   0   1   |
             */
            matrix4x4 mat;
            mat[0][0] = _scale.x_;
            mat[1][1] = _scale.y_;
            mat[2][2] = _scale.z_;
            mat[3][3] = 1.0f;

            return mat;
        }

        /**
         * @brief Get the model matrix given translation, rotation and scale.
         * 
         * @param _translation translation of the object in 3D space
         * @param _euler_angles rotation of the object in 3D space
         * @param _scale the scale of the object in 3D space
         * @return matrix4x4 the Homogeneous matrix generated
         */
        static matrix4x4 model_matrix(const vector3 &_translation,
                                      const vector3 &_euler_angles,
                                      const vector3 &_scale) {
            return translation_matrix(_translation) *
                   rotation_matrix(_euler_angles) *
                   scale_matrix(_scale);
        }

        /**
         * @brief Get the view matrix given the cameras orientation and position
         * 
         * @param _forward the forward vector
         * @param _up the up vector
         * @param _position the position in 3D space
         * @return matrix4x4 the view matrix in 3D space
         */
        static matrix4x4 view_matrix(const vector3 &_forward, const vector3 &_up, const vector3 &_position) {
            vector3 right = _forward.cross(_up);
            /**
             * translation * orientation matrix  = homogeneous matrix
             * [ 1 0 0 T1 ][ R11  R12  R13  0 ]    [ R11 R12 R13 T1 ]
             * [ 0 1 0 T2 ][ R21  R22  R23  0 ]  = [ R21 R22 R23 T2 ]
             * [ 0 0 1 T3 ][ R31  R32  R33  0 ]    [ R31 R32 R33 T3 ]
             * [ 0 0 0 1  ][  0    0    0   1 ]    [  0   0   0  1  ]
             *
             * [  Right   trans ]
             * [  Up      trans ]
             * [ -Forward trans ]
             * [  0   0   0   1  ]
             */
            matrix4x4 mat;
            mat[3][3] = 1.0f;

            // Column 0
            mat[0][0] = right.x_;
            mat[0][1] = _up.x_;
            mat[0][2] = -_forward.x_;

            // Column 1
            mat[1][0] = right.y_;
            mat[1][1] = _up.y_;
            mat[1][2] = -_forward.y_;

            // Column 2
            mat[2][0] = right.z_;
            mat[2][1] = _up.z_;
            mat[2][2] = -_forward.z_;

            // Column 3
            mat[3][0] = _position.x_;
            mat[3][1] = _position.y_;
            mat[3][2] = _position.z_;

            return mat;
        }

        static matrix4x4 perspective_matrix(float _aspect_ratio, float _fov, float _near_plane, float _far_plane) {
            matrix4x4 mat = matrix4x4::identity();

            /**
             * AR  = Aspect Ratio
             * FOV = Field of View
             * N   = Near Plane
             * F   = Far Plane
             *
             * Perspective Matrix:
             * | 1/(tan(FOV/2) * AR)         0                0              0       |
             * |          0            1/tan(FOV/2)           0              0       |
             * |          0                  0          (N+F)/(N-F)   (2*F*N)/(N-F)  |
             * |          0                  0               -1              0       |
             * */

            float b = 1.0f / std::tan(0.5f * _fov);
            float a = b / _aspect_ratio;
            float c = (_near_plane + _far_plane) / (_near_plane - _far_plane);
            float d = (2.0f * _near_plane * _far_plane) / (_near_plane - _far_plane);

            mat[0][0] = a;
            mat[1][1] = b;
            mat[2][2] = c;
            mat[2][3] = -1.0f;
            mat[3][2] = d;

            return mat;
        }

        static matrix4x4 orthographic_matrix(float _left, float _right, float _top, float _bottom,
                                             float _near_plane, float _far_plane) {
            matrix4x4 mat = matrix4x4::identity();

            mat[0][0] = 2.0f / (_right - _left);
            mat[1][1] = 2.0f / (_top - _bottom);
            mat[2][2] = 2.0f / (_near_plane - _far_plane);

            mat[3][0] = (_right + _left) / (_left - _right);
            mat[3][1] = (_top + _bottom) / (_bottom - _top);
            mat[3][2] = (_far_plane + _near_plane) / (_near_plane - _far_plane);

            mat[3][3] = 1.0f;

            return mat;
        }

        matrix4x4 orthographic_matrix(float _aspect_ratio, float _ortho_size, float _near_plane, float _far_plane) {
            // [https://en.wikipedia.org/wiki/Orthographic_projection]
            // Most tutorial take the absolute position of the viewing box in the world as the input.
            // So there is a need to translate the box back to the origin.
            // We do not need to do so as we do our calculation assuming that we are already at the origin.

            matrix4x4 mat = matrix4x4::identity();

            float top = _ortho_size * 0.5f;
            float bottom = -top;
            float right = top * _aspect_ratio;
            float left = bottom * _aspect_ratio;

            mat[0][0] = 2.0f / (right - left);
            mat[1][1] = 2.0f / (top - bottom);
            mat[2][2] = 2.0f / (_far_plane - _near_plane);

            return mat;
        }
    };
}