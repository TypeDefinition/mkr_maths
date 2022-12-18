#pragma once

#include <optional>
#include "maths/vector3.h"

namespace mkr {
    /**
     * A line represented by the formula p + λd.
     */
    class line {
    public:
        static const line x_axis;
        static const line y_axis;
        static const line z_axis;

        /**
         * A point on the line.
         */
        vector3 point_;
        /**
         * The direction of the line.
         */
        vector3 direction_;

        /**
         * Constructs the line.
         * @param _point A point on the line.
         * @param _direction The direction of the line.
         */
        line(const vector3& _point = vector3::zero, const vector3& _direction = vector3::x_axis);

        bool operator==(const line& _line) const;

        /**
         * Checks if 2 lines are parallel.
         * @param _line The line to compare to.
         * @return Returns true if the 2 line are parallel, else returns false.
         */
        [[nodiscard]] bool is_parallel(const line& _line) const;

        /**
         * Checks if this line is parallel to a vector.
         * @param _vector The vector to compare to.
         * @return Returns true if this line is parallel to the vector, else returns false.
         */
        [[nodiscard]] bool is_parallel(const vector3& _vector) const;

        /**
         * Checks if 2 lines are perpendicular.
         * @param _line The line to compare to.
         * @return Returns true if the 2 line are perpendicular, else returns false.
         */
        [[nodiscard]] bool is_perpendicular(const line& _line) const;

        /**
         * Checks if this line is perpendicular to a vector.
         * @param _vector The vector to compare to.
         * @return Returns true if this line is perpendicular to a vector, else returns false.
         */
        [[nodiscard]] bool is_perpendicular(const vector3& _vector) const;

        /**
         * Returns the angle between 2 lines.
         * @param _line The line to find the angle with.
         * @return The angle between 2 lines.
         */
        [[nodiscard]] float angle_between(const line& _line) const;

        /**
         * Returns the angle between this line and a vector.
         * @param _vector The vector to find the angle with.
         * @return The angle between this line and a vector.
         */
        [[nodiscard]] float angle_between(const vector3& _vector) const;

        /**
         * Checks if a point lines on this line.
         * @param _point The point to check.
         * @return Returns true if the point lies on this line, else returns false.
         */
        [[nodiscard]] bool contains(const vector3& _point) const;

        /**
         * Returns the point on this line closest to a given point.
         * @param _point The given point.
         * @return The point on this line closest to a given point.
         */
        [[nodiscard]] vector3 closest_point(const vector3& _point) const;
    };
}