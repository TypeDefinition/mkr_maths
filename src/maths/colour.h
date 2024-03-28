#pragma once

namespace mkr {
    /**
     * Represents a colour with the components RGBA.
     * Each component is usually within the range of 0 to 1.
     */
    class colour {
    public:
        static constexpr colour red() { return colour{1.0f, 0.0f, 0.0f, 1.0f}; }
        static constexpr colour green() { return colour{0.0f, 1.0f, 0.0f, 1.0f}; }
        static constexpr colour blue() { return colour{0.0f, 0.0f, 1.0f, 1.0f}; }
        static constexpr colour yellow() { return colour{1.0f, 1.0f, 0.0f, 1.0f}; }
        static constexpr colour cyan() { return colour{0.0f, 1.0f, 1.0f, 1.0f}; }
        static constexpr colour magenta() { return colour{1.0f, 0.0f, 1.0f, 1.0f}; }
        static constexpr colour black() { return colour{0.0f, 0.0f, 0.0f, 1.0f}; }
        static constexpr colour white() { return colour{1.0f, 1.0f, 1.0f, 1.0f}; }
        static constexpr colour grey() { return colour{0.25f, 0.25f, 0.25f, 1.0f}; }
        static constexpr colour dark_grey() { return colour{0.1f, 0.1f, 0.1f, 1.0f}; }

        /// Red component.
        float r_;
        /// Green component.
        float g_;
        /// Blue component.
        float b_;
        /// Alpha component.
        float a_;

        /**
         * Constructs the colour.
         * @param _r The red component of the colour.
         * @param _g The green component of the colour.
         * @param _b The blue component of the colour.
         * @param _a The alpha component of the colour.
         */
        constexpr colour(float _r = 1.0f, float _g = 1.0f, float _b = 1.0f, float _a = 1.0f)
                : r_(_r), g_(_g), b_(_b), a_(_a) {}

        bool operator==(const colour& _rhs) const;

        bool operator!=(const colour& _rhs) const;

        colour operator+(const colour& _rhs) const;

        colour& operator+=(const colour& _rhs);

        colour operator-(const colour& _rhs) const;

        colour& operator-=(const colour& _rhs);

        colour operator*(const colour& _rhs) const;

        colour& operator*=(const colour& _rhs);

        colour operator*(float _scalar) const;

        colour& operator*=(float _scalar);

        friend colour operator*(float _scalar, const colour& _colour);
    };
}
