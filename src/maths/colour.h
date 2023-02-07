#pragma once

namespace mkr {
    /**
     * Represents a colour with the components RGBA.
     * Each component is usually within the range of 0 to 1.
     */
    class colour {
    public:
        static const colour red;
        static const colour green;
        static const colour blue;
        static const colour yellow;
        static const colour cyan;
        static const colour magenta;
        static const colour black;
        static const colour white;
        static const colour grey;
        static const colour dark_grey;

        static constexpr colour constexpr_red() { return colour{1.0f, 0.0f, 0.0f, 1.0f}; }
        static constexpr colour constexpr_green() { return colour{0.0f, 1.0f, 0.0f, 1.0f}; }
        static constexpr colour constexpr_blue() { return colour{0.0f, 0.0f, 1.0f, 1.0f}; }
        static constexpr colour constexpr_yellow() { return colour{1.0f, 1.0f, 0.0f, 1.0f}; }
        static constexpr colour constexpr_cyan() { return colour{0.0f, 1.0f, 1.0f, 1.0f}; }
        static constexpr colour constexpr_magenta() { return colour{1.0f, 0.0f, 1.0f, 1.0f}; }
        static constexpr colour constexpr_black() { return colour{0.0f, 0.0f, 0.0f, 1.0f}; }
        static constexpr colour constexpr_white() { return colour{1.0f, 1.0f, 1.0f, 1.0f}; }
        static constexpr colour constexpr_grey() { return colour{0.25f, 0.25f, 0.25f, 1.0f}; }
        static constexpr colour constexpr_dark_grey() { return colour{0.1f, 0.1f, 0.1f, 1.0f}; }

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
