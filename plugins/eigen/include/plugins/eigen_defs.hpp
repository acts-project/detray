#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <any>
#include <tuple>
#include <cmath>

namespace detray
{
    // eigen definitions
    namespace eigen
    {
        // Primitives definition
        using scalar = float;
        using vector2 = Eigen::Matrix<scalar, 2, 1>;
        using point2 = vector2;
        using point2pol = Eigen::Matrix<scalar, 2, 1>;
        using point2cyl = Eigen::Matrix<scalar, 2, 1>;
        using vector3 = Eigen::Matrix<scalar, 3, 1>;
        using point3 = vector3;

        /** Transform wrapper class to ensure standard API within differnt plugins
         * 
         **/
        struct transform3
        {

            using context = std::any;

            Eigen::Transform<scalar, 3, Eigen::Affine> _data =
                Eigen::Transform<scalar, 3, Eigen::Affine>::Identity();

            using matrix44 = Eigen::Transform<scalar, 3, Eigen::Affine>::MatrixType;

            /** The contextual transform interface, ignored for the moment
             **/
            auto contextual(const context & /*ignored*/) const
            {
                return _data;
            }

            /** Contructor with arguments: t, z, x, ctx
             * 
             * @param t the translation (or origin of the new frame)
             * @param z the z axis of the new frame, normal vector for planes
             * @param x the x axis of the new frame
             * 
             **/
            transform3(const vector3 &t, const vector3 &z, const vector3 &x, const context & /*ctx*/)
            {
                auto y = z.cross(x);

                auto &matrix = _data.matrix();
                matrix.block<3, 1>(0, 0) = x;
                matrix.block<3, 1>(0, 1) = y;
                matrix.block<3, 1>(0, 2) = z;
                matrix.block<3, 1>(0, 3) = t;
            }

            /** Constructor with arguments: translation
             *
             * @param t is the transform
             **/
            transform3(const vector3 &t, const context & /*ctx*/)
            {
                auto &matrix = _data.matrix();
                matrix.block<3, 1>(0, 3) = t;
            }

            /** Constructor with arguments: matrix 
             * 
             * @param mat is the full 4x4 matrix 
             **/
            transform3(const matrix44 &m, const context & /*ctx*/)
            {
                _data.matrix() = m;
            }
        };

    } // namespace eigen

    // Getter methdos
    namespace getter
    {
        /** This method retrieves phi from a vector, vector base with rows > 2
         * 
         * @param v the input vector 
         **/
        template <typename derived_type>
        eigen::scalar phi(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::phi() required rows >= 2.");
            return std::atan2(v[0], v[1]);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <typename derived_type>
        eigen::scalar theta(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::theta() required rows >= 3.");
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the pseudo-rapidity from a vector or vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <typename derived_type>
        eigen::scalar eta(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::eta() required rows >= 3.");
            return std::atanh(v[2] / v.norm());
        }

        /** This method retrieves the perpenticular magnitude of a vector with rows >= 2
         * 
         * @param v the input vector 
         **/
        template <typename derived_type>
        eigen::scalar perp(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::perp() required rows >= 2.");
            return std::sqrt(v[0] * v[0] + v[1] * v[1]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         * 
         * @param v the input vector 
         **/
        template <typename derived_type>
        eigen::scalar norm(const Eigen::MatrixBase<derived_type> &v)
        {
            return v.norm();
        }

        /** This method retrieves the rotation of a transform
         * 
         * @param trf the input transform 
         * @param ctx the context object
         * 
         * @note this is a contextual method
         **/
        auto rotation(const eigen::transform3 &trf, const eigen::transform3::context &ctx)
        {
            return trf.contextual(ctx).matrix().block<3, 3>(0, 0).eval();
        }

      /** This method retrieves the translation of a transform
         * 
         * @param trf the input transform 
         * @param ctx the context object
         * 
         * @note this is a contextual method
         **/
        auto translation(const eigen::transform3 &trf, const eigen::transform3::context &ctx)
        {
            return trf.contextual(ctx).matrix().block<3, 1>(0, 3).eval();
        }

       /** This method retrieves the 4x4 matrix of a transform
         * 
         * @param trf the input transform 
         * @param ctx the context object
         * 
         * @note this is a contextual method
         **/
        auto matrix(const eigen::transform3 &trf, const eigen::transform3::context &ctx)
        {
            return trf.contextual(ctx).matrix();
        }


        /** This method retrieves a column from a matrix
         * 
         * @param m the input matrix 
         **/
        template <typename kCOLS, typename kROWS, typename derived_type>
        auto block(const Eigen::MatrixBase<derived_type> &m, unsigned int col, unsigned int row)
        {
            return m.template block<kCOLS,kROWS>(col, row);
        }


    } // namespace getter

    // Non-contextual vector transfroms
    namespace vector
    {

        template <typename derived_type>
        auto normalize(const Eigen::MatrixBase<derived_type> &v)
        {
            return v.normalized();
        }

        template <typename derived_type>
        auto dot(const Eigen::MatrixBase<derived_type> &a, const Eigen::MatrixBase<derived_type> &b)
        {
            return a.dot(b);
        }

        template <typename derived_type>
        auto cross(const Eigen::MatrixBase<derived_type> &a, const Eigen::MatrixBase<derived_type> &b)
        {
            return a.cross(b);
        }
    } // namespace vector

    // Contextual and non-contextual transform operations
    namespace transform
    {

        /** This method transform from a point from the local 3D cartesian frame to the global 3D cartesian frame
         * 
         * @note this is a contextual method 
         * */
        template <typename derived_type>
        const auto lpoint3_to_gpoint3(const eigen::transform3 &trf, const Eigen::MatrixBase<derived_type> &v, const eigen::transform3::context &ctx)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::lpoint3_to_gpoint3(v) requires a (3,1) matrix");
            return trf.contextual(ctx) * v;
        }

        /** This method transform from a vector from the local 3D cartesian frame to the global 3D cartesian frame
        * 
        * @note this is a contextual method 
        * */
        template <typename derived_type>
        const auto lvector3_to_gvector3(const eigen::transform3 &trf, const Eigen::MatrixBase<derived_type> &v, const eigen::transform3::context &ctx)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::lvector3_to_gvector3(v) requires a (3,1) matrix");
            return (trf.contextual(ctx).linear() * v).eval();
        }

        /** This method transform from a point from the global 3D cartesian frame to the local 3D cartesian frame
         * 
         * @note this is a contextual method 
         * */
        template <typename derived_type>
        const auto gpoint3_to_lpoint3(const eigen::transform3 &trf, const Eigen::MatrixBase<derived_type> &v, const eigen::transform3::context &ctx)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::gpoint3_to_lpoint3(v) requires a (3,1) matrix");
            return (trf.contextual(ctx).inverse() * v).eval();
        }

        /** This method transform from a vector from the global 3D cartesian frame to the global 3D cartesian frame
        * 
        * @note this is a contextual method 
        * */
        template <typename derived_type>
        const auto gvector3_to_lvector3(const eigen::transform3 &trf, const Eigen::MatrixBase<derived_type> &v, const eigen::transform3::context &ctx)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::gvector3_to_lvector3(v) requires a (3,1) matrix");
            return (trf.contextual(ctx).inverse().linear() * v).eval();
        }

        /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
        * */
        template <typename derived_type>
        const auto point3_to_point2(const Eigen::MatrixBase<derived_type> &v)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::point3_to_point2(v) requires a (3,1) matrix");
            return v.template segment<2>(0);
        }

        /** This method transform from a point from 2D or 3D cartesian frame to a 2D polar point */
        template <typename derived_type>
        eigen::point2pol point_to_point2pol(const Eigen::MatrixBase<derived_type> &v)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows >= 2 and cols == 1, "transform::point_to_point2pol(v) requires a (>2,1) matrix");
            return {getter::perp(v), getter::phi(v)};
        }

        /** This method transform from a point from 2 3D cartesian frame to a 2D cylindrical point */
        template <typename derived_type>
        const eigen::point2cyl point3_to_point2cyl(const Eigen::MatrixBase<derived_type> &v)
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
            static_assert(rows == 3 and cols == 1, "transform::point3_to_point2cyl(v) requires a a (3,1) matrix");
            return {getter::perp(v) * getter::phi(v), v[2]};
        }

    } // namespace transform

} // namespace detray
