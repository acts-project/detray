#pragma once

#include "SMatrix.h"
#include "SVector.h"

#include <any>
#include <cmath>

namespace detray
{
    // eigen definitions
    namespace smatrix
    {
        // Primitives definition
        using scalar = float;
        using vector2 = SVector<scalar, 2>;
        using point2 = vector2;
        using point2pol = SVector<scalar, 2>;
        using point2cyl = SVector<scalar, 2>;
        using vector3 = SVector<scalar, 3>;
        using point3 = SVector<scalar, 3>;

        // Context definition
        using context = std::any;

        /** Transform wrapper class to ensure standard API within differnt plugins
         **/
        struct transform3
        {

            ROOT::Math::SMatrix<scalar, 4, 4> _data = ROOT::Math : SMatrix<scalar, 3, 3>(SMatrixIdentity);

            using matrix44 = decltype(_data);
            using matrix33 = ROOT::Math::SMatrix<scalar, 3, 3>;

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
                _data.Place_in_col<smatrix::scalar, 3>(x, 0, 0);
                _data.Place_in_col<smatrix::scalar, 3>(y, 0, 1);
                _data.Place_in_col<smatrix::scalar, 3>(z, 0, 2);
                _data.Place_in_col<smatrix::scalar, 3>(t, 0, 3);
            }

            /** Constructor with arguments: translation
             *
             * @param t is the transform
             **/
            transform3(const vector3 &t, const context & /*ctx*/)
            {
                _data.Place_in_col<3, 1>(t, 0, 3);
            }

            /** Constructor with arguments: matrix 
             * 
             * @param mat is the full 4x4 matrix 
             **/
            transform3(const matrix44 &m, const context & /*ctx*/)
            {
                _data = m;
            }
        };

    } // namespace smatrix

    // Getter methdos
    namespace getter
    {
        /** This method retrieves phi from a vector, vector base with rows > 2
         * 
         * @param v the input vector 
         **/
        template <typename kROWS>
        smatrix::scalar phi(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v) noexcept
        {
            static_assert(kROWS >= 2, "vector::phi() required rows >= 2.");
            return std::atan2(v[0], v[1]);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <typename kROWS>
        smatrix::scalar theta(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v) noexcept
        {
            static_assert(kROWS >= 2, "vector::theta() required rows >= 3.");
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the pseudo-rapidity from a vector or vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <typename kROWS>
        smatrix::scalar eta(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v) noexcept
        {
            static_assert(kROWS >= 2, "vector::eta() required rows >= 3.");
            return std::atanh(v[2] / v.norm());
        }

        /** This method retrieves the perpenticular magnitude of a vector with rows >= 2
         * 
         * @param v the input vector 
         **/
        template <typename kROWS>
        smatrix::scalar perp(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v) noexcept
        {
            static_assert(kROWS >= 2, "vector::perp() required rows >= 2.");
            return std::sqrt(v[0] * v[0] + v[1] * v[1]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         * 
         * @param v the input vector 
         **/
        template <typename kROWS>
        smatrix::scalar norm(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v)
        {
            return ROOT::Math::Mag(v);
        }

        /** This method retrieves rotation of a transform
         * 
         * @param trf the input transform 
         * @param ctx the context object
         * 
         * @note this is a contextual method
         **/
        auto rotation(const smatrix::transform3 &trf, const smatrix::context &ctx)
        {
            return trf.contextual(ctx).Sub<smatrix::matrix33>(0, 0);
        }

        auto translation(const smatrix::transform3 &trf, const smatrix::context &ctx)
        {
            return trf.contextual(ctx).Sub<smatrix::vector3>(0, 3);
        }

        auto matrix(const smatrix::transform3 &trf, const smatrix::context &ctx)
        {
            return trf.contextual(ctx);
        }

    } // namespace getter

    // Non-contextual vector transfroms
    namespace vector
    {

        template <typename kROWS>
        auto normalize(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v)
        {
            return v.Unit();
        }

        template <typename kROWS>
        auto dot(const ROOT::Math::SVector<smatrix::scalar, kROWS> &a, const ROOT::Math::SVector < smatrix::scalar, kROWS> &b)
        {
            return a.dot(b);
        }

        template <typename kROWS>
        auto cross(const ROOT::Math::SVector<smatrix::scalar, kROWS> &a, const ROOT::Math::SVector<smatrix::scalar, kROWS> &b)
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
        const auto lpoint3_to_gpoint3(const smatrix::transform3 &trf, const smatrix::point3 &v, const smatrix::context &ctx)
        {
            auto cmatrix = getter::matrix(trf, ctx);
            return cmatrix.Sub<transform3::matrix33>(cmatrix, 0,0) * v + cmatrix.Sub<smatrix::vector3>(cmatrix, 0, 3);
        }

        /** This method transform from a vector from the local 3D cartesian frame to the global 3D cartesian frame
        * 
        * @note this is a contextual method 
        * */
        const auto lvector3_to_gvector3(const smatrix::transform3 &trf, const smatrix::vector3 &v, const smatrix::context &ctx)
        {
            static_assert(kROWS == 3, "transform::lvector3_to_gvector3(v) requires a (3,1) matrix");
            return (trf.contextual(ctx).linear() * v).eval();
        }

        /** This method transform from a point from the global 3D cartesian frame to the local 3D cartesian frame
         * 
         * @note this is a contextual method 
         * */
        const auto gpoint3_to_lpoint3(const smatrix::transform3 &trf, const smatrix::point3 &v, const smatrix::context &ctx)
        {
            return (trf.contextual(ctx).inverse() * v).eval();
        }

        /** This method transform from a vector from the global 3D cartesian frame to the global 3D cartesian frame
        * 
        * @note this is a contextual method 
        * */
        const auto gvector3_to_lvector3(const smatrix::transform3 &trf, const smatrix::vector3 &v, const smatrix::context &ctx)
        {
            return (trf.contextual(ctx).inverse().linear() * v).eval();
        }

        /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
        * */
        const auto point3_to_point2(const smatrix::point3 &v)
        {
            static_assert(kROWS == 3, "transform::point3_to_point2(v) requires a (3,1) matrix");
            return v.template segment<2>(0);
        }

        /** This method transform from a point from 2D or 3D cartesian frame to a 2D polar point */
        template <typename kROWS>
        smatrix::point2pol point_to_point2pol(const ROOT::Math::SVector<smatrix::scalar, kROWS> &v)
        {
            static_assert(kROWS >= 2, "transform::point_to_point2pol(v) requires a (>2,1) matrix");
            return {getter::perp(v), getter::phi(v)};
        }

        /** This method transform from a point from 2 3D cartesian frame to a 2D cylindrical point */
        const smatrix::point2cyl point3_to_point2cyl(const smatrix::point3 &v)
        {
            static_assert(kROWS, "transform::point3_to_point2cyl(v) requires a a (3,1) matrix");
            return {getter::perp(v) * getter::phi(v), v[2]};
        }

    } // namespace transform

} // namespace detray
