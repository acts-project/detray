#pragma once

#include "Math/SMatrix.h"
#include "Math/SVector.h"

#include <any>
#include <array>
#include <tuple>
#include <cmath>

#ifdef DETRAY_CUSTOM_SCALARTYPE
using detray_scalar = DETRAY_CUSTOM_SCALARTYPE;
#else
using detray_scalar = double;
#endif

#define __plugin smatrix

namespace detray
{
    using scalar = detray_scalar;

    using namespace ROOT::Math;

    // smatrix getter methdos
    namespace getter
    {
        /** This method retrieves phi from a vector, vector base with rows > 2
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto phi(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM >= 2, "vector::phi() required rows >= 2.");
            return std::atan2(v[0], v[1]);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto theta(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM > 2, "vector::theta() required rows > 2.");
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto norm(const SVector<scalar, kDIM> &v)
        {
            return ROOT::Math::Dot(v, v);
        }

        /** This method retrieves the pseudo-rapidity from a vector or vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto eta(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM > 2, "vector::eta() required rows > 2.");
            return std::atanh(v[2] / norm(v));
        }

        /** This method retrieves the perpenticular magnitude of a vector with rows >= 2
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto perp(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM >= 2, "vector::perp() required rows >= 2.");
            return std::sqrt(v[0] * v[0] + v[1] * v[1]);
        }
        
        /** This method retrieves a column from a matrix
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS, typename matrix_type>
        auto vector(const matrix_type &m, unsigned int row, unsigned int col)
        {
            return m.template SubCol<SVector<scalar, kROWS> > (col, row);
        }

        /** This method retrieves a column from a matrix
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
        auto block(const matrix_type &m, unsigned int row, unsigned int col)
        {
            return m.template Sub <SMatrix<scalar, kROWS, kCOLS> >(row, col);
        }

    } // namespace getter

    // eigen definitions
    namespace smatrix
    {
        /** Transform wrapper class to ensure standard API within differnt plugins
         * 
         **/
        struct transform3
        {
            using vector3 = SVector<scalar, 3>;
            using point3 = vector3;
            using context = std::any;

            SMatrix<scalar, 4, 4> _data = ROOT::Math::SMatrixIdentity();
            SMatrix<scalar, 4, 4> _data_inv = ROOT::Math::SMatrixIdentity();

            using matrix44 = decltype(_data);
            using matrix33 = SMatrix<scalar, 3, 3>;

            /** Contructor with arguments: t, z, x, ctx
             * 
             * @param t the translation (or origin of the new frame)
             * @param z the z axis of the new frame, normal vector for planes
             * @param x the x axis of the new frame
             * 
             **/
            transform3(const vector3 &t, const vector3 &z, const vector3 &x, const context & /*ctx*/)
            {
                auto y = Cross(z, x);
                _data(0, 0) = x[0];
                _data(1, 0) = x[1];
                _data(2, 0) = x[2];
                _data(0, 1) = y[0];
                _data(1, 1) = y[1];
                _data(2, 1) = y[2];
                _data(0, 2) = z[0];
                _data(1, 2) = z[1];
                _data(2, 2) = z[2];
                _data(0, 3) = t[0];
                _data(1, 3) = t[1];
                _data(2, 3) = t[2];
                
                int ifail = 0;
                _data_inv = _data.Inverse(ifail);
                // This should be an exception, "transform3 could not be initialized. Matrix not invertible."
                assert(ifail == 0);
            }

            /** Constructor with arguments: translation
             *
             * @param t is the translation
             **/
            transform3(const vector3 &t, const context & /*ctx*/)
            {
                _data(0, 3) = t[0];
                _data(1, 3) = t[1];
                _data(2, 3) = t[2];
                
                int ifail = 0;
                _data_inv = _data.Inverse(ifail);
                // This should be an exception, "transform3 could not be initialized. Matrix not invertible."
                assert(ifail == 0);
            }

            /** Constructor with arguments: matrix 
             * 
             * @param m is the full 4x4 matrix 
             **/
            transform3(const matrix44 &m, const context & /*ctx*/)
            {
                _data = m;
                
               int ifail = 0;
                _data_inv = _data.Inverse(ifail);
                // This should be an exception, "transform3 could not be initialized. Matrix not invertible."
                assert(ifail == 0);
            }

            /** Constructor with arguments: matrix as std::aray of scalar
             * 
             * @param ma is the full 4x4 matrix asa 16 array
             **/
            transform3(const std::array<scalar, 16> &ma, const context & /*ctx*/)
            {

                _data(0, 0) = ma[0];
                _data(1, 0) = ma[4];
                _data(2, 0) = ma[8];
                _data(3, 0) = ma[12];
                _data(0, 1) = ma[1];
                _data(1, 1) = ma[5];
                _data(2, 1) = ma[9];
                _data(3, 1) = ma[13];
                _data(0, 2) = ma[2];
                _data(1, 2) = ma[6];
                _data(2, 2) = ma[10];
                _data(3, 2) = ma[14];
                _data(0, 3) = ma[3];
                _data(1, 3) = ma[7];
                _data(2, 3) = ma[11];
                _data(3, 3) = ma[15];
                
                int ifail = 0;
                _data_inv = _data.Inverse(ifail);
                // This should be an exception, "transform3 could not be initialized. Matrix not invertible."
                assert(ifail == 0);
            }

            /** Default contructors */
            transform3(const transform3 &rhs) = default;
            ~transform3() = default;

            /** This method retrieves the rotation of a transform
             * 
             * @param ctx the context object
             * 
             * @note this is a contextual method
             **/
            auto rotation(const context & /*ctx*/) const
            {
                return (_data.Sub<SMatrix<scalar, 3, 3> >(0, 0));
            }

            /** This method retrieves the translation of a transform
             * 
             * @param ctx the context object
             * 
             * @note this is a contextual method
             **/
            auto translation(const context & /*ctx*/) const
            {
                return (_data.SubCol<SVector<scalar, 3> >(3, 0));
            }

            /** This method retrieves the 4x4 matrix of a transform
             * 
             * @param ctx the context object
             * 
             * @note this is a contextual method
             **/
            auto matrix(const context & /*ctx*/) const
            {
                return _data;
            }
                     
            /** This method retrieves the translation of a transform
             * 
             * @param ctx the context object
             * 
             * @note this is a contextual method
             **/
            auto translation_inv(const context & /*ctx*/) const
            {
                return (_data_inv.SubCol<SVector<scalar, 3> >(3, 0));
            }
            
            /** This method retrieves the rotation of a transform
             * 
             * @param ctx the context object
             * 
             * @note this is a contextual method
             **/
            auto rotation_inv(const context & /*ctx*/) const
            {
                return (_data_inv.Sub<SMatrix<scalar, 3, 3> >(0, 0));
            }

            /** This method transform from a point from the local 3D cartesian frame to the global 3D cartesian frame
             * 
             * @note this is a contextual method 
             **/
            const point3 point_to_global(const point3 &v, const smatrix::transform3::context & ctx) const
            {
                return translation(ctx) + (rotation_inv(ctx) * v);
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame
             * 
             * @note this is a contextual method 
             **/
            const point3 point_to_local(const point3 &v, const smatrix::transform3::context & ctx) const
            {
                return translation_inv(ctx) + (rotation_inv(ctx) * v);
            }

            /** This method transform from a vector from the local 3D cartesian frame to the global 3D cartesian frame
             * 
             * @note this is a contextual method 
             **/
            const point3 vector_to_global(const vector3 &v, const smatrix::transform3::context & ctx) const
            {
                return rotation(ctx) * v;
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame
             * 
             * @note this is a contextual method 
             **/
            const point3 vector_to_local(const vector3 &v, const smatrix::transform3::context & ctx) const
            {
                return rotation_inv(ctx) * v;
            }
        };

        /** Non-contextual local frame projection into a cartesian coordinate frame
         */
        struct cartesian2
        {
            using point2 = SVector<scalar, 2>;

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame,
              * including the contextual transform into the local 3D frame
              * 
              * @tparam the type of the surface from which also point3 and context type can be deduced
              * 
              */
            template <typename surface_type>
            const auto operator()(const surface_type &s,
                                  const typename surface_type::transform3::point3 &p,
                                  const typename surface_type::transform3::context &ctx) const
            {
                return operator()(s.transform().point_to_local(p, ctx));
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             */
            template<typename point3_type> const auto operator()(const point3_type &v) const
            {
                return v.template Sub<SVector<scalar, 2> >(0);
            }
        };

        /** Non-contextual local frame projection into a polar coordinate frame
         **/
        struct polar2
        {
            using point2 = SVector<scalar, 2>;

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame,
              * including the contextual transform into the local 3D frame
              * 
              * @tparam the type of the surface from which also point3 and context type can be deduced
              * 
              */
            template <typename surface_type>
            const auto operator()(const surface_type &s,
                                  const typename surface_type::transform3::point3 &p,
                                  const typename surface_type::transform3::context &ctx) const
            {
                return operator()(s.transform().point_to_local(p, ctx));
            }

            /** This method transform from a point from 2D or 3D cartesian frame to a 2D polar point */
            template <typename point3_type>
            const auto operator()(const point3_type &v) const
            {
                return point2{getter::perp(v), getter::phi(v)};
            }
        };

        /** Non-contextual local frame projection into a polar coordinate frame
         **/
        struct cylindrical2
        {
            using point2 = SVector<scalar, 2>;

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame,
              * including the contextual transform into the local 3D frame
              * 
              * @tparam the type of the surface from which also point3 and context type can be deduced
              * 
              */
            template <typename surface_type>
            const auto operator()(const surface_type &s,
                                  const typename surface_type::transform3::point3 &p,
                                  const typename surface_type::transform3::context &ctx) const
            {
                return operator()(s.transform().point_to_local(p, ctx));
            }

            /** This method transform from a point from 2 3D cartesian frame to a 2D cylindrical point */
            template <typename point3_type>
            const auto operator()(const point3_type &v) const
            {
                return point2{getter::perp(v) * getter::phi(v), v[2]};
            }
        };

    } // namespace eigen

    // Non-contextual vector transfroms
    namespace vector
    {

        /** Get a normalized version of the input vector
         * 
         * @tparam derived_type is the matrix template
         * 
         * @param v the input vector
         **/
        template <typename vector_type>
        auto normalize(const vector_type &v)
        {
            return ROOT::Math::Unit(v);
        }

        /** Dot product between two input vectors
         * 
         * @tparam derived_type_lhs is the first matrix (epresseion) template
         * @tparam derived_type_rhs is the second matrix (epresseion) template
         * 
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return the scalar dot product value 
         **/        
        template <typename vector3_type, typename vecexpr3_type>
        auto dot(const vector3_type &a, const vecexpr3_type &b)
        {
            return ROOT::Math::Dot(a, b);
        }

        /** Cross product between two input vectors
         * 
         * @tparam derived_type_lhs is the first matrix (epresseion) template
         * @tparam derived_type_rhs is the second matrix (epresseion) template
         *           
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return a vector (expression) representing the cross product
         **/
        template <typename vector3_type>
        auto cross(const vector3_type &a, const vector3_type &b)
        {
            return ROOT::Math::Cross(a, b);
        }

    } // namespace vector

} // namespace detray
