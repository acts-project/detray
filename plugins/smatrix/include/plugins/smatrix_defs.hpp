#pragma once

#include "utils/containers.hpp"

#include "Math/SMatrix.h"
#include "Math/SVector.h"

#include <any>
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
        template <typename vec_expr_type>
        auto phi(const vec_expr_type &v) noexcept
        {
            //static_assert(kDIM >= 2, "vector::perp() required rows >= 2.");
            scalar element0 = v.apply(0);
            scalar element1 = v.apply(1);
            return std::atan2(element1, element0);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto theta(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM > 2, "vector::theta() required rows >= 3.");
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto norm(const SVector<scalar, kDIM> &v)
        {
            return std::sqrt(ROOT::Math::Dot(v, v));
        }

        /** This method retrieves the pseudo-rapidity from a vector or vector base with rows >= 3
         * 
         * @param v the input vector 
         **/
        template <unsigned int kDIM>
        auto eta(const SVector<scalar, kDIM> &v) noexcept
        {
            static_assert(kDIM > 2, "vector::eta() required rows >= 3.");
            return std::atanh(v[2] / norm(v));
        }

        /** This method retrieves the perpenticular magnitude of a vector with rows >= 2
         * 
         * @param v the input vector 
         **/
        template <typename vec_expr_type>
        auto perp(const vec_expr_type &v) noexcept
        {
            //static_assert(kDIM >= 2, "vector::perp() required rows >= 2.");
            scalar element0 = v.apply(0);
            scalar element1 = v.apply(1);
            return std::sqrt(element0 * element0 + element1 * element1);
        }

        /** This method retrieves a column from a matrix
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS, typename matrix_type>
        auto vector(const matrix_type &m, unsigned int row, unsigned int col)
        {
            return m.template SubCol<SVector<scalar, kROWS>>(col, row);
        }

        /** This method retrieves a column from a matrix
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
        auto block(const matrix_type &m, unsigned int row, unsigned int col)
        {
            return m.template Sub<SMatrix<scalar, kROWS, kCOLS>>(row, col);
        }

    } // namespace getter

    // eigen definitions
    namespace smatrix
    {
        using vector3 = SVector<scalar, 3>;
        using point3 = vector3;
        using vector2 = SVector<scalar, 2>;

        /** Transform wrapper class to ensure standard API within differnt plugins
         * 
         **/
        struct transform3
        {

            SMatrix<scalar, 4, 4> _data = ROOT::Math::SMatrixIdentity();
            SMatrix<scalar, 4, 4> _data_inv = ROOT::Math::SMatrixIdentity();

            using matrix44 = decltype(_data);
            using matrix33 = SMatrix<scalar, 3, 3>;

            /** Contructor with arguments: t, z, x
             * 
             * @param t the translation (or origin of the new frame)
             * @param z the z axis of the new frame, normal vector for planes
             * @param x the x axis of the new frame
             * 
             **/
            transform3(const vector3 &t, const vector3 &z, const vector3 &x)
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
            transform3(const vector3 &t)
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
            transform3(const matrix44 &m)
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
            transform3(const std::array<scalar, 16> &ma)
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
            transform3() = default;
            transform3(const transform3 &rhs) = default;
            ~transform3() = default;

            /** Equality operator */
            bool operator==(const transform3 &rhs) const
            {
                return _data == rhs._data;
            }

            /** This method retrieves the rotation of a transform */
            auto rotation() const
            {
                return (_data.Sub<SMatrix<scalar, 3, 3>>(0, 0));
            }

            /** This method retrieves the translation of a transform */
            auto translation() const
            {
                return (_data.SubCol<SVector<scalar, 3>>(3, 0));
            }

            /** This method retrieves the 4x4 matrix of a transform */
            auto matrix() const
            {
                return _data;
            }

            /** This method transform from a point from the local 3D cartesian frame to the global 3D cartesian frame */
            const point3 point_to_global(const point3 &v) const
            {
                SVector<scalar, 4> vector_4 = SVector<scalar, 4>();
                vector_4.Place_at(v, 0);
                vector_4[3] = static_cast<scalar>(1);
                return SVector<scalar, 4>(_data * vector_4).template Sub<SVector<scalar, 3>>(0);
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame */
            const point3 point_to_local(const point3 &v) const
            {
                SVector<scalar, 4> vector_4 = SVector<scalar, 4>();
                vector_4.Place_at(v, 0);
                vector_4[3] = static_cast<scalar>(1);
                return SVector<scalar, 4>(_data_inv * vector_4).template Sub<SVector<scalar, 3>>(0);
            }

            /** This method transform from a vector from the local 3D cartesian frame to the global 3D cartesian frame */
            const point3 vector_to_global(const vector3 &v) const
            {
                SVector<scalar, 4> vector_4 = SVector<scalar, 4>();
                vector_4.Place_at(v, 0);
                return SVector<scalar, 4>(_data * vector_4).template Sub<SVector<scalar, 3>>(0);
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame */
            const point3 vector_to_local(const vector3 &v) const
            {
                SVector<scalar, 4> vector_4 = SVector<scalar, 4>();
                vector_4.Place_at(v, 0);
                return SVector<scalar, 4>(_data_inv * vector_4).template Sub<SVector<scalar, 3>>(0);
            }
        };

        /** Local frame projection into a cartesian coordinate frame */
        struct cartesian2
        {
            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param v the point in local frame
             * 
             * @return a local point2
             */
            template <typename point3_type>
            const auto operator()(const point3_type &v) const
            {
                return v.template Sub<SVector<scalar, 2>>(0);
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame 
             * 
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             * 
             * @return a local point2
             **/
            const auto operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

        /** Local frame projection into a polar coordinate frame
         **/
        struct polar2
        {
            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param v the point in local frame
             * 
             * @return a local point2
             */
            template <typename point3_type>
            const auto operator()(const point3_type &v) const
            {
                return point2{getter::perp(v), getter::phi(v)};
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame 
             * 
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             * 
             * @return a local point2
             **/
            const auto operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

        /** Local frame projection into a polar coordinate frame
         **/
        struct cylindrical2
        {
            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param v the point in local frame
             * 
             * @return a local point2
             */
            template <typename point3_type>
            const auto operator()(const point3_type &v) const
            {
                return point2{getter::perp(v) * getter::phi(v), v[2]};
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame 
             * 
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             * 
             * @return a local point2
             **/
            const auto operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

    } // namespace smatrix

    // Vector transfroms
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
        template <typename vector3_type, typename vecexpr3_type>
        auto cross(const vecexpr3_type &a, const vector3_type &b)
        {
            return ROOT::Math::Cross(a, b);
        }

    } // namespace vector

} // namespace detray
