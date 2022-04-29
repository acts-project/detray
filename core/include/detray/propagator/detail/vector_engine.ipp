/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

template <typename scalar_t>
typename detray::detail::vector_engine<scalar_t>::point3
detray::detail::vector_engine<scalar_t>::pos(
    const detray::free_vector& free_vec) const {
    return {matrix_operator().element(free_vec, e_free_pos0, 0),
            matrix_operator().element(free_vec, e_free_pos1, 0),
            matrix_operator().element(free_vec, e_free_pos2, 0)};
}

template <typename scalar_t>
void detray::detail::vector_engine<scalar_t>::set_pos(
    detray::free_vector& free_vec, const point3& pos) {
    matrix_operator().element(free_vec, e_free_pos0, 0) = pos[0];
    matrix_operator().element(free_vec, e_free_pos1, 0) = pos[1];
    matrix_operator().element(free_vec, e_free_pos2, 0) = pos[2];
}

template <typename scalar_t>
typename detray::detail::vector_engine<scalar_t>::vector3
detray::detail::vector_engine<scalar_t>::dir(
    const detray::free_vector& free_vec) const {
    return {matrix_operator().element(free_vec, e_free_dir0, 0),
            matrix_operator().element(free_vec, e_free_dir1, 0),
            matrix_operator().element(free_vec, e_free_dir2, 0)};
}

template <typename scalar_t>
void detray::detail::vector_engine<scalar_t>::set_dir(
    detray::free_vector& free_vec, const vector3& dir) {
    matrix_operator().element(free_vec, e_free_dir0, 0) = dir[0];
    matrix_operator().element(free_vec, e_free_dir1, 0) = dir[1];
    matrix_operator().element(free_vec, e_free_dir2, 0) = dir[2];
}

template <typename scalar_t>
typename detray::detail::vector_engine<scalar_t>::point3
detray::detail::vector_engine<scalar_t>::local(
    const detray::bound_vector& bound_vec) const {
    return {matrix_operator().element(bound_vec, e_bound_loc0, 0),
            matrix_operator().element(bound_vec, e_bound_loc1, 0), 0};
}

template <typename scalar_t>
typename detray::detail::vector_engine<scalar_t>::vector3
detray::detail::vector_engine<scalar_t>::dir(
    const detray::bound_vector& bound_vec) const {

    const auto& phi = matrix_operator().element(bound_vec, e_bound_phi, 0);
    const auto& theta = matrix_operator().element(bound_vec, e_bound_theta, 0);
    const auto cosTheta = std::cos(theta);
    const auto sinTheta = std::sin(theta);

    return {std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta};
}

template <typename scalar_t>
detray::bound_vector
detray::detail::vector_engine<scalar_t>::free_to_bound_vector(
    const transform3& trf, const detray::free_vector& free_vec) const {

    const point3 pos = this->pos(free_vec);
    const vector3 dir = this->dir(free_vec);

    const auto local = trf.point_to_local(pos);

    bound_vector bound_vec;
    getter::element(bound_vec, e_bound_loc0, 0) = local[0];
    getter::element(bound_vec, e_bound_loc1, 0) = local[1];
    getter::element(bound_vec, e_bound_phi, 0) = getter::phi(dir);
    getter::element(bound_vec, e_bound_theta, 0) = getter::theta(dir);
    getter::element(bound_vec, e_bound_time, 0) =
        matrix_operator().element(free_vec, e_free_time, 0);
    getter::element(bound_vec, e_bound_qoverp, 0) =
        matrix_operator().element(free_vec, e_free_qoverp, 0);

    return bound_vec;
}

template <typename scalar_t>
detray::free_vector
detray::detail::vector_engine<scalar_t>::bound_to_free_vector(
    const transform3& trf, const detray::bound_vector& bound_vec) const {

    const vector3 local = this->local(bound_vec);

    const auto pos = trf.point_to_global(local);

    const vector3 dir = this->dir(bound_vec);

    free_vector free_vec;
    getter::element(free_vec, e_free_pos0, 0) = pos[0];
    getter::element(free_vec, e_free_pos1, 0) = pos[1];
    getter::element(free_vec, e_free_pos2, 0) = pos[2];
    getter::element(free_vec, e_free_time, 0) =
        matrix_operator().element(bound_vec, e_bound_time, 0);
    getter::element(free_vec, e_free_dir0, 0) = dir[0];
    getter::element(free_vec, e_free_dir1, 0) = dir[1];
    getter::element(free_vec, e_free_dir2, 0) = dir[2];
    getter::element(free_vec, e_free_qoverp, 0) =
        matrix_operator().element(bound_vec, e_bound_qoverp, 0);

    return free_vec;
}