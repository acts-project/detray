#include "algebra/definitions/vc_array.hpp"

#define ALGEBRA_PLUGIN vc_array

namespace detray {

    using scalar = algebra::scalar;

    template <typename value_type, unsigned int kDIM>
    using darray = algebra::array_t<value_type, kDIM>;

    template <typename value_type>
    using dvector = algebra::vector_t<value_type>;

    template <typename key_type, typename value_type>
    using dmap = algebra::map_t<key_type, value_type>;

    template< class... types>
    using dtuple = algebra::tuple_t<types ...>;

    namespace getter = algebra::getter;
    namespace vector = algebra::vector;

} //namespace detray
