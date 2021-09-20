
#include "algebra/definitions/eigen.hpp"

#define ALGEBRA_PLUGIN eigen

namespace detray {

    using scalar = algebra::scalar;

    template <typename value_type, unsigned int kDIM>
    using darray = std::array<value_type, kDIM>;

    template <typename value_type>
    using dvector = vecmem::vector<value_type>;

    template <typename key_type, typename value_type>
    using dmap = algebra::map_s<key_type, value_type>;

    template< class... types>
    using dtuple = algebra::tuple_s<types ...>;

    namespace getter = algebra::getter;
    namespace vector = algebra::vector;

} //namespace detray
