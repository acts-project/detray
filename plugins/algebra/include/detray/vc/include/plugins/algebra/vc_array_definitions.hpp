#include "algebra/definitions/vc_array.hpp"
#include "vecmem/containers/jagged_vector.hpp"

#define ALGEBRA_PLUGIN vc_array
#define __vc__

namespace detray {

using scalar = algebra::scalar;

template <typename value_type, unsigned int kDIM>
using darray = std::array<value_type, kDIM>;

template <typename value_type>
using dvector = algebra::vector_v<value_type>;

template <typename value_type>
using djagged_vector = algebra::vector_v<algebra::vector_v<value_type>>;

template <typename key_type, typename value_type>
using dmap = algebra::map_s<key_type, value_type>;

template <class... types>
using dtuple = algebra::tuple_s<types...>;

using algebra::operator*;
using algebra::operator/;
using algebra::operator+;
using algebra::operator-;

namespace getter = algebra::getter;
namespace vector = algebra::vector;

}  // namespace detray
