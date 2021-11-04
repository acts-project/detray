
#include <array>
#include <map>
#include <tuple>

#include "algebra/array_cmath.hpp"
#include "vecmem/containers/jagged_vector.hpp"
#include "vecmem/containers/vector.hpp"

#define __plugin algebra::array
#define ALGEBRA_PLUGIN array

namespace detray {

using scalar = algebra::scalar;

template <typename value_type, unsigned int kDIM>
using darray = std::array<value_type, kDIM>;

template <typename value_type>
using dvector = vecmem::vector<value_type>;

template <typename value_type>
using djagged_vector = vecmem::jagged_vector<value_type>;

template <typename key_type, typename value_type>
using dmap = std::map<key_type, value_type>;

template <class... types>
using dtuple = std::tuple<types...>;

using algebra::operator*;
using algebra::operator+;
using algebra::operator-;

namespace getter = algebra::getter;
namespace vector = algebra::vector;

}  // namespace detray
