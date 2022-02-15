
#include <array>
#include <map>
#include <tuple>

#include "algebra/vc_vc.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "vecmem/containers/jagged_vector.hpp"
#include "vecmem/containers/vector.hpp"

#define __plugin algebra::vc
#define ALGEBRA_PLUGIN vc_array

namespace detray {

using scalar = DETRAY_CUSTOM_SCALARTYPE;

template <typename value_type, std::size_t kDIM>
using darray = std::array<value_type, kDIM>;

template <typename value_type>
using dvector = vecmem::vector<value_type>;

template <typename value_type>
using djagged_vector = vecmem::jagged_vector<value_type>;

template <typename key_type, typename value_type>
using dmap = std::map<key_type, value_type>;

template <class... types>
using dtuple = std::tuple<types...>;

namespace getter = algebra::getter;
namespace vector = algebra::vector;

template <typename T, std::size_t ROWS, std::size_t COLS>
using matrix = __plugin::matrix_type<T, ROWS, COLS>;
template <typename T, std::size_t N>
using sym_matrix = matrix<T, N, N>;

}  // namespace detray
