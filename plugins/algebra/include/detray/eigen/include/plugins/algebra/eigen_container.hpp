
#include "algebra/containers/eigen.hpp"

namespace detray {

    using scalar = algebra::scalar;
    template <typename value_type, unsigned int kDIM>
    using darray = algebra::array_t<value_type, kDIM>;
    template <typename value_type>
    using dvector = algebra::vector_t<value_type>;

    namespace getter = algebra::getter;
    namespace vector = algebra::vector;

} //namespace detray
