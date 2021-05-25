
#include "algebra/containers/array.hpp"

namespace detray {

    using scalar = algebra::scalar;
    template <typename value_type, unsigned int kDIM>
    using darray = algebra::array_t<value_type, kDIM>;
    template <typename value_type>
    using dvector = algebra::vector_t<value_type>;

    using algebra::operator*;
    using algebra::operator+;
    using algebra::operator-;

    namespace getter = algebra::getter;
    namespace vector = algebra::vector;

} //namespace detray
