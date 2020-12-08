# detray

Detray is part of the ACTS project (R&D line for parallelization), the ACTS project can be found: https://github.com/acts-project/acts.

This is a C++17 header only library for detector surface intersections using different algebra plugin libraries. It follows the navigation and propagation concept of ACTS, however, with an attempt to create
a geometry without polymorphic inheritance structure.

### Repository structure

The code for defining detector surfaces is within the `core` directory, the algebra plugins are in the `plugin` directory, and tests are all within the `test` library.
