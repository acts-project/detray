/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

// This tests the construction of a detector class
TEST(detector, detector_kernel) {

    using namespace detray;

    vecmem::host_memory_resource host_mr;

    using detector_t =
        detector<detector_registry::default_detector, covfie::field>;
    using mask_ids = typename detector_t::masks::id;
    using material_ids = typename detector_t::materials::id;

    detector_t::geometry_context ctx0{};

    detector_t::transform_container trfs;
    detector_t::surface_container surfaces = {};
    detector_t::mask_container masks(host_mr);
    detector_t::material_container materials(host_mr);

    /// Surface 0
    point3 t0{0., 0., 0.};
    trfs.emplace_back(ctx0, t0);
    masks.template emplace_back<mask_ids::e_rectangle2>(empty_context{}, 0UL,
                                                        -3.f, 3.f);
    materials.template emplace_back<material_ids::e_slab>(empty_context{},
                                                          gold<scalar>(), 3.);

    /// Surface 1
    point3 t1{1., 0., 0.};
    trfs.emplace_back(ctx0, t1);
    masks.template emplace_back<mask_ids::e_annulus2>(
        empty_context{}, 0UL, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f);
    materials.template emplace_back<material_ids::e_slab>(
        empty_context{}, tungsten<scalar>(), 12.);

    /// Surface 2
    point3 t2{2., 0., 0.};
    trfs.emplace_back(ctx0, t2);
    masks.template emplace_back<mask_ids::e_trapezoid2>(empty_context{}, 0UL,
                                                        1.f, 2.f, 3.f);
    materials.template emplace_back<material_ids::e_rod>(
        empty_context{}, aluminium<scalar>(), 4.);

    detector_t d(host_mr);

    auto &v =
        d.new_volume(volume_id::e_cylinder, {0., 10., -5., 5., -M_PI, M_PI});
    d.add_objects_per_volume(ctx0, v, surfaces, masks, trfs, materials);
}

// This tests the construction of a detector class
TEST(detector, volume_builder) {

    using namespace detray;

    vecmem::host_memory_resource host_mr;

    using detector_t =
        detector<detector_registry::default_detector, covfie::field>;

    detector_t d(host_mr);

    EXPECT_TRUE(d.volumes().size() == 0);

    volume_builder<detector_t> vbuilder{};
    vbuilder.init_vol(d, volume_id::e_cylinder,
                      {0., 10., -5., 5., -M_PI, M_PI});
    const auto &vol = d.volumes().back();

    EXPECT_TRUE(d.volumes().size() == 1);
    EXPECT_EQ(vol.index(), 0);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);
}
