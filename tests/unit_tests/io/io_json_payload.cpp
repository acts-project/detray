/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/json/json_serializers.hpp"

// GTest include(s)
#include <gtest/gtest.h>

/// This tests the json io for the general file header information
TEST(io, json_header_payload) {

    detray::header_payload<bool> h;
    h.common.version = "v0.0.1";
    h.common.detector = "test_detector";
    h.common.tag = "test_file";
    h.common.date = "01.01.2023";

    nlohmann::ordered_json j;
    j["header"] = h;

    detray::header_payload<bool> ph = j["header"];

    EXPECT_EQ(h.common.version, ph.common.version);
    EXPECT_EQ(h.common.detector, ph.common.detector);
    EXPECT_EQ(h.common.tag, ph.common.tag);
    EXPECT_EQ(h.common.date, ph.common.date);
}

/// This tests the json io for a single index link
TEST(io, single_link_payload) {
    detray::single_link_payload sl;
    sl.link = 3u;

    nlohmann::ordered_json j;
    j["single_link"] = sl;

    detray::single_link_payload psl = j["single_link"];

    EXPECT_EQ(sl.link, psl.link);
}

/// This tests the json io for a transform3
TEST(io, json_algebra_payload) {

    detray::transform_payload p;
    p.tr = {100.f, 200.f, 300.f};
    p.rot = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};

    nlohmann::ordered_json j;
    j["transform"] = p;

    detray::transform_payload pt = j["transform"];

    EXPECT_EQ(p.tr, pt.tr);
    EXPECT_EQ(p.rot, pt.rot);
}

/// This tests the json io for a grid axis
TEST(io, json_axis_payload) {

    detray::axis_payload ea;
    ea.binning = detray::n_axis::binning::e_regular;
    ea.bounds = detray::n_axis::bounds::e_circular;
    ea.label = detray::n_axis::label::e_phi;
    ea.edges = {-detray::constant<detray::real_io>::pi,
                detray::constant<detray::real_io>::pi};
    ea.bins = 10UL;

    nlohmann::ordered_json je;
    je["axis"] = ea;

    detray::axis_payload pea = je["axis"];

    EXPECT_EQ(ea.binning, pea.binning);
    EXPECT_EQ(ea.bounds, pea.bounds);
    EXPECT_EQ(ea.edges, pea.edges);

    EXPECT_EQ(ea.bins, pea.bins);

    detray::axis_payload va;
    va.binning = detray::n_axis::binning::e_irregular;
    va.bounds = detray::n_axis::bounds::e_closed;
    va.label = detray::n_axis::label::e_r;
    va.edges = {0.f, 1.f, 4.f, 5.f, 8.f, 10.f};
    va.bins = va.edges.size() - 1UL;

    nlohmann::ordered_json jv;
    jv["axis"] = va;

    detray::axis_payload pva = jv["axis"];

    EXPECT_EQ(va.binning, pva.binning);
    EXPECT_EQ(va.bounds, pva.bounds);
    EXPECT_EQ(va.label, pva.label);
    EXPECT_EQ(va.edges, pva.edges);
    EXPECT_EQ(va.bins, pva.bins);
}

/// This tests the json io for a grid bin
TEST(io, json_bin_payload) {

    detray::grid_bin_payload b;
    b.loc_index = std::vector<unsigned int>{1u, 0u, 2u};
    b.content = std::vector<std::size_t>{0u, 1u, 2u, 3u};

    nlohmann::ordered_json jbin;
    jbin["bin"] = b;

    detray::grid_bin_payload pb = jbin["bin"];

    EXPECT_EQ(b.loc_index.size(), pb.loc_index.size());
    EXPECT_EQ(b.content.size(), pb.content.size());
}

/// This tests the json io for a grid
TEST(io, json_grid_payload) {

    std::vector<detray::grid_bin_payload> bins = {
        {{0u, 1u}, {0u, 2u}}, {{1u, 1u}, {1u, 2u}}, {{2u, 1u}, {2u, 2u}}};

    detray::axis_payload a0{
        detray::n_axis::binning::e_regular, detray::n_axis::bounds::e_circular,
        detray::n_axis::label::e_phi, 3u,
        std::vector<detray::real_io>{-detray::constant<detray::real_io>::pi,
                                     detray::constant<detray::real_io>::pi}};

    detray::axis_payload a1{
        detray::n_axis::binning::e_regular, detray::n_axis::bounds::e_closed,
        detray::n_axis::label::e_r, 2u, std::vector<detray::real_io>{0.f, 2.f}};

    detray::grid_payload g;
    g.acc_link = {detray::grid_payload::grid_type::polar2_grid, 12u};
    g.volume_link = {2u};
    g.axes = {a0, a1};
    g.bins = bins;

    nlohmann::ordered_json j;
    j["grid"] = g;

    detray::grid_payload pg = j["grid"];

    EXPECT_EQ(g.acc_link.type, pg.acc_link.type);
    EXPECT_EQ(g.acc_link.index, pg.acc_link.index);
    EXPECT_EQ(g.axes.size(), pg.axes.size());
    EXPECT_EQ(g.bins.size(), pg.bins.size());
}

/// This tests the json io for a surface mask
TEST(io, json_mask_payload) {

    detray::single_link_payload sl;
    sl.link = 3u;

    detray::mask_payload m;
    m.shape = detray::mask_payload::mask_shape::cylinder3;
    m.volume_link = sl;
    m.boundaries = {10.f, 100.f};

    nlohmann::ordered_json j;
    j["mask"] = m;

    detray::mask_payload pm = j["mask"];

    EXPECT_EQ(m.shape, pm.shape);
    EXPECT_EQ(m.volume_link.link, pm.volume_link.link);
    EXPECT_EQ(m.boundaries, pm.boundaries);
}

/// This tests the json io for a surface material link
TEST(io, json_material_link_payload) {

    detray::material_link_payload m;
    m.type = detray::material_link_payload::type_id::slab;
    m.index = 2u;

    nlohmann::ordered_json j;
    j["material"] = m;

    detray::material_link_payload pm = j["material"];

    EXPECT_EQ(m.type, pm.type);
    EXPECT_EQ(m.index, pm.index);
}

/// This tests the json payload for a surface (descriptor + data)
TEST(io, json_surface_payload) {

    detray::surface_payload s;

    detray::transform_payload t;
    t.tr = {100.f, 200.f, 300.f};
    t.rot = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};

    detray::mask_payload m;
    detray::single_link_payload sl;
    sl.link = 1u;
    m.shape = detray::mask_payload::mask_shape::trapezoid2;
    m.volume_link = sl;
    m.boundaries = {10.f, 20.f, 34.f, 1.4f};

    detray::material_link_payload mat;
    mat.type = detray::material_link_payload::type_id::slab;
    mat.index = 2u;

    s.transform = t;
    s.mask = m;
    s.type = detray::surface_id::e_passive;
    s.material = mat;

    nlohmann::ordered_json j;
    j["surface"] = s;

    detray::surface_payload ps = j["surface"];

    EXPECT_EQ(s.transform.tr, ps.transform.tr);
    EXPECT_EQ(s.transform.rot, ps.transform.rot);

    EXPECT_EQ(s.mask.shape, ps.mask.shape);
    EXPECT_EQ(s.mask.volume_link.link, ps.mask.volume_link.link);
    EXPECT_EQ(s.mask.boundaries, ps.mask.boundaries);

    EXPECT_EQ(s.type, ps.type);

    EXPECT_EQ(s.material.value().type, ps.material.value().type);
    EXPECT_EQ(s.material.value().index, ps.material.value().index);
}

/// This tests the json io for a surface material link
TEST(io, acc_links_payload) {

    detray::acc_links_payload l;
    l.type = detray::acc_links_payload::type_id::cylinder2_grid;
    l.index = 2u;

    nlohmann::ordered_json j;
    j["acc_link"] = l;

    detray::acc_links_payload pl = j["acc_link"];

    EXPECT_EQ(l.type, pl.type);
    EXPECT_EQ(l.index, pl.index);
}

/// This tests the json payload for a volume (descriptor + data (transform,
/// surfaces)
TEST(io, json_volume_payload) {

    detray::transform_payload t;
    t.tr = {100.f, 200.f, 300.f};
    t.rot = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};

    detray::single_link_payload sl;
    sl.link = 1u;

    detray::acc_links_payload al;
    al.type = detray::acc_links_payload::type_id::cylinder2_grid;
    al.index = 2u;

    detray::surface_payload s;

    detray::mask_payload m;
    m.shape = detray::mask_payload::mask_shape::trapezoid2;
    m.volume_link = sl;
    m.boundaries = {10.f, 20.f, 34.f, 1.4f};

    detray::material_link_payload mat;
    mat.type = detray::material_link_payload::type_id::slab;
    mat.index = 2u;

    s.transform = t;
    s.mask = m;
    s.type = detray::surface_id::e_portal;
    s.material = mat;

    detray::volume_payload v;
    v.name = "volume";
    v.type = detray::volume_id::e_cylinder;
    sl.link = 2u;
    v.index = sl;
    v.transform = t;
    v.surfaces = {s};
    v.acc_links = {al};

    nlohmann::ordered_json j;
    j["volume"] = v;

    detray::volume_payload pv = j["volume"];

    EXPECT_EQ(v.name, pv.name);
    EXPECT_EQ(v.index.link, pv.index.link);
    EXPECT_EQ(v.transform.tr, v.transform.tr);
    EXPECT_EQ(v.transform.rot, v.transform.rot);
    EXPECT_EQ(v.type, pv.type);
    EXPECT_EQ(v.surfaces.size(), pv.surfaces.size());
    EXPECT_EQ(v.acc_links->size(), pv.acc_links->size());
}

/// This tests the json io for a material slab/rod
TEST(io, json_material_slab_payload) {

    detray::material_slab_payload m;
    m.mat_link = {detray::material_slab_payload::type::slab, 21u};
    m.thickness = 1.2f;
    m.mat.params = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f};

    nlohmann::ordered_json j;
    j["material"] = m;

    detray::material_slab_payload pm = j["material"];

    EXPECT_EQ(m.mat_link.type, pm.mat_link.type);
    EXPECT_EQ(m.mat_link.index, pm.mat_link.index);
    EXPECT_EQ(m.thickness, pm.thickness);
    EXPECT_EQ(m.mat.params, pm.mat.params);
}

/// This tests the json io for a material slab
TEST(io, json_detector_payload) {

    detray::detector_payload d;
    d.volumes = {detray::volume_payload{}, detray::volume_payload{}};

    nlohmann::ordered_json j;
    j["detector"] = d;

    detray::detector_payload pd = j["detector"];

    EXPECT_EQ(d.volumes.size(), pd.volumes.size());
}
