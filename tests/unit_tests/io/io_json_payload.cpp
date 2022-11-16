/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/io/json_io.hpp"

TEST(io, json_algebra_payload) {

    detray::transform_payload p;
    p.tr = {100., 200., 300.};
    p.rot = {1., 0., 0., 0., 1., 0., 0., 0., 1};

    nlohmann::json j;
    j["transform"] = p;

    detray::transform_payload pt = j["transform"];

    EXPECT_EQ(p.tr, pt.tr);
    EXPECT_EQ(p.rot, pt.rot);
}

TEST(io, json_axis_payload) {

    detray::axis_payload ea;
    ea.type = detray::axis_payload::axis_type::equidistant;
    ea.bracket = detray::axis_payload::axis_bracket::closed;
    ea.lookup = detray::axis_payload::axis_lookup::phi;
    ea.borders = {-M_PI, M_PI};
    ea.bins = 10u;

    nlohmann::json je;
    je["axis"] = ea;

    detray::axis_payload pea = je["axis"];

    EXPECT_EQ(ea.type, pea.type);
    EXPECT_EQ(ea.bracket, pea.bracket);
    EXPECT_EQ(ea.borders, pea.borders);

    EXPECT_EQ(ea.bins, pea.bins);

    detray::axis_payload va;
    va.type = detray::axis_payload::axis_type::variable;
    va.bracket = detray::axis_payload::axis_bracket::bound;
    va.lookup = detray::axis_payload::axis_lookup::r;
    va.borders = {0, 1, 4, 5, 8, 10};
    va.bins = va.borders.size() - 1;

    nlohmann::json jv;
    jv["axis"] = va;

    detray::axis_payload pva = jv["axis"];

    EXPECT_EQ(va.type, pva.type);
    EXPECT_EQ(va.bracket, pva.bracket);
    EXPECT_EQ(va.lookup, pva.lookup);
    EXPECT_EQ(va.borders, pva.borders);
    EXPECT_EQ(va.bins, pva.bins);
}

TEST(io, json_grid_payload) {

    std::vector<std::vector<unsigned int>> entries = {{0, 1}, {0, 2}, {1, 1},
                                                      {1, 2}, {2, 1}, {2, 2}};

    detray::axis_payload a0{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::closed,
                            detray::axis_payload::axis_lookup::phi,
                            std::vector<detray::real_io>{-M_PI, M_PI}, 3u};

    detray::axis_payload a1{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::bound,
                            detray::axis_payload::axis_lookup::r,
                            std::vector<detray::real_io>{0, 2u}, 2u};

    detray::grid_payload g;
    g.axes = {a0, a1};
    g.entries = entries;

    nlohmann::json j;
    j["grid"] = g;

    detray::grid_payload pg = j["grid"];

    EXPECT_EQ(g.axes.size(), pg.axes.size());
    EXPECT_EQ(g.entries, pg.entries);
}

TEST(io, single_object_payload) {
    detray::single_object_payload so;
    so.link = 3u;

    nlohmann::json j;
    j["single_link"] = so;

    detray::single_object_payload pso = j["single_link"];

    EXPECT_EQ(so.link, pso.link);
}

TEST(io, grid_objects_payload) {
    detray::grid_objects_payload go;

    std::vector<std::vector<unsigned int>> entries = {{0, 1}, {0, 2}, {1, 1},
                                                      {1, 2}, {2, 1}, {2, 2}};

    detray::axis_payload a0{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::closed,
                            detray::axis_payload::axis_lookup::phi,
                            std::vector<detray::real_io>{-M_PI, M_PI}, 3u};

    detray::axis_payload a1{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::bound,
                            detray::axis_payload::axis_lookup::r,
                            std::vector<detray::real_io>{0, 2u}, 2u};

    detray::grid_payload g;
    g.axes = {a0, a1};
    g.entries = entries;

    go.grid = g;

    nlohmann::json j;
    j["links"] = go;
    detray::grid_objects_payload pgo = j["links"];

    EXPECT_EQ(go.grid.axes.size(), pgo.grid.axes.size());
    EXPECT_EQ(go.grid.entries, pgo.grid.entries);

    detray::transform_payload p;
    p.tr = {100., 200., 300.};
    p.rot = {1., 0., 0., 0., 1., 0., 0., 0., 1};

    detray::grid_objects_payload got;
    got.transform = p;

    got.grid = g;

    j["links_t"] = got;
    detray::grid_objects_payload pgot = j["links_t"];

    EXPECT_EQ(got.grid.axes.size(), pgot.grid.axes.size());
    EXPECT_EQ(got.grid.entries, pgot.grid.entries);
    EXPECT_EQ(got.transform.value().tr, pgot.transform.value().tr);
    EXPECT_EQ(got.transform.value().rot, pgot.transform.value().rot);
}

TEST(io, json_links_payload) {

    std::vector<std::vector<unsigned int>> entries = {{0, 1}, {0, 2}, {1, 1},
                                                      {1, 2}, {2, 1}, {2, 2}};

    detray::axis_payload a0{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::closed,
                            detray::axis_payload::axis_lookup::phi,
                            std::vector<detray::real_io>{-M_PI, M_PI}, 3u};

    detray::axis_payload a1{detray::axis_payload::axis_type::equidistant,
                            detray::axis_payload::axis_bracket::bound,
                            detray::axis_payload::axis_lookup::r,
                            std::vector<detray::real_io>{0, 2u}, 2u};

    detray::grid_payload g;
    g.axes = {a0, a1};
    g.entries = entries;

    detray::grid_objects_payload go;
    go.grid = g;

    detray::single_object_payload so;
    so.link = 3u;

    detray::links_payload l;
    l.grid_links = go;
    l.single_links = {so};

    nlohmann::json j;
    j["links"] = l;

    detray::links_payload pl = j["links"];

    EXPECT_EQ(l.single_links.size(), pl.single_links.size());
    EXPECT_EQ(l.grid_links.value().grid.axes.size(),
              pl.grid_links.value().grid.axes.size());
    EXPECT_EQ(l.grid_links.value().grid.entries,
              pl.grid_links.value().grid.entries);
}

TEST(io, json_mask_payload) {

    detray::mask_payload m;
    m.type = detray::mask_payload::mask_type::cylinder3;
    m.boundaries = {10., 100.};

    nlohmann::json j;
    j["mask"] = m;

    detray::mask_payload pm = j["mask"];

    EXPECT_EQ(m.type, pm.type);
    EXPECT_EQ(m.boundaries, pm.boundaries);
}

TEST(io, json_material_slab_payload) {

    detray::material_slab_payload m;
    m.slab = {1., 2., 3., 4., 5.};

    nlohmann::json j;
    j["material"] = m;

    detray::material_slab_payload pm = j["material"];

    EXPECT_EQ(m.slab, pm.slab);
}

TEST(io, json_surface_payload) {

    detray::surface_payload s;

    detray::transform_payload t;
    t.tr = {100., 200., 300.};
    t.rot = {1., 0., 0., 0., 1., 0., 0., 0., 1};

    detray::mask_payload m;
    m.type = detray::mask_payload::mask_type::trapezoid2;
    m.boundaries = {10., 20., 34., 1.4};

    detray::material_slab_payload mat;
    mat.slab = {1., 2., 3., 4., 5.};

    s.transform = t;
    s.mask = m;
    s.material = mat;

    nlohmann::json j;
    j["surface"] = s;

    detray::surface_payload ps = j["surface"];

    EXPECT_EQ(s.transform.tr, ps.transform.tr);
    EXPECT_EQ(s.transform.rot, ps.transform.rot);

    EXPECT_EQ(s.mask.type, ps.mask.type);
    EXPECT_EQ(s.mask.boundaries, ps.mask.boundaries);

    EXPECT_EQ(s.material.slab, ps.material.slab);
}

TEST(io, json_portal_payload) {

    detray::surface_payload s;

    detray::transform_payload t;
    t.tr = {100., 200., 300.};
    t.rot = {1., 0., 0., 0., 1., 0., 0., 0., 1};

    detray::mask_payload m;
    m.type = detray::mask_payload::mask_type::trapezoid2;
    m.boundaries = {10., 20., 34., 1.4};

    detray::material_slab_payload mat;
    mat.slab = {1., 2., 3., 4., 5.};

    s.transform = t;
    s.mask = m;
    s.material = mat;

    detray::single_object_payload so;
    so.link = 3u;

    detray::links_payload l;
    l.single_links = {so};

    detray::portal_payload p;
    p.surface = s;
    p.volume_links = l;

    nlohmann::json j;
    j["portal"] = p;

    detray::portal_payload pp = j["portal"];

    EXPECT_EQ(p.surface.transform.tr, pp.surface.transform.tr);
    EXPECT_EQ(p.surface.transform.rot, pp.surface.transform.rot);
    EXPECT_EQ(p.volume_links.single_links.size(),
              pp.volume_links.single_links.size());
}

TEST(io, json_volume_bounds_payload) {

    detray::volume_bounds_payload vb;
    vb.type = detray::volume_bounds_payload::volume_bounds_type::cylindrical;
    vb.values = {0., 100., 120.};

    nlohmann::json j;
    j["volume_bounds"] = vb;

    detray::volume_bounds_payload pvb = j["volume_bounds"];

    EXPECT_EQ(vb.type, pvb.type);
    EXPECT_EQ(vb.values, pvb.values);
}

TEST(io, json_volume_payload) {

    detray::transform_payload t;
    t.tr = {100., 200., 300.};
    t.rot = {1., 0., 0., 0., 1., 0., 0., 0., 1};

    detray::volume_bounds_payload vb;
    vb.type = detray::volume_bounds_payload::volume_bounds_type::cylindrical;
    vb.values = {0., 100., 120.};

    detray::single_object_payload so;
    so.link = 0u;

    detray::links_payload l;
    l.single_links = {so};

    detray::surface_payload s;

    detray::mask_payload m;
    m.type = detray::mask_payload::mask_type::trapezoid2;
    m.boundaries = {10., 20., 34., 1.4};

    detray::material_slab_payload mat;
    mat.slab = {1., 2., 3., 4., 5.};

    s.transform = t;
    s.mask = m;
    s.material = mat;

    detray::portal_payload p;
    p.surface = s;
    p.volume_links = l;

    detray::volume_payload v;
    v.name = "volume";
    v.transform = t;
    v.volume_bounds = vb;
    v.portals = {p};
    v.surfaces = {s};
    v.surface_links = l;

    nlohmann::json j;
    j["volume"] = v;

    detray::volume_payload pv = j["volume"];

    EXPECT_EQ(v.name, pv.name);
    EXPECT_EQ(v.surfaces.size(), pv.surfaces.size());
    EXPECT_EQ(v.portals.size(), pv.portals.size());
}

TEST(io, json_detector_payload) {

    detray::detector_payload d;
    d.name = "detector";
    d.volumes = {detray::volume_payload{}, detray::volume_payload{}};

    nlohmann::json j;
    j["detector"] = d;

    detray::detector_payload pd = j["detector"];

    EXPECT_EQ(d.name, pd.name);
    EXPECT_EQ(d.volumes.size(), pd.volumes.size());
}
