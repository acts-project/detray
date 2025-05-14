# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Example file that will generate a toy detector metadata
# --------------------------------------------------------------------------

from impl import metadata
from impl import Shape, Material

""" Types that are typically needed for calorimeters """


def add_calorimeter_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):

    metadata.add_sensitive(Shape.RECTANGLE)
    metadata.add_sensitive(Shape.TRAPEZOID)

    metadata.add_portal(Shape.CONCENTRIC_CYLINDER)
    metadata.add_portal(Shape.RING)

    if use_mat_maps:
        metadata.add_material(Material.CYLINDER_MAP3D)
    if use_homogeneous_mat:
        metadata.add_material(Material.VOLUME)

    # TODO: Add acceleration structures (e.g. Frustum navigation)
