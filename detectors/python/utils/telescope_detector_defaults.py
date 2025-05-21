# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Example file that will generate a toy detector metadata
# --------------------------------------------------------------------------

from impl import metadata
from impl import Shape, Material, Accelerator

""" Types that are typically needed for telescope detectors """


def add_telescope_detector_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):

    # Sensitive and portal shapes
    metadata.add_sensitive(Shape.RECTANGLE)
    metadata.add_portal(Shape.RECTANGLE)

    # Map the material for the support structures onto the portals
    if use_mat_maps:
        metadata.add_material(Material.RECTANGLE_MAP2D)

    # Sensitive material
    if use_homogeneous_mat:
        metadata.add_material(Material.SLAB)

    # Acceleration struct for portals and passives
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "portal")
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    # Acceleration struct for telescope detector volumes
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "volume")
