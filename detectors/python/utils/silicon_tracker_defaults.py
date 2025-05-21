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

""" Types that are typically needed for silicon tracker detectors """


def add_silicon_tracker_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):

    # Barrel Detector
    metadata.add_sensitive(Shape.RECTANGLE)
    metadata.add_accel_struct(Accelerator.CONCENTRIC_CYLINDER_GRID2D, "sensitive")
    if use_mat_maps:
        metadata.add_material(Material.CONCENTIRC_CYLINDER_MAP2D)

    # Endcap Detector
    metadata.add_sensitive(Shape.TRAPEZOID)
    metadata.add_accel_struct(Accelerator.DISC_GRID2D, "sensitive")
    if use_mat_maps:
        metadata.add_material(Material.DISC_MAP2D)

    # Slabs can be used for both barrel and endcap surface shapes
    if use_homogeneous_mat:
        metadata.add_material(Material.SLAB)

    # Cylindrical volume portals (barrel and endcap)
    metadata.add_portal(Shape.CONCENTRIC_CYLINDER)
    metadata.add_portal(Shape.RING)

    # Acceleration struct for portals and passives
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "portal")
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    # Volume accelerator for layered cylindrical detectors
    metadata.add_accel_struct(Accelerator.CYLINDER_GRID3D, "volume")
