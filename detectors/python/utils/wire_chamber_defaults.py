# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Example file that will generate a toy detector metadata
# --------------------------------------------------------------------------

from impl import metadata
from impl import Shape, Material, SurfaceAccelerator

""" Types that are typically needed for wirechamber-like """


def add_wire_chamber_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):

    # Sensitive shapes
    metadata.add_sensitive(Shape.DRIFT_CELL)
    metadata.add_sensitive(Shape.STRAW_TUBE)

    # Cylindrical volume portals (barrel and endcap)
    metadata.add_portal(Shape.CONCENTRIC_CYLINDER)
    metadata.add_portal(Shape.RING)

    # Surface acceleration structure for the wires
    metadata.add_accel_structure(SurfaceAccelerator.CONCENTRIC_CYLINDER_GRID2D)

    if use_mat_maps:
        # Sensitive material
        metadata.add_material(Material.ROD)
        # Map the material for the support structures
        metadata.add_material(Material.CONCENTIRC_CYLINDER_MAP2D)
        metadata.add_material(Material.DISC_MAP2D)
        # Experimetal: map the all of the material above into 3D bins
        # metadata.add_material(Material.CYLINDER3D)

    if use_homogeneous_mat:
        # Sensitive material
        metadata.add_material(Material.ROD)
        # Model the gas content
        metadata.add_material(Material.VOLUME)
