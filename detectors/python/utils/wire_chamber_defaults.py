# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Function that will add types commonly needed for wire chamber detectors
# --------------------------------------------------------------------------

from impl import metadata
from impl import Shape, Material, Accelerator

import logging

""" Types that are typically needed for wirechamber-like detectors """


def add_wire_chamber_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if metadata in add_wire_chamber_defaults.clients:
        return

    add_wire_chamber_defaults.clients.append(metadata)

    logger = logging.getLogger(__name__)
    logger.info("Define wire chamber types:")

    # Sensitive shapes
    logger.info("-> adding sensitive types")
    metadata.add_sensitive(Shape.DRIFT_CELL)
    metadata.add_sensitive(Shape.STRAW_TUBE)

    # Cylindrical volume portals (barrel and endcap)
    logger.info("-> adding portal types")
    metadata.add_portal(Shape.CONCENTRIC_CYLINDER)
    metadata.add_portal(Shape.RING)

    # Acceleration struct for portals and passives
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "portal", is_default=True)
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    # Surface acceleration structure for the wires
    metadata.add_accel_struct(Accelerator.CONCENTRIC_CYLINDER_GRID2D, "sensitive")

    if use_mat_maps:
        logger.info("-> requested material map types")
        # Map the material for the support structures
        metadata.add_material(Material.CONCENTIRC_CYLINDER_MAP2D)
        metadata.add_material(Material.DISC_MAP2D)
        # Experimetal: map the all of the material above into 3D bins
        # metadata.add_material(Material.CYLINDER3D)

    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        # Model the gas content
        metadata.add_material(Material.RAW)

    # Sensitive material
    metadata.add_material(Material.ROD)

    # Volume accelerator for layered cylindrical detectors
    logger.info("-> adding detector volume acceleration structure")
    metadata.add_accel_struct(Accelerator.CYLINDER_GRID3D, "volume", is_default=True)


add_wire_chamber_defaults.clients = []
