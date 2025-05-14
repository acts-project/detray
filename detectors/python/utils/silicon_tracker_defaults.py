# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Function that will add types commonly needed for silicon trackers
# --------------------------------------------------------------------------

from impl import metadata
from impl import Shape, Material, Accelerator

import logging

""" Types that are typically needed for silicon tracker detectors """


def add_silicon_tracker_defaults(
    metadata: metadata, use_mat_maps=False, use_homogeneous_mat=False
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if metadata in add_silicon_tracker_defaults.clients:
        return

    add_silicon_tracker_defaults.clients.append(metadata)

    logger = logging.getLogger(__name__)
    logger.info("Define silicon tracker types:")

    # Cylindrical volume portals (barrel and endcap)
    logger.info("-> adding portal types")
    metadata.add_portal(Shape.CONCENTRIC_CYLINDER)
    metadata.add_portal(Shape.RING)

    # Acceleration struct for portals and passives
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "portal", is_default=True)
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    # Barrel Detector
    logger.info("-> adding barrel section types")
    metadata.add_sensitive(Shape.RECTANGLE)
    metadata.add_accel_struct(Accelerator.CONCENTRIC_CYLINDER_GRID2D, "sensitive")
    if use_mat_maps:
        metadata.add_material(Material.CONCENTIRC_CYLINDER_MAP2D)

    # Endcap Detector
    logger.info("-> adding endcap section types")
    metadata.add_sensitive(Shape.TRAPEZOID)
    metadata.add_accel_struct(Accelerator.DISC_GRID2D, "sensitive")
    if use_mat_maps:
        metadata.add_material(Material.DISC_MAP2D)

    # Slabs can be used for both barrel and endcap surface shapes
    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        metadata.add_material(Material.SLAB)

    # Volume accelerator for layered cylindrical detectors
    logger.info("-> adding detector volume acceleration structure")
    metadata.add_accel_struct(Accelerator.CYLINDER_GRID3D, "volume", is_default=True)


add_silicon_tracker_defaults.clients = []
