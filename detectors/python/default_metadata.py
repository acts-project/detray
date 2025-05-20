# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape, Material, SurfaceAccelerator
from utils import (
    add_calorimeter_defaults,
    add_silicon_tracker_defaults,
    add_telescope_detector_defaults,
    add_wire_chamber_defaults,
)

from itk_metadata import add_itk_types
from odd_metadata import add_odd_types

# --------------------------------------------------------------------------
# Generate the default metadata type (can represent all detectors)
# --------------------------------------------------------------------------


def __main__():
    # Collect the types required for a detector called "test_det"
    md = metadata("default_new")

    # Make sure all of the defaults are added
    add_silicon_tracker_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_calorimeter_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_telescope_detector_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_wire_chamber_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)

    # Make sure all detectors can be represented by this metadata
    add_odd_types(md)
    add_itk_types(md)

    #
    # Add some special types (e.g. for the test detectors)
    #

    # Special surface types
    # md.add_sensitive(Shape.UNBOUNDED)  # Always hit the surface
    # md.add_sensitive(Shape.UNMASKED)

    # Add more material types
    # Material maps on non-concentric cylinders
    # md.add_material(Material.CYLINDER_MAP2D)

    # Add surface acceleration structures
    # md.add_accel_struct(SurfaceAccelerator.CYLINDER_GRID3D)
    # md.add_accel_struct(SurfaceAccelerator.CUBOID_GRID3D)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
