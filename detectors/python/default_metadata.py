# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Type, Algebra, Shape, Material, Accelerator
from utils import (
    add_logging_options,
    parse_logging_options,
    add_calorimeter_defaults,
    add_silicon_tracker_defaults,
    add_telescope_detector_defaults,
    add_wire_chamber_defaults,
)

from itk_metadata import add_itk_types
from odd_metadata import add_odd_types

import logging

# --------------------------------------------------------------------------
# Generate the default metadata type (can represent all detectors)
# --------------------------------------------------------------------------


def __main__():
    # Commandline options
    parser = add_logging_options()
    parse_logging_options(parser.parse_args())

    logger = logging.getLogger(__name__)

    # Collect the types required for a detector that can hold all detector types
    md = metadata("default")

    # Specify a particular algebra plugin (otherwise left as template param.)
    # md.set_algebra_plugin(Algebra.ARRAY, Type.SINGLE)

    # Make sure all of the defaults are added
    add_silicon_tracker_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_calorimeter_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_telescope_detector_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_wire_chamber_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)

    # Make sure all detectors can be represented by this metadata
    add_odd_types(md)
    add_itk_types(md)

    # Add an acceleration struct for general 2D cylinders
    md.add_accel_struct(Accelerator.CYLINDER_GRID2D)
    md.add_material(Material.CYLINDER_MAP2D)

    #
    # Add some special types (e.g. for detector R&D)

    # Special surface types
    # md.add_sensitive(Shape.UNBOUNDED)  # Always hit the surface
    # md.add_sensitive(Shape.UNMASKED)

    # Add more material types
    # Material maps on non-concentric cylinders
    # md.add_material(Material.CYLINDER_MAP2D)

    # Add surface acceleration structures
    # md.add_accel_struct(Accelerator.CYLINDER_GRID3D)
    # md.add_accel_struct(Accelerator.CUBOID_GRID3D)

    # Make sure a default acceleration struct is chosen that can be used in all
    # detector types
    logger.info("Overwrite default acceleration structures to most generic type")
    md.set_default_accel_struct(Accelerator.BRUTE_FORCE, "portal")
    md.set_default_accel_struct(Accelerator.CYLINDER_GRID3D, "volume")

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
