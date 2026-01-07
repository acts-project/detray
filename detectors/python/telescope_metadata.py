# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape, Material
from utils import (
    add_logging_options,
    parse_logging_options,
    add_telescope_detector_defaults,
)

import logging

# --------------------------------------------------------------------------
# Generate the metadata type for the detray telescope test detector
# --------------------------------------------------------------------------


def __main__():
    # Commandline options
    parser = add_logging_options()
    parse_logging_options(parser.parse_args())
    logger = logging.getLogger(__name__)

    # Collect the types required for a telescope test detector
    md = metadata("telescope_test_detector")

    # Cuboid volume based detector shape
    add_telescope_detector_defaults(md)

    # Add more geometric shape types for telescope test detector
    logger.info("-> adding additional telescope sensitive types")
    md.add_sensitive(Shape.ANNULUS)
    md.add_sensitive(Shape.CONCENTRIC_CYLINDER)
    md.add_sensitive(Shape.CYLINDER2D)
    md.add_sensitive(Shape.RING)
    md.add_sensitive(Shape.TRAPEZOID)
    md.add_sensitive(Shape.DRIFT_CELL)
    md.add_sensitive(Shape.STRAW_TUBE)
    md.add_sensitive(Shape.UNBOUNDED)  # Always hit the surface

    # Add more material types
    md.add_material(Material.ROD)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
