# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape
from utils import (
    add_logging_options,
    parse_logging_options,
    add_silicon_tracker_defaults,
)

import logging

# --------------------------------------------------------------------------
# Generate the Open Data Detector metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ACTS Open Data Detector (ODD) """


def add_odd_types(metadata: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ACTS Open Data Detector (ODD):")

    # Beampipe passive surface
    logger.info("-> adding ODD beampipe")
    metadata.add_passive(Shape.CYLINDER2D)

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=metadata, use_homogeneous_mat=False, use_mat_maps=True
    )


def __main__():
    # Commandline options
    parser = add_logging_options()
    parse_logging_options(parser.parse_args())

    md = metadata("open_data_detector")

    # Specify a particular algebra plugin (otherwise left as template param.)
    # md.set_algebra_plugin(Algebra.ARRAY, Type.SINGLE)

    add_odd_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
