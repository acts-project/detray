# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape
from utils import add_logging_options, parse_logging_options, add_wire_chamber_defaults

import logging

# --------------------------------------------------------------------------
# Generate the wire chamber metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the wire chamber test detector """


def add_wire_chamber_types(metadata: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the wire chamber detector:")

    # Add default types for wire chambers
    add_wire_chamber_defaults(
        metadata=metadata, use_homogeneous_mat=True, use_mat_maps=False
    )

    # Beampipe passive surface
    logger.info("-> adding wire chamber beampipe")
    metadata.add_passive(Shape.CYLINDER2D)


def __main__():
    # Commandline options
    parser = add_logging_options()
    parse_logging_options(parser.parse_args())

    md = metadata("wire_chamber_test_detector")

    add_wire_chamber_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
