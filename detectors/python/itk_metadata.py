# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import detray

from detray.detectors import metadata, metadata_generator
from detray.detectors import Shape
from detray.detectors import (
    add_logging_options,
    parse_logging_options,
    add_silicon_tracker_defaults,
)

import argparse
import logging
import sys

# --------------------------------------------------------------------------
# Generate the ITk metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ATLAS ITk detector """


def add_itk_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ATLAS Inner Tracker (ITk):")

    # Strip stereo annulus shape
    logger.info("-> adding ITk strip detecor custom shape")
    md.add_sensitive(Shape.ANNULUS)

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(metadata=md, use_mat_maps=True)

    # Beampipe passive surface
    logger.info("-> adding ITk beampipe")
    md.add_passive(Shape.CYLINDER2D)

    logger.info("Done")


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    detray.detectors.add_logging_options(parser)
    detray.detectors.parse_logging_options(parser.parse_args())

    md = metadata("itk")

    add_itk_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
