# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape
from utils import add_silicon_tracker_defaults

# --------------------------------------------------------------------------
# Generate the ITk metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ATLAS ITk detector """


def add_itk_types(metadata: metadata):
    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(metadata=metadata, use_mat_maps=True)

    # Strip stereo annulus shape
    metadata.add_sensitive(Shape.ANNULUS)

    # Beampipe passive surface
    metadata.add_passive(Shape.CYLINDER2D)


def __main__():
    md = metadata("itk")

    add_itk_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
