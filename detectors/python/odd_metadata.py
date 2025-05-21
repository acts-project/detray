# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape
from utils import add_silicon_tracker_defaults

# --------------------------------------------------------------------------
# Generate the Open Data Detector metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ACTS Open Data Detector (ODD) """


def add_odd_types(metadata: metadata):
    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=metadata, use_homogeneous_mat=True, use_mat_maps=True
    )

    # Beampipe passive surface
    metadata.add_passive(Shape.CYLINDER2D)


def __main__():
    md = metadata("odd_test")

    add_odd_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
