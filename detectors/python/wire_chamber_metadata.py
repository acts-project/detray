# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape
from utils import add_wire_chamber_defaults

# --------------------------------------------------------------------------
# Generate the wire chamber metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the wire chamber test detector """


def add_wire_chamber_types(metadata: metadata):
    # Add default types for wire chambers
    add_wire_chamber_defaults(
        metadata=metadata, use_homogeneous_mat=True, use_mat_maps=False
    )

    # Beampipe passive surface
    metadata.add_passive(Shape.CYLINDER2D)


def __main__():
    md = metadata("wire_chamber_test_detector")

    add_wire_chamber_types(md)

    # Dump the metadata to header file
    metadata_generator(md)


if __name__ == "__main__":
    __main__()
