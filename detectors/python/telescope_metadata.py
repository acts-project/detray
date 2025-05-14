# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from impl import metadata, metadata_generator
from impl import Shape, Material
from utils import add_telescope_detector_defaults

# --------------------------------------------------------------------------
# Generate the metadata type for the detray test telescope detector
# --------------------------------------------------------------------------


def __main__():
    # Collect the types required for a detector called "test_det"
    md = metadata("telescope")

    # Cuboid volume based detector shape
    add_telescope_detector_defaults(md)

    # Add more geometric shape types for telescope test detector
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
