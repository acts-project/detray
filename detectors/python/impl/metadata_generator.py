# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project includes
from .type_helpers import link
from .definitions import (
    Type,
    Algebra,
    Shape,
    Material,
    SurfaceAccelerator,
    VolumeAccelerator,
)

# Python includes
from datetime import datetime
import logging


""" Class that represents the c++ metadata type with all of its fields """


class metadata:

    def __init__(self, detector_name):
        self.det_name = detector_name
        self.id_base = Type.UINT_8
        self.nav_link = link(link_type="single", data_type=Type.UINT_LEAST_16)
        self.mask_link = link(link_type="range", data_type=Type.UINT_32)
        self.material_link = link(link_type="single", data_type=Type.UINT_32)
        self.surface_types = []
        self.shapes = []
        self.materials = []
        self.acceleration_structs = [SurfaceAccelerator.BRUTE_FORCE]
        self.volume_finder = []

    # Register a new geometric shape for a portal surface
    def add_portal(self, shape: Shape):
        self.add_shape(shape)
        if "portal" not in self.surface_types:
            self.surface_types.append("portal")

    # Register a new geometric shape for a sensitive surface
    def add_sensitive(self, shape: Shape):
        self.add_shape(shape)
        if "sensitive" not in self.surface_types:
            self.surface_types.append("sensitive")

    # Register a new geometric shape for a passive surface
    def add_passive(self, shape: Shape):
        self.add_shape(shape)
        if "passive" not in self.surface_types:
            self.surface_types.append("passive")

    # Register a new geometric shape for the detector
    def add_shape(self, shape: Shape):
        if not shape in self.shapes:
            self.shapes.append(shape)

    # Register a new material type
    def add_material(self, mat: Material):
        if mat not in self.materials:
            self.materials.append(mat)

    # Register a new surface acceleration structure
    def add_accel_structure(self, accel: SurfaceAccelerator):
        if accel not in self.acceleration_structs:
            self.acceleration_structs.append(accel)


""" Class that accumulates detector type data and writes a metadat header """


class metadata_generator:

    def __init__(
        self,
        md: metadata,
        log_level=logging.INFO,
    ):
        # Internal state duriing header generation
        self.file = None
        self.logger = logging
        self.indent = 0

        # Save the type specifiers
        self.shape_specifiers = []
        self.material_specifiers = []
        self.accel_specifiers = []

        # Write log to terminal
        self.logger.basicConfig(
            format=("%(levelname)s (%(module)s): %(message)s"), level=log_level
        )

        # Dump the header
        self.__generate(md)

    # Write the current metadata to file
    def __generate(
        self, md: metadata, src_dir="../detectors/include/detray/detectors/"
    ):

        self.file = open(f"{src_dir}{md.det_name}_metadata.hpp", "w+")

        # Beginning of the header (copyright, includes etc.)
        self.__preamble(md.det_name, md.shapes, md.materials, md.acceleration_structs)

        # Basic typedefs
        self.__typedef("algebra_type", "algebra_t")
        self.__typedef("scalar_t", "dscalar<algebra_type>")
        self.__lines(1)
        self.__typedef("nav_link", md.nav_link.data_type)

        # Transform store
        self.__lines(1)
        self.__put(
            f"\
template <template <typename...> class vector_t = dvector>\n\
{self.__tabs()}using transform_store =\n\
{self.__tabs()}    single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;"
        )

        # Mask types
        assert md.shapes, "Define at least one geometric shape"
        self.__lines(2)
        for shape in md.shapes:
            self.__declare_mask(shape, "algebra_type", "nav_link", shape.param)

        self.__lines(2)
        self.__make_type_enum("mask_id", self.shape_specifiers, md.id_base)
        self.__lines(2)

        # Mask Store
        self.__declare_multi_store("mask_store", "mask_id", self.shape_specifiers)

        # Material types (optional)
        if md.materials:
            self.__lines(2)
            for mat in md.materials:
                self.__declare_material(mat, "algebra_type")

            self.__lines(1)
            # Add an option for 'no material'
            self.__make_type_enum(
                "material_id", self.material_specifiers, md.id_base, ["none"]
            )
        self.__lines(2)

        # Material Store
        self.__declare_multi_store(
            "material_store", "material_id", self.material_specifiers, is_regular=False
        )

        # Acceleration Structures (optional)
        if md.acceleration_structs:
            self.__lines(2)
            for acc in md.acceleration_structs:
                self.__declare_accel(acc, "algebra_type")

            self.__lines(1)
            # Add an option for 'no material'
            self.__make_type_enum("accel_id", self.accel_specifiers, md.id_base)
        self.__lines(2)

        # Accelerator Store
        self.__declare_multi_store(
            "accelerator_store", "accel_id", self.accel_specifiers, is_regular=False
        )

        self.__finish()

    # Add a string to the header file
    def __put(self, string):
        self.file.write(f"{self.__tabs()}{string}")

    # Add empty lines
    def __lines(self, n):
        self.file.write("\n" * n)

    # Add the current indent
    def __tabs(self):
        return "\t" * self.indent

    # Write C++ typedef
    def __typedef(self, name, type):
        self.__put(f"using {name} = {type};\n")

    # Write C++ template parameter list
    def __template_list(self, params=[]):
        return "" if not params else f"template <{','.join(params)}>"

    # Get the the class name from the full type (which has namespace, template params etc.)
    def __name_from_specifier(self, specifier):
        # Strip template parameter list
        tokens = specifier.split("<")
        tp_str = tokens[0]

        tokens = []

        # Strip namespace
        tokens = tp_str.split(":")

        return tokens[-1]

    # Beginning of the header
    def __preamble(self, det_name, shapes, materials, accel_structs):
        # Write the copyright statement
        year = datetime.now().year
        copy_right = f"\
/** Detray library, part of the ACTS project (R&D line)\n\
 *\n\
 * (c) {year} CERN for the benefit of the ACTS project\n\
 *\n\
 * Mozilla Public License Version 2.0\n\
 */"
        self.__put(copy_right)
        self.__lines(2)
        self.__add_header_includes(shapes, materials, accel_structs)
        self.__lines(2)
        self.__put("namespace detray {")
        self.__lines(2)
        # Write the C++ metatdata struct definition
        metadata_template_params = ["concepts::algeba algebra_type"]
        struct_def = f"{self.__template_list(metadata_template_params)}\nstruct {det_name}_metadata {{"
        self.__put(struct_def)
        self.indent = self.indent + 1
        self.__lines(2)

    # Close the header file
    def __finish(self):
        self.__lines(2)
        self.indent = self.indent - 1
        self.__put(f"}};\n\n}} // namespace detray")
        self.file.close()

    # Write a full mask type
    def __declare_mask(self, shape, algebra, link, shape_params={}):
        type_specifier = f"{self.__name_from_specifier(shape.specifier)}_t"

        # Distinguish line types
        if type_specifier == "line_square_t" or type_specifier == "line_circular_t":
            type_specifier = (
                "drift_cell_t" if type_specifier == "line_square_t" else "straw_tube_t"
            )

        self.shape_specifiers.append(type_specifier)

        # Prepate extra template parameters
        params = [str(v).lower() for v in shape_params.values()]
        template_params = "" if not shape_params else f"<{','.join(params)}>"

        self.__put(
            f"using {type_specifier} = mask<{shape.specifier}{template_params}, {algebra}, {link}>;\n"
        )

    # Write a full material type
    def __declare_material(self, mat, algebra):
        type_specifier = f"{self.__name_from_specifier(mat.specifier)}_t"

        if mat is Material.SLAB:
            self.__put(f"using {type_specifier} = slab<scalar_t>;\n")
        elif mat is Material.ROD:
            self.__put(f"using {type_specifier} = rod<scalar_t>;\n")
        elif mat is Material.VOLUME:
            self.__put(f"using volume_material_t = material<scalar_t>;\n")
        else:
            shape_specifier = mat.param["shape"].specifier
            shape_type = self.__name_from_specifier(shape_specifier)
            type_specifier = f"{shape_type}_map_t"
            template_list = f"template <typename container_t>\n"
            self.__put(
                f"{template_list}{self.__tabs()}using {type_specifier} = material_map<{algebra}, {shape_specifier}, container_t>;\n"
            )

        self.material_specifiers.append(type_specifier)

    # Write a full acceleration structure type
    def __declare_accel(self, acc, algebra):
        type_specifier = f"{self.__name_from_specifier(acc.specifier)}_t"

        if type_specifier.endswith("grid_t"):
            shape_specifier = acc.param["shape"].specifier
            shape_type = self.__name_from_specifier(shape_specifier)
            type_specifier = f"{shape_type}_grid_t"
            template_list = f"template <typename container_t>\n"
            self.__put(
                f"{template_list}{self.__tabs()}using {type_specifier} = surface_grid_t<axes<{shape_specifier}>, surface_type, container_t>;\n"
            )

        self.accel_specifiers.append(type_specifier)

    # Generate the shape IDs
    def __make_type_enum(self, specifier, items, base_type, extra_items=[]):
        self.__put(f"enum class {specifier} : {base_type} {{\n")

        self.indent = self.indent + 1
        for i, v in enumerate(items + extra_items):
            # Remove the trailing '_t' suffix
            specifier = f"e_{v[:-2]}" if v.endswith("_t") else f"e_{v}"
            self.__put(f"{specifier} = {i}u,\n")

        self.indent = self.indent - 1

        self.__put(f"}};")

    # Identify the required dependencies and add the includes to the header
    def __add_header_includes(self, shapes, materials, accel_structs):
        header_str = '#pragma once\n\n// Project include(s)\n\
#include "detray/core/detail/multi_store.hpp"\n\
#include "detray/core/detail/single_store.hpp"\n\
#include "detray/definitions/detail/algebra.hpp"\n\
#include "detray/definitions/detail/containers.hpp"\n\
#include "detray/definitions/detail/indexing.hpp"\n\
#include "detray/geometry/detail/surface_descriptor.hpp"\n\
#include "detray/geometry/mask.hpp\n'
        self.__put(header_str)

        # Shape class names
        shape_names = set(self.__name_from_specifier(s.specifier) for s in shapes)

        # Material class names
        mat_names = set(self.__name_from_specifier(m.specifier) for m in materials)
        # Add shapes for material maps that have not been added, yet
        if "material_map" in mat_names:
            shape_names.update(
                (
                    self.__name_from_specifier(m.param["shape"].specifier)
                    if m.param
                    else None
                )
                for m in materials
            )

        # Acceleration structure class names
        accel_names = set(
            self.__name_from_specifier(a.specifier) for a in accel_structs
        )
        if "grid" in accel_names:
            accel_names.remove("grid")
            accel_names.add("surface_grid")
            # Add shapes for surface grids that have not been added, yet
            shape_names.update(
                (
                    self.__name_from_specifier(a.param["shape"].specifier)
                    if a.param
                    else None
                )
                for a in accel_structs
            )

        # Add the corresponding includes
        shape_names.remove(None)
        [self.__put(f'#include "detray/geometry/shapes/{n}.hpp\n') for n in shape_names]

        [self.__put(f'#include "detray/materials/{n}.hpp\n') for n in mat_names]

        [
            self.__put(f'#include "detray/navigation/accelerators/{n}.hpp\n')
            for n in accel_names
        ]

    # Generate the type declaration for the multi_stores
    def __declare_multi_store(
        self, specifier, id_name, types, context="empty_context", is_regular=True
    ):
        tenplate_params = (
            "template<typename...> class vector_t = dvector"
            if is_regular
            else "typename container_t = host_container_types"
        )
        self.__put(f"{self.__template_list([tenplate_params])}\n")
        self.__put(f"using {specifier} =\n")

        if is_regular:
            self.__put(
                f"\tregular_multi_store<{id_name}, {context}, dtuple, vector_t, {','.join(types)}>;"
            )
        else:
            type_list = []
            if "brute_force_finder_t" in types:
                types.remove("brute_force_finder_t")
                type_list.append("brute_force_collection<surface_type, container_t>")

            type_list = type_list + [
                (
                    f"grid_collection<{t}<container_t>>"
                    if (t.endswith("map_t") or t.endswith("grid_t"))
                    else f"typename container_t::template vector_type<{t}>"
                )
                for t in types
            ]

            self.__put(
                f"\tmulti_store<{id_name}, {context}, dtuple, {','.join(type_list)}>"
            )
