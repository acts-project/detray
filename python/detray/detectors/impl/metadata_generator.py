# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project includes
from .type_helpers import link
from .definitions import (
    Type,
    Algebra,
    Shape,
    Material,
    Accelerator,
)

# Python includes
from collections import Counter
from datetime import datetime
import itertools
import logging
import numbers
import os
from typing import Optional


""" Class that represents the c++ metadata struct with all of its types """


class metadata:

    def __init__(
        self,
        detector_name,
    ):
        self.logger = logging.getLogger(__name__)
        self.det_name = detector_name
        # Custom linear algebra implementation
        self.algebra = Algebra.ANY
        # Only relevant, if algebra type is not 'ANY'
        self.precision = None
        # Base type for type ID enums
        self.id_base = Type.UINT_LEAST_8
        # Volume link in the masks (resolve geometry connectivity in navigation)
        self.nav_link = link(link_type="single", data_type=Type.UINT_LEAST_16)
        # Surface link to its mask(s) (type ID, index range)
        self.mask_link = link(link_type="range", data_type=Type.UINT_32)
        # Surface link to its material (type ID, index)
        self.material_link = link(link_type="single", data_type=Type.UINT_32)
        # Surface types: passive, portal, sensitive etc.
        self.surface_types = []
        # Surface shapes: rectangle, ring etc.
        self.shapes = {}
        # Material types: slabs, material maps etc.
        self.materials = {}
        # Set brute force finder as default for portals/passives and volumes
        self.default_surface_accel = Accelerator.BRUTE_FORCE
        self.default_volume_accel = Accelerator.BRUTE_FORCE
        # Surface and volume acceleration structures
        self.acceleration_structs = {}

        self.logger.info(f'Detector: "{self.det_name}"')

    # Choose an algebra-plugin
    def set_algebra_plugin(self, plugin: Algebra, precision: Optional[Type] = None):
        self.algebra = plugin

        if self.algebra is not Algebra.ANY:
            if precision is not Type.SINGLE or precision is not Type.DOUBLE:
                self.logger.warning(
                    f'Incorrect precision "{precision}" for algebra type. Using "float" instead'
                )
                self.precision = Type.SINGLE
            else:
                self.precision = precision
            self.logger.info(f"Algebra plugin: {self.algebra}<{self.precision}>")
        elif precision is not None:
            self.logger.warning(
                f'Precision "{precision}" will be ignored for generic algebra type'
            )

    # Register a new geometric shape for a portal surface at a given position
    # in the mask store (type_id)
    def add_portal(self, shape: Shape, type_id: int = -1):
        self.add_shape(shape, type_id)
        if "portal" not in self.surface_types:
            self.surface_types.append("portal")

    # Register a new geometric shape for a sensitive surface at a given position
    # in the mask store (type_id)
    def add_sensitive(self, shape: Shape, type_id: int = -1):
        self.add_shape(shape, type_id)
        if "sensitive" not in self.surface_types:
            self.surface_types.append("sensitive")

    # Register a new geometric shape for a passive surface at a given position
    # in the mask store (type_id)
    def add_passive(self, shape: Shape, type_id: int = -1):
        self.add_shape(shape, type_id)
        if "passive" not in self.surface_types:
            self.surface_types.append("passive")

    # Register a new geometric shape for the detector at a given position
    # in the mask store (type_id)
    def add_shape(self, shape: Shape, type_id: int):
        if shape not in itertools.chain(*self.shapes.values()):
            self.logger.debug(f'--> surface shape "{shape.specifier}"')

            # Automatically enumerate the shape if no type ID is given
            i = type_id if type_id >= 0 else len(self.shapes)
            self.shapes.setdefault(i, []).append(shape)

    # Register a new material type
    def add_material(self, mat: Material, type_id: int = -1):
        if mat not in itertools.chain(*self.materials.values()):
            if mat is Material.SLAB or mat is Material.ROD or mat is Material.RAW:
                self.logger.debug(f'--> material type "{mat.specifier}"')
            else:
                shape_secifier = mat.param["shape"].specifier
                self.logger.debug(
                    f'--> material type "{mat.specifier}<{shape_secifier}>"'
                )

            # Automatically enumerate the material if no type ID is given
            i = type_id if type_id >= 0 else len(self.materials)
            self.materials.setdefault(i, []).append(mat)

    # Register a new acceleration structure for a geometric object type
    def add_accel_struct(
        self,
        accel: Accelerator,
        obj_type: str = "sensitive",
        type_id: int = -1,
        value_type: str = "surface",
        is_default: bool = False,
    ):
        print(f"{type_id}, {value_type}")
        # Make sure volume acceleration structures have indices as values
        value_type = "index" if "volume" in obj_type else value_type
        chosen = False
        if obj_type not in self.acceleration_structs.keys():
            self.acceleration_structs[obj_type] = [(accel, type_id, value_type)]
            chosen = True
        elif (accel, type_id, value_type) not in self.acceleration_structs[obj_type]:
            chosen = True
            self.acceleration_structs[obj_type].append((accel, type_id, value_type))

        if chosen:
            shape_secifier = (
                "" if "grid" not in accel.specifier else accel.param["shape"].specifier
            )
            accel_type = (
                f"{accel.specifier}"
                if "grid" not in accel.specifier
                else f"{accel.specifier}<{shape_secifier}>"
            )

            self.logger.debug(f"--> accel. struct ({obj_type}): {accel_type}")

        if is_default:
            self.set_default_accel_struct(accel, obj_type, value_type)

    # Mark an acceleration struct as default for the given object type
    def set_default_accel_struct(
        self,
        accel: Accelerator,
        obj_type: str,
        type_id: int = -1,
        value_type: str = "surface",
    ):
        # Make sure volume acceleration structures have indices as values
        value_type = "index" if "volume" in obj_type else value_type

        shape_secifier = (
            "" if "grid" not in accel.specifier else accel.param["shape"].specifier
        )
        accel_type = (
            f"{accel.specifier}"
            if "grid" not in accel.specifier
            else f"{accel.specifier}<{shape_secifier}>"
        )

        if obj_type == "portal":
            self.logger.debug(
                f"--> setting default surface accel. struct: {accel_type}"
            )
            self.default_surface_accel = accel
        elif obj_type == "volume":
            self.logger.debug(f"--> setting default volume accel. struct: {accel_type}")
            self.default_volume_accel = accel
        else:
            self.logger.warning(
                f'Cannot set default acceleration structure for geometry object type "{obj_type}"'
            )

        # Make sure the requested default exists
        if obj_type not in self.acceleration_structs.keys() or (
            accel,
            value_type,
        ) not in [(a, v) for a, i, v in self.acceleration_structs[obj_type]]:
            # print(obj_type)
            # print(self.acceleration_structs)
            self.logger.warning(
                f"Requested default acceleration structure ({obj_type}, {accel_type}) not defined in metadata: Adding it now..."
            )
            self.add_accel_struct(accel, obj_type, type_id, value_type)


""" Class that accumulates detector type data and writes a metadat header """


class metadata_generator:

    def __init__(
        self,
        md: metadata,
    ):
        # Internal state during header generation
        self.file = None
        self.logger = logging.getLogger(__name__)
        self.indent = 0

        # Save the type specifiers
        self.shape_specifiers = {}
        self.material_specifiers = {}
        self.accel_specifiers = {}

        # Common header section
        root_dir = '#include "detray'
        self.__common_includes = f'#pragma once\n\n// Project include(s)\n\
{root_dir}/core/detail/multi_store.hpp"\n\
{root_dir}/core/detail/single_store.hpp"\n\
{root_dir}/definitions/algebra.hpp"\n\
{root_dir}/definitions/containers.hpp"\n\
{root_dir}/definitions/indexing.hpp"\n\
{root_dir}/geometry/mask.hpp"\n\
{root_dir}/geometry/surface_descriptor.hpp"\n'

        # Dump the header
        self.__generate(md)

    # Write the given metadata to file
    def __generate(
        self, md: metadata, src_dir="../detectors/include/detray/detectors/"
    ):
        det_name = md.det_name.replace("_detector", "")
        filename = f"{os.path.abspath(src_dir)}/{det_name}_metadata.hpp"
        self.logger.debug(f'Open "{filename}"')
        self.file = open(filename, "w+")

        self.logger.info("Generating metadata...")

        # Beginning of the header (copyright, includes etc.)
        self.__preamble(md)

        # Basic typedefs
        self.__local_typedefs(md)

        # Transform store
        self.logger.info(" -> Transforms")
        self.__declare_transform_store()

        # Mask store
        self.logger.info(" -> Masks")
        self.__declare_mask_store(md)

        # Mask store
        self.logger.info(" -> Material")
        self.__declare_material_store(md)

        # Surface type to be used in surface accelerator types
        # (The volume type is automatically assembled in the detector class)
        self.__lines(2)
        self.__declare_surface_descriptor(md)

        # Acceleration structure store
        self.logger.info(" -> Acceleration Structures")
        self.__declare_accel_store(md)

        # Definition of geometric object types in a detector volume
        self.__declare_geometry_objects(md)

        # Finish (write header file to disk)
        self.__finish()

        self.logger.info("Done!")

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
    def __preamble(self, md: metadata):
        copy_right = f"\
/** Detray library, part of the ACTS project (R&D line)\n\
 *\n\
 * (c) {datetime.now().year} CERN for the benefit of the ACTS project\n\
 *\n\
 * Mozilla Public License Version 2.0\n\
 */"
        self.__put(copy_right)
        self.__lines(2)
        self.__add_header_includes(md.shapes, md.materials, md.acceleration_structs)
        self.__lines(2)
        self.__put("namespace detray {")
        self.__lines(2)
        # Write the C++ metatdata struct definition
        template_params = (
            f"{self.__template_list([str(md.algebra)])}\n"
            if md.algebra == Algebra.ANY
            else ""
        )
        struct_def = f"{template_params}struct {md.det_name}_metadata {{"
        self.__put(struct_def)
        self.indent = self.indent + 1
        self.__lines(2)

    # Identify the required dependencies and add the includes to the header
    def __add_header_includes(self, shapes, materials, accel_structs):
        # Add the common includes
        self.__put(self.__common_includes)

        # Set of shape class names
        shape_names = {
            self.__name_from_specifier(s.specifier)
            for s in itertools.chain(*shapes.values())
        }

        # Correct the header name for the line surfaces
        add_line = False
        if "line_circular" in shape_names:
            shape_names.remove("line_circular")
            add_line = True
        if "line_square" in shape_names:
            shape_names.remove("line_square")
            add_line = True

        if add_line and "line" not in shape_names:
            shape_names.add("line")

        # Set of material class names
        mat_names = {
            self.__name_from_specifier(m.specifier)
            for m in itertools.chain(*materials.values())
        }
        # Add shapes for material maps that have not been added, yet
        if "material_map" in mat_names:
            shape_names.update(
                (
                    self.__name_from_specifier(m.param["shape"].specifier)
                    if m.param
                    else None
                )
                for m in itertools.chain(*materials.values())
            )

        # Set of acceleration structure class names
        accel_names = {
            self.__name_from_specifier(a.specifier)
            for accels in accel_structs.values()
            for a, i, v in accels
        }
        if "spatial_grid" in accel_names:
            # Add shapes for spatial grids that have not been added, yet
            shape_names.update(
                (
                    self.__name_from_specifier(a.param["shape"].specifier)
                    if a.param
                    else None
                )
                for accels in accel_structs.values()
                for a, i, v in accels
            )

        # Add the corresponding includes
        root_dir = '#include "detray'
        shape_names.discard(None)
        [self.__put(f'{root_dir}/geometry/shapes/{n}.hpp"\n') for n in shape_names]

        [self.__put(f'{root_dir}/material/{n}.hpp"\n') for n in mat_names]

        [
            self.__put(f'{root_dir}/navigation/accelerators/{n}.hpp"\n')
            for n in accel_names
        ]

    # Add some mandatory local typedefs of the metadata class
    def __local_typedefs(self, md: metadata):
        self.__typedef(
            "algebra_type",
            (
                "algebra_t"
                if md.algebra == Algebra.ANY
                else f"{md.algebra}<{md.precision}>"
            ),
        )
        self.__typedef("scalar_t", "dscalar<algebra_type>")
        self.__lines(1)
        self.__typedef("nav_link", md.nav_link.data_type)

    # Declare the transform store type
    def __declare_transform_store(self):
        self.__lines(1)
        self.__put(
            f"\
template <template <typename...> class vector_t = dvector>\n\
{self.__tabs()}using transform_store =\n\
{self.__tabs()}    single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;"
        )

    # Generate the IDs enums for the multi store
    def __declare_type_enum(self, specifier, items, base_type, extra_items={}):
        self.__put(f"enum class {specifier} : {base_type} {{\n")

        print(items)
        print(extra_items)

        self.indent = self.indent + 1
        for i, values in itertools.chain(items.items(), extra_items.items()):
            for v in values:
                # Special value was passed
                item = v[0] if isinstance(v, list) else v
                tmp_value = v[1] if isinstance(v, list) else i
                value = f"{i}u" if isinstance(tmp_value, numbers.Number) else tmp_value

                print(value)
                # Remove the trailing '_t' suffix
                sub_specifier = f"e_{item[:-2]}" if item.endswith("_t") else f"e_{item}"
                self.__put(f"{sub_specifier} = {value},\n")

        self.indent = self.indent - 1

        self.__put("};")

    # Add the streaming operator for an enum
    def __define_enum_stream_op(self, specifier, items, extra_items={}):
        self.__put(
            f"DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os, {specifier} id) {{\n"
        )

        self.indent = self.indent + 1

        self.__put("switch (id) {\n")

        self.indent = self.indent + 1
        for i, values in itertools.chain(items.items(), extra_items.items()):
            for v in values:
                # Special value was passed
                item = v[0] if isinstance(v, list) else v
                value = v[1] if isinstance(v, list) else f"{i}u"
                # Remove the trailing '_t' suffix
                sub_specifier = f"e_{item[:-2]}" if item.endswith("_t") else f"e_{item}"

                self.__put(f"case {specifier}::{sub_specifier}:\n")
                self.indent = self.indent + 1
                self.__put(
                    f'os << "{value if isinstance(v, list) else sub_specifier}";\n'
                )
                self.__put("break;\n")
                self.indent = self.indent - 1

        # Add the default case
        self.__put("default:\n")
        self.indent = self.indent + 1
        self.__put('os << "invalid";\n')
        self.indent = self.indent - 1

        self.indent = self.indent - 1

        self.__put("}\n")
        self.__put("return os;\n")

        self.indent = self.indent - 1

        self.__put("};")

    # Generate the type declaration for the multi_stores
    def __declare_multi_store(
        self, specifier, id_name, types, context="empty_context", is_regular=True
    ):
        template_params = (
            "template<typename...> class vector_t = dvector"
            if is_regular
            else "typename container_t = host_container_types"
        )
        self.__put(f"{self.__template_list([template_params])}\n")
        self.__put(f"using {specifier} =\n")

        if is_regular:
            self.__put(
                f"\tregular_multi_store<{id_name}, {context}, dtuple, vector_t, {','.join(itertools.chain(*types.values()))}>;"
            )
        else:
            # The brute force searcher needs to be at the beginning of the store
            type_list = []
            flattened_types = [t for t in itertools.chain(*types.values())]

            if "surface_brute_force_t" in flattened_types:
                flattened_types = [
                    t for t in flattened_types if t != "surface_brute_force_t"
                ]
                type_list.append("brute_force_collection<surface_type, container_t>")

            # Have the volume acceleration structure at the end of the store
            volume_brute_force = ""
            if "volume_brute_force_t" in flattened_types:
                flattened_types = [
                    t for t in flattened_types if t != "volume_brute_force_t"
                ]
                volume_brute_force = "brute_force_collection<dindex, container_t>"

            type_list = type_list + [
                (
                    f"grid_collection<{t}<container_t>>"
                    if (t.endswith("map_t") or t.endswith("grid_t"))
                    else f"typename container_t::template vector_type<{t}>"
                )
                for t in flattened_types
            ]

            if volume_brute_force:
                type_list.append(volume_brute_force)

            self.__put(
                f"\tmulti_store<{id_name}, {context}, dtuple, {','.join(type_list)}>;"
            )

    # Declare the surface descriptor type with all links
    def __declare_surface_descriptor(self, md: metadata):
        self.__put("using transform_link = typename transform_store<>::single_link;\n")

        link_type = (
            "single_link" if md.mask_link.link_type == "single" else "range_link"
        )
        mask_link = f"typename mask_store<>::{link_type}"
        self.__put(f"using mask_link = {mask_link};\n")

        link_type = (
            "single_link" if md.material_link.link_type == "single" else "range_link"
        )
        material_link = f"typename material_store<>::{link_type}"
        self.__put(f"using material_link = {material_link};\n")

        self.__put(
            "using surface_type = surface_descriptor<mask_link, material_link, transform_link, nav_link>;"
        )

    # Write a full mask type
    def __declare_mask(self, shape, type_id, algebra, link, shape_params={}):
        type_specifier = f"{self.__name_from_specifier(shape.specifier)}_t"

        # Distinguish line types
        if type_specifier == "line_square_t" or type_specifier == "line_circular_t":
            type_specifier = (
                "drift_cell_t" if type_specifier == "line_square_t" else "straw_tube_t"
            )

        self.shape_specifiers.setdefault(type_id, []).append(type_specifier)

        # Prepate extra template parameters
        params = [str(v).lower() for v in shape_params.values()]
        template_params = "" if not shape_params else f"<{','.join(params)}>"

        self.__put(
            f"using {type_specifier} = mask<{shape.specifier}{template_params}, {algebra}, {link}>;\n"
        )

    # Declare the mask types and mask store for the detector
    def __declare_mask_store(self, md: metadata):
        # Mask types
        assert md.shapes, "Define at least one geometric shape"

        self.__lines(2)
        for type_id, shapes in md.shapes.items():
            for shape in shapes:
                self.__declare_mask(
                    shape, type_id, "algebra_type", "nav_link", shape.param
                )

        self.__lines(1)
        self.__declare_type_enum("mask_id", self.shape_specifiers, md.id_base)
        self.__lines(2)
        self.__define_enum_stream_op("mask_id", self.shape_specifiers)
        self.__lines(2)

        # Mask Store
        self.__declare_multi_store("mask_store", "mask_id", self.shape_specifiers)

    # Write a full material type
    def __declare_material(self, mat, type_id, algebra):
        type_specifier = f"{self.__name_from_specifier(mat.specifier)}_t"

        if mat is Material.SLAB:
            self.__put(f"using {type_specifier} = material_slab<scalar_t>;\n")
        elif mat is Material.ROD:
            self.__put(f"using {type_specifier} = material_rod<scalar_t>;\n")
        elif mat is Material.RAW:
            type_specifier = "raw_material_t"
            self.__put(f"using {type_specifier} = material<scalar_t>;\n")
        else:
            shape_specifier = mat.param["shape"].specifier
            shape_type = self.__name_from_specifier(shape_specifier)
            type_specifier = f"{shape_type}_map_t"
            template_list = "template <typename container_t>\n"
            self.__put(
                f"{template_list}{self.__tabs()}using {type_specifier} = material_map<{algebra}, {shape_specifier}, container_t>;\n"
            )

        self.material_specifiers.setdefault(type_id, []).append(type_specifier)

    # Declare the material types and material store for the detector
    def __declare_material_store(self, md: metadata):
        # Material types
        if md.materials:
            self.__lines(2)
            for type_id, materials in md.materials.items():
                for mat in materials:
                    self.__declare_material(mat, type_id, "algebra_type")

            self.__lines(1)
            # Add an option for 'no material'
            self.__declare_type_enum(
                "material_id",
                self.material_specifiers,
                md.id_base,
                {len(md.materials): ["none"]},
            )
            self.__lines(2)
            self.__define_enum_stream_op(
                "material_id", self.material_specifiers, {len(md.materials): ["none"]}
            )
        self.__lines(2)

        # Material Store
        self.__declare_multi_store(
            "material_store", "material_id", self.material_specifiers, is_regular=False
        )

    # Write a full acceleration structure type
    def __declare_accel(self, obj_type, acc, type_id: int, value_type: str):
        type_specifier = f"{self.__name_from_specifier(acc.specifier)}_t"

        if type_specifier.endswith("grid_t"):
            # First spatial grid declaration?
            if not any(
                sc.endswith("grid_t")
                for sc in itertools.chain(*self.accel_specifiers.values())
            ):
                template_list = "template <typename axes_t, typename bin_entry_t, typename container_t>\n"
                self.__put(
                    f"{template_list}{self.__tabs()}using spatial_grid_t = spatial_grid<algebra_type, axes_t,bins::dynamic_array<bin_entry_t>, simple_serializer, container_t, false>;"
                )
                self.__lines(2)

            shape_specifier = acc.param["shape"].specifier
            type_specifier = f"{self.__name_from_specifier(shape_specifier)}_grid_t"
            entry_type = "surface_type" if value_type == "surface" else "dindex"

            template_list = "template <typename container_t>\n"
            self.__put(
                f"{template_list}{self.__tabs()}using {obj_type}_{type_specifier} = spatial_grid_t<axes<{shape_specifier}>, {entry_type}, container_t>;\n"
            )

        self.accel_specifiers.setdefault(type_id, []).append(
            f"{obj_type}_{type_specifier}"
        )

    # Declare the accleration structure types and accel store for the detector
    def __declare_accel_store(self, md: metadata):
        # Acceleration Structures
        assert (
            len(md.acceleration_structs) > 2
        ), "Define at least one default surface(portal) and one default volume acceleration structure"

        assert (
            "portal" in md.acceleration_structs.keys()
        ), "Define at least one portal acceleration structure"

        assert (
            "volume" in md.acceleration_structs.keys()
        ), "Define at least one volume acceleration structure"

        if not md.acceleration_structs:
            return

        self.__lines(2)
        # Make sure an acceleration struct is declared only once,
        # even if it is used for multiple surface types
        unique_accel = []
        for geo_obj, accels in md.acceleration_structs.items():
            obj_type = "volume" if "volume" in geo_obj else "surface"
            for acc, type_id, value_type in accels:
                if (obj_type, acc, type_id) not in unique_accel:
                    self.__declare_accel(obj_type, acc, type_id, value_type)
                    unique_accel.append((obj_type, acc, type_id, value_type))

        # Add options for the default acceleration structs
        self.__lines(1)

        # Set enum item for default surface acceleration structure
        surface_default = ""
        if "grid" in md.default_surface_accel.specifier:
            shape_specifier = md.default_surface_accel.param["shape"].specifier
            surface_default = (
                f"e_surface_{self.__name_from_specifier(shape_specifier)}_grid"
            )
        else:
            surface_default = f"e_surface_{self.__name_from_specifier(md.default_surface_accel.specifier)}"

        # Set enum item for default volume acceleration structure
        volume_default = ""
        if "grid" in md.default_volume_accel.specifier:
            shape_specifier = md.default_volume_accel.param["shape"].specifier
            volume_default = (
                f"e_volume_{self.__name_from_specifier(shape_specifier)}_grid"
            )
        else:
            volume_default = f"e_volume_{self.__name_from_specifier(md.default_volume_accel.specifier)}"

        extra_items = {
            surface_default: ["surface_default"],
            volume_default: ["volume_default"],
        }
        self.__declare_type_enum(
            "accel_id",
            self.accel_specifiers,
            md.id_base,
            extra_items=extra_items,
        )

        self.__lines(2)
        self.__define_enum_stream_op(
            "accel_id",
            self.accel_specifiers,
        )

        # Accelerator Store
        self.__lines(2)
        self.__declare_multi_store(
            "accelerator_store", "accel_id", self.accel_specifiers, is_regular=False
        )

    # Declare the geometry object types (e.g. passive, portal, sensitive)
    # that are distinguishable to the volume: use enum values to link a
    # geometry object type category to the corresponding accel. struct
    def __declare_geometry_objects(self, md: metadata):

        # Helper to determine if an accel. struct is one of the defaults
        def is_default_accel(value_type, accel):
            is_surface_default = (
                value_type == "surface" and accel is md.default_surface_accel.specifier
            )
            is_volume_default = (
                value_type == "index" and accel is md.default_volume_accel.specifier
            )

            return is_surface_default or is_volume_default

        # Enumerate the different geometric objects and assign them to the accel
        # link of the corresponding acceleration structure, e.g.
        # geo_id:  value_type:        accel_link [accel_id, accel_idx]
        # ---------------------------------------------------------------------
        # 0:       portal/passive:  [brute_force, 0] <- two in one accel. struct
        # 1:       sensitive:       [disk_grid, 14]  <- distinct accel. struct
        # 2:       volume:          [brute_force, 2] <- distinct accel. struct

        # Find which geometric objects go together into which accel. struct

        print(f"Input:\n{md.acceleration_structs}\n")

        accel_to_types = {}
        for t, accels in md.acceleration_structs.items():
            for a, i, v in accels:
                specifier = f"{a.specifier}"

                if (specifier, i, v) not in accel_to_types:
                    accel_to_types[specifier, i, v] = set()

                accel_to_types[specifier, i, v].add(t)

        print(f"Accel to types:\n{accel_to_types}\n")

        # If specific type IDs are required, merge the categories accordingly
        type_id_dict = {}
        for acc, type_id, value_type in accel_to_types.keys():
            if type_id > -1:
                type_id_dict.setdefault((type_id, value_type), []).append(acc)

        print(f"type ID dict:\n{type_id_dict}\n")

        # Find a representative for each category and merge the obj type sets
        accel_to_types_merged = {}
        for acc, type_id, value_type in accel_to_types.keys():
            representative = acc
            if (
                not is_default_accel(value_type, acc)
                and (type_id, value_type) in type_id_dict.keys()
                and len(type_id_dict[type_id, value_type]) > 1
            ):
                tmp_repr = type_id_dict[type_id, value_type][0]
                representative = (
                    type_id_dict[type_id, value_type][1]
                    if is_default_accel(value_type, tmp_repr)
                    else tmp_repr
                )

            accel_to_types_merged.setdefault(
                (representative, value_type, type_id), set()
            ).update(accel_to_types[acc, type_id, value_type])

        print(f"Accel to types (merged):\n{accel_to_types_merged}\n")

        # Check duplicate sets (make the set hashable by converting to tuple)
        counted_dict = Counter([tuple(v) for v in accel_to_types_merged.values()])
        duplicates = {
            key: value
            for key, value in accel_to_types_merged.items()
            if counted_dict[tuple(value)] > 1
        }

        print(f"Counted: {counted_dict}")
        print(f"duplicates: {duplicates}")

        # Get the list of keys for every duplicate
        flipped_duplicates = {}
        for key, value in duplicates.items():
            flipped_duplicates.setdefault(tuple(value), []).append(key)

        print(f"flipped_duplicates: {flipped_duplicates}")

        # Find best representative key for the duplicates and remove the others
        removal_keys = []
        for obj_set, keys in flipped_duplicates.items():
            # Keep the key with the smallest required type ID
            type_id_counter = 10000
            for acc, value_type, type_id in keys:
                if is_default_accel(value_type, acc):
                    continue
                # Remove the non-specific type ID of "-1" in any case
                if type_id == -1 or type_id >= type_id_counter:
                    removal_keys.append((acc, value_type, type_id))
                else:
                    type_id_counter = type_id

        for rk in removal_keys:
            accel_to_types_merged.pop(rk, None)

        print(f"Accel to types (deduplicated):\n{accel_to_types_merged}\n")

        # Sort, in order to find the largest geo object sets
        accel_to_types_sorted = dict(
            sorted(
                accel_to_types_merged.items(),
                key=lambda item: len(item[1]),
                reverse=True,
            )
        )

        print(f"Accel to types (sorted):\n{accel_to_types_sorted}\n")

        # Link slot to geometry object types
        accel_link_types = {}

        # Add a new link slot, if there are previously unknown object types
        def add_link(key, new_objects):
            new_link_types = [
                l
                for l in new_objects
                if l not in itertools.chain(*accel_link_types.values())
            ]

            if new_link_types:
                accel_link_types[key] = new_link_types

        for (accel, value_type, type_id), obj_set in accel_to_types_sorted.items():
            # Only add the default methods where no other accelerator is
            # available (see below)
            if is_default_accel(value_type, accel):
                continue

            add_link(key=type_id, new_objects=obj_set)

        print(f"Final accel distinct link types:\n{accel_link_types}\n")

        # Add the default link for the leftover surface and volume types
        surface_default_list = [
            v
            for k, v in accel_to_types_sorted.items()
            if "surface" in k and md.default_surface_accel.specifier in k
        ]
        add_link(
            key="surface_default",
            new_objects=itertools.chain.from_iterable(surface_default_list),
        )

        volume_default_list = [
            v
            for k, v in accel_to_types_sorted.items()
            if "index" in k and md.default_volume_accel.specifier in k
        ]
        add_link(
            key="volume_default",
            new_objects=itertools.chain.from_iterable(volume_default_list),
        )

        print(f"final accel links categories: {accel_link_types}")

        # Write the id enum for the types of distinct geometric objects
        self.__lines(2)
        self.__put(f"enum geo_objects : {md.id_base} {{\n")

        self.indent = self.indent + 1

        # The surface accelerator that contains the portals has to be in slot 0:
        (portal_key, portal_group) = [
            (k, v) for (k, v) in accel_link_types.items() if "portal" in v
        ][0]

        print(portal_key)
        print(portal_group)

        if not portal_group:
            self.logger.error(
                "No acceleration structure link define for portal surfaces!"
            )
            return

        for gid in portal_group:
            self.__put(f"e_{gid} = 0u,\n")

        # Loop over the other link groups and enumerate them
        i = 1
        for key, link_group in accel_link_types.items():
            if key == portal_key:
                continue
            for gid in link_group:
                self.__put(f"e_{gid} = {i}u,\n")
                i = i + 1

        self.__put(f"e_size = {len(accel_link_types)}u,\n")
        self.__put("e_all = e_size,\n")

        self.indent = self.indent - 1

        self.__put("};")

        # Add streaming operator for geo_objects enum
        # print(list(itertools.chain.from_iterable(accel_link_types.values())))
        # self.__lines(2)
        # self.__define_enum_stream_op(
        #    "geo_objects",
        #    list(itertools.chain.from_iterable(accel_link_types.values())),
        #    extra_items={"size": ["e_size"]},
        # )

        # Add the surface acceleration struct link to the volume descriptor
        self.__lines(2)
        self.__put(
            "using object_link_type = dmulti_index<dtyped_index<accel_id, dindex>, geo_objects::e_size>;"
        )

    # Close the header file
    def __finish(self):
        self.__lines(1)
        self.indent = self.indent - 1
        self.__put("};\n\n} // namespace detray")
        self.file.close()

        self.logger.debug("Wrote file to disk")
