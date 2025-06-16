# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# python includes
import copy
import json

"""Append a surface to the volume and save the new surface index"""


def __append_surface(volume, surface, new_idx, idx_dict):
    # If the 'mask' entry has not been updated, remove it now
    if "mask" in surface:
        surface["masks"] = [surface["mask"]]
        del surface["mask"]

    # Map the original surface index to the new one after merging
    idx_dict[surface["index_in_coll"]] = new_idx
    surface["index_in_coll"] = new_idx

    # Add the updated surface to the volume
    volume["surfaces"].append(surface)

    # Update and return the current surface index
    new_idx += 1
    return new_idx


""" Add a new value to a dictionary of lists. Append if the key exists"""


def __add_to_dict(key, value, dictionary):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]


""" Read geometry data from json and merge surfaces """


def merge_surfaces(logging, in_json):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Geometry: Converted header")

    out_json["data"] = {}
    out_json["data"]["volumes"] = []

    # Map old to new indices per volume
    sf_index_per_volume = {}
    # Total number of surfaces after merging
    n_surfaces = 0

    # Go through each volume
    for volume in in_json["data"]["volumes"]:

        sf_index_per_volume[volume["index"]] = {}
        sf_index_dict = sf_index_per_volume[volume["index"]]

        # Save the input surfaces
        old_surfaces = volume["surfaces"]
        volume["surfaces"] = []

        # Sort the portals by shape for merging
        shape_dict = {}

        # Go through the surfaces
        sf_idx = 0
        for surface in old_surfaces:
            # Do NOT merge sensitive or passive surfaces
            if surface["type"] == 1 or surface["type"] == 2:
                sf_idx = __append_surface(volume, surface, sf_idx, sf_index_dict)
                continue

            # Add portals to shape dict
            shape_id = surface["mask"]["shape"]
            __add_to_dict(shape_id, surface, shape_dict)

        # Create new surfaces of same shape by mergin masks into list

        # Identify the rings that can be merged to a disc, by comparing the z-position (portals do not have rotations or x, y translation)
        ring_dict = {}
        for ring in shape_dict[6]:
            z_pos = ring["transform"]["translation"][2]
            __add_to_dict(z_pos, ring, ring_dict)

        for rings in ring_dict.values():
            # Nothing to merge, just append to output
            if len(rings) == 1:
                sf_idx = __append_surface(volume, rings[0], sf_idx, sf_index_dict)
                continue

            # The merged disc
            disc = copy.deepcopy(rings[0])

            # Remove old mask entry and replace by masks list
            if "mask" not in disc:
                logging.error(f"Geometry: surface already converted {disc}")

            del disc["mask"]
            disc["masks"] = []

            for ring in rings:
                disc["masks"].append(ring["mask"])

            # Copy disc into volume
            sf_idx = __append_surface(volume, disc, sf_idx, sf_index_dict)

            # Add also the other surfaces to the dict
            for i in range(1, len(rings)):
                sf_index_dict[rings[i]["index_in_coll"]] = sf_idx - 1

        # Identify the rings that can be merged to a disc, by comparing the z-position (portals do not have rotations or x, y translation)
        cyl_dict = {}
        for cyl in shape_dict[4]:
            rad = cyl["mask"]["boundaries"][0]
            __add_to_dict(rad, cyl, cyl_dict)

        for sub_cyls in cyl_dict.values():
            # Nothing to merge, just append to output
            if len(sub_cyls) == 1:
                sf_idx = __append_surface(volume, sub_cyls[0], sf_idx, sf_index_dict)
                continue

            # The merged disc
            cyl = copy.deepcopy(sub_cyls[0])

            # Remove old mask entry and replace by masks list
            if "mask" not in cyl:
                logging.error(f"Geometry: surface already converted {cyl}")

            del cyl["mask"]
            cyl["masks"] = []

            # Remove translation
            cyl["transform"]["translation"] = [0.0, 0.0, 0.0]

            for sub_cyl in sub_cyls:
                # Apply the translation of the surface to the mask z-boundaries
                z_shift = sub_cyl["transform"]["translation"][2]

                sub_cyl["mask"]["boundaries"][1] += z_shift
                sub_cyl["mask"]["boundaries"][2] += z_shift

                cyl["masks"].append(sub_cyl["mask"])

            # Copy disc into volume
            sf_idx = __append_surface(volume, cyl, sf_idx, sf_index_dict)

            # Add also the other surfaces to the dict
            for i in range(1, len(sub_cyls)):
                sf_index_dict[sub_cyls[i]["index_in_coll"]] = sf_idx - 1

        out_json["data"]["volumes"].append(volume)
        n_surfaces += sf_idx

    logging.info(f"Geometry: Converted {len(out_json["data"]["volumes"])} volumes")

    # Copy the volume search data structure
    out_json["data"]["volume_grid"] = in_json["data"]["volume_grid"]

    logging.info(f"Geometry: Converted volume acceleration structure")

    # Update the metadata
    out_json["header"]["surface_count"] = n_surfaces
    out_json["header"]["volume_count"] = len(out_json["data"]["volumes"])

    return out_json, sf_index_per_volume


""" Update the surface indices in accelerator grids after surface merging"""


def update_grids(logging, in_json, sf_index_per_volume):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Grids: Converted header")

    out_json["data"] = {}
    out_json["data"]["grids"] = []

    # Go through each volume
    for grids in in_json["data"]["grids"]:
        # Go through all grids in a volume
        for grid in grids["grid_data"]:
            volume_idx = grid["owner_link"]
            sf_idx_map = sf_index_per_volume[volume_idx]

            # Update the surface indices in the grid bins
            for grid_bin in grid["bins"]:
                for i, sf_idx in enumerate(grid_bin["content"]):
                    grid_bin["content"][i] = sf_idx_map[sf_idx]

        out_json["data"]["grids"].append(grids)

    # Update metadata (helps debugging)
    out_json["header"]["grid_count"] = len(out_json["data"]["grids"])

    logging.info(f"Grids: Converted {len(out_json["data"]["grids"])} grids")

    return out_json


""" Update the surface owner indices of the material maps after surface merging"""


def update_material(logging, in_json, sf_index_per_volume):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Material: Converted header")

    if out_json["header"]["common"]["tag"] != "material_maps":
        logging.error("Material: Only material maps conversion is available!")
        return in_json

    out_json["data"] = {}
    out_json["data"]["grids"] = []

    # Go through each volume
    for grids in in_json["data"]["grids"]:
        sf_idx_map = sf_index_per_volume[grids["volume_link"]]

        # Make sure that each surface gets only one material map
        surface_map = []

        # Go through all grids in a volume
        i = 0
        new_grid_data = []
        while i < len(grids["grid_data"]):
            grid = grids["grid_data"][i]

            old_sf_idx = grid["owner_link"]
            new_sf_idx = sf_idx_map[old_sf_idx]

            grid["owner_link"] = new_sf_idx
            if new_sf_idx not in surface_map:
                surface_map.append(new_sf_idx)
                new_grid_data.append(grid)
            else:
                # Make sure, only redundant data is removed
                for new_grid in new_grid_data:
                    if new_grid["owner_link"] == new_sf_idx:
                        assert len(grid["axes"]) == len(new_grid["axes"])
                        nbins0 = grid["axes"][0]["bins"]
                        new_nbins0 = new_grid["axes"][0]["bins"]

                        nbins1 = grid["axes"][1]["bins"]
                        new_nbins1 = new_grid["axes"][1]["bins"]

                        assert nbins0 == new_nbins0
                        assert nbins1 == new_nbins1

                        # Volume material grids
                        nbins2 = 1
                        if len(grid["axes"]) > 2:
                            nbins2 = grid["axes"][2]["bins"]
                            assert nbins2 == new_grid["axes"][2]["bins"]

                        for j in range(0, nbins0 * nbins1 * nbins2):
                            assert (
                                grid["bins"][j]["content"][0]["material"]
                                == new_grid["bins"][j]["content"][0]["material"]
                            )

            i += 1

        grids["grid_data"] = new_grid_data

        assert len(grids["grid_data"]) <= len(surface_map)
        out_json["data"]["grids"].append(grids)

    # Update metadata (helps debugging)
    out_json["header"]["grid_count"] = len(out_json["data"]["grids"])

    logging.info(f"Material: Converted {len(out_json["data"]["grids"])} grids")

    return out_json
