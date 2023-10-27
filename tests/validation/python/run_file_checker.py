# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import argparse
import json
import os
import sys
from jsonschema import validate

from json_schemas.geometry import geometry_schema
from json_schemas.homogeneous_material import homogeneous_material_schema
from json_schemas.surface_grids import surface_grid_schema

def __main__():
#----------------------------------------------------------------arg parsing

    parser = argparse.ArgumentParser(description = "Detray File Validation")

    parser.add_argument("--geometry_file",
                        help=("Input geometry json file."),
                        default = "", type=str)
    parser.add_argument("--homogeneous_material_file",
                        help=("Input homogeneous material json file."),
                        default = "", type=str)
    parser.add_argument("--grid_file",
                        help=("Surface grid json file."),
                        default = "", type=str)

    args = parser.parse_args()

    # Check input json files
    filename_dict = {}

    geo_file = args.geometry_file
    if not geo_file == "":
        if not os.path.isfile(geo_file):
            print(f"Geometry file does not exist! ({geo_file})")
            sys.exit(1)
        else:
            filename_dict[geo_file] = geometry_schema

    mat_file = args.homogeneous_material_file
    if not mat_file == "":
        if not os.path.isfile(mat_file):
            print(f"Homogeneous material file does not exist! ({mat_file})")
            sys.exit(1)
        else:
            filename_dict[mat_file] = homogeneous_material_schema

    grid_file = args.grid_file
    if not grid_file == "":
        if not os.path.isfile(grid_file):
            print(f"Surface grid file does not exist! ({grid_file})")
            sys.exit(1)
        else:
            filename_dict[grid_file] = surface_grid_schema

#------------------------------------------------------------------------run

    for filename, schema in filename_dict.items():
        with open(filename) as file:
            try:
                input_json = json.load(file)
            except json.decoder.JSONDecodeError:
                print(f"Invalid json file: {filename}")
            else:
                validate(instance=input_json, schema=schema)
                print(f"{filename}: OK")

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

#-------------------------------------------------------------------------------
