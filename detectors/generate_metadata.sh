# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Print help message
print_help()
{
   echo "Syntax: generate_metadata.sh [-m|b|o|v|vv|h]"
   echo "Options:"
   echo "  m     Full path and filename of the python metadata script"
   echo "        in the \"detray/detectors/python\" directory."
   echo "  b     Where to find the detray python setup script."
   echo "  o     Output directory."
   echo "  v     Verbose logging."
   echo "  vv    Debug level logging."
   echo "  h     Print this help."
   echo
}

# Verbosity level of the generator script
log_lvl=0

# Parse options
while getopts ":hvb:vv:m:o:" opt; do
    case $opt in
        h)
            echo "Generate custom detray Detector Metadata"
            echo
            print_help
            exit
        ;;
        v)
            ((log_lvl++))
        ;;
        b)
            build_dir=$OPTARG
        ;;
        m)
            metadata_generator=$OPTARG
        ;;
        o)
            out_dir=$OPTARG
        ;;
        \?)
         echo
         echo "ERROR: Invalid option $opt! Usage:"
         echo
         print_help
         exit;;
    esac
done

python_command="python $metadata_generator"

if [[ -z "$metadata_generator" ]];then
    echo
    echo "ERROR: No detector metadata scrip supplied! Usage:"
    echo
    print_help
    exit
fi

# Configure verbosity
if [[ $log_lvl == 1 ]]; then
    python_command="$python_command -v"
elif [[ $log_lvl == 2 ]]; then
    python_command="$python_command -vv"
fi

# Add the option for the custom output location
if [[ ! -z "$out_dir" ]]; then
    python_command="$python_command -o $out_dir"
fi

# Set uo the detray python package
source $build_dir/python/setup.sh

# Run the generation of the requested metadata
eval $python_command
