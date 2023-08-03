#!/bin/bash

set -e

# Obtain the path to the script
SCRIPT_PATH="$(dirname "$(realpath "$0")")"
# Obtain the root directory path by removing /src from the script path
ROOT_PATH="${SCRIPT_PATH%/src}"

PYTHON_SCRIPT="$ROOT_PATH/src/visualization/difference_visualization.py"
FILES_PATH="$ROOT_PATH/../feelppdb/exampleShadingMask/np_1/shadingMasks"
DESTINATION="$ROOT_PATH/results/Differences"

python3 $PYTHON_SCRIPT --dir_path $FILES_PATH --destination $DESTINATION