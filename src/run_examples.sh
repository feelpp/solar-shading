#!/bin/bash

set -e

# Obtain the path to the script
SCRIPT_PATH="$(dirname "$(realpath "$0")")"
# Obtain the root directory path by removing /src from the script path
ROOT_PATH="${SCRIPT_PATH%/src}"

cd "$ROOT_PATH/build/default/src"

read -p "Enter the number: " NUMBER

EXECUTABLE="./feelpp_ss_example_ShadingMasks_comparison${NUMBER}"

$EXECUTABLE --config-file "$ROOT_PATH/src/cases/exampleShadingMask/exampleShadingMask.cfg"

PYTHON_SCRIPT="$ROOT_PATH/src/visualization/shadingMask_visualization.py"
FILES_PATH="$ROOT_PATH/../feelppdb/exampleShadingMask/np_1/shadingMasks"
DESTINATION="$ROOT_PATH/results/ShadingMasks"

python3 $PYTHON_SCRIPT --dir_path $FILES_PATH --destination $DESTINATION