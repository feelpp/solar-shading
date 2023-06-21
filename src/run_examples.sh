#!/bin/bash

set -e

# Obtain the path to the script
SCRIPT_PATH="$(dirname "$(realpath "$0")")"
# Obtain the root directory path by removing /src from the script path
ROOT_PATH="${SCRIPT_PATH%/src}"

cd "$ROOT_PATH/build/default/src"

./feelpp_ss_example_ShadingMasks_comparison --config-file "$ROOT_PATH/src/cases/exampleShadingMask/exampleShadingMask.cfg"

PYTHON_SCRIPT="$ROOT_PATH/src/visualization/shadingMask_visualization.py"
FILES_PATH="$ROOT_PATH/../feelppdb/exampleShadingMask/np_1/shadingMasks"
DESTINATION="$ROOT_PATH/results/ShadingMasks"

for FILE in "$FILES_PATH"/*.csv; do
    NAME=$(basename "$FILE")
    python3 $PYTHON_SCRIPT --file_path $FILE --dir_path $FILES_PATH --destination $DESTINATION
done
