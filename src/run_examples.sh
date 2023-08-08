#!/bin/bash

set -e

# Obtain the path to the script
SCRIPT_PATH="$(dirname "$(realpath "$0")")"
# Obtain the root directory path by removing /src from the script path
ROOT_PATH="${SCRIPT_PATH%/src}"

cd "$ROOT_PATH/build/default/src"

echo "Select the config file: "
echo "1. exampleShadingMask"
echo "2. 19buildingsStrasbourg"
read -p "Enter the number (1 or 2): " NUMBER1

case "$NUMBER1" in
  1)
    CONFIG_FILE="$ROOT_PATH/src/cases/exampleShadingMask/exampleShadingMask.cfg"
    ;;
  2)
    CONFIG_FILE="$ROOT_PATH/src/cases/19buildingsStrasbourg/19buildingsStrasbourg.cfg"
    ;;
  *)
    echo "Invalid input! Enter 1 or 2."
    exit 1
esac

echo "Select the compilation flags to be used:"
echo "- 1 : -Ofast -march=native"
echo "- 2 : -O3 -march=native"
echo "- 3 : -O2 -march=native"
echo "- 4 : -Ofast -march=native -funroll-loops"
echo "- 5 : -O3 -march=native -funroll-loops"
echo "- 6 : -O2 -march=native -funroll-loops"
echo "- 7 : -Ofast -march=native -funroll-loops -funsafe-math-optimizations"
read -p "Enter the number (1 to 7): " NUMBER

EXECUTABLE="./example_ShadingMasks_comparison${NUMBER}"

echo "Running the command : $EXECUTABLE --config-file $CONFIG_FILE"

$EXECUTABLE --config-file $CONFIG_FILE

PYTHON_SCRIPT="$ROOT_PATH/src/visualization/shadingMask_visualization.py"
FILES_PATH="$ROOT_PATH/../feelppdb/exampleShadingMask/np_1/shadingMasks"
DESTINATION="$ROOT_PATH/results/ShadingMasks"

python3 $PYTHON_SCRIPT --dir_path $FILES_PATH --destination $DESTINATION
