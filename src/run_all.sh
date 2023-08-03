#!/bin/bash

set -e

flags=( 
    "-Ofast -march=native"
    "-O3 -march=native"
    "-O2 -march=native"
    "-Ofast -march=native -funroll-loops"
    "-O3 -march=native -funroll-loops"
    "-O2 -march=native -funroll-loops"
    "-Ofast -march=native -funroll-loops -funsafe-math-optimizations"
)

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

declare -A best_times_map
declare -A best_flags_map

for method in "RNG STD" "RNG XOSHIRO" "RNG PCG" "RNG MIX" "RNG EIGEN" "RNG EIGEN CAST" "RNG EIGEN VEC"; do
    best_times_map["$method"]=1000
    best_flags_map["$method"]=""
done

for ((NUMBER=1; NUMBER<${#flags[@]}+1; NUMBER++))
do
    EXECUTABLE="./feelpp_ss_example_ShadingMasks_comparison${NUMBER}"

    echo "Running the command : $EXECUTABLE --config-file $CONFIG_FILE"
    output=$( $EXECUTABLE --config-file $CONFIG_FILE )

    for method in "${!best_times_map[@]}"; do
        execution_time=$(echo "$output" | grep -oP "(?<=Elapsed time for shading mask computation with $method:)[0-9.]+")

        if (( $(echo "$execution_time < ${best_times_map[$method]}" | bc -l) ))
        then
            best_times_map["$method"]=$execution_time
            best_flags_map["$method"]=${flags[$NUMBER]}
        fi
    done
done

# Display the summary
printf "\n\n"
printf "==============================================================================\n"
printf "%-25s %-10s %s\n" "METHOD" "TIME (seconds)" "FLAGS"
printf "==============================================================================\n"
for method in "${!best_times_map[@]}"; do
    printf "%-25s %-15s %s\n" "$method:" "${best_times_map[$method]} s" "'${best_flags_map[$method]}'"
done
printf "==============================================================================\n"