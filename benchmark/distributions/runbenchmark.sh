#!/bin/bash

set -e

# Compiler flags array
declare -a flagsArr=("-O1" "-O2" "-O3" "-Ofast" "-O1 -march=native" "-O2 -march=native" "-O3 -march=native" "-Ofast -march=native" "-O1 -march=native -funroll-loops" "-O2 -march=native -funroll-loops" "-O3 -march=native -funroll-loops" "-Ofast -march=native -funroll-loops" "-O1 -march=native -fdisable-tree-cunrolli" "-O2 -march=native -fdisable-tree-cunrolli" "-O3 -march=native -fdisable-tree-cunrolli" "-Ofast -march=native -fdisable-tree-cunrolli")

# Loop through the array of flags
for flag in "${flagsArr[@]}"
do
    # Build the C++ code with the current flag
    g++ -Wno-enum-compare $flag -o uniform uniform.cpp # -I/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/include -L/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

    # Execute the C++ code with the current flag as a command-line argument
    ./uniform $flag

    # Remove the executable
    rm uniform
done

# Execute the Python script for plotting
python3 plotresults.py

# Ask user for confirmation to delete all lines except the first one from the CSV file
read -p "Do you want to delete all lines except the first one from the CSV file? [Y/n] " answer

case $answer in
    [Yy]* )
        sed -i '2,$d' comparison.csv
        echo "All lines except the first one from the CSV file have been deleted."
        ;;
    * )
        echo "Keeping all lines in the CSV file."
        ;;
esac
