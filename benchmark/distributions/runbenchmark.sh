#!/bin/bash

set -e

# Ask for the .cpp file name
echo "Enter the name of the .cpp file (without the .cpp extension):"
read cpp_file_name

# Compiler flags array
declare -a flagsArr=("-O1" "-O2" "-O3" "-Ofast" "-O1 -march=native" "-O2 -march=native" "-O3 -march=native" "-Ofast -march=native" "-O1 -march=native -funroll-loops" "-O2 -march=native -funroll-loops" "-O3 -march=native -funroll-loops" "-Ofast -march=native -funroll-loops" "-O1 -march=native -fdisable-tree-cunrolli" "-O2 -march=native -fdisable-tree-cunrolli" "-O3 -march=native -fdisable-tree-cunrolli" "-Ofast -march=native -fdisable-tree-cunrolli")

# Loop through the array of flags
for flag in "${flagsArr[@]}"
do
    # Build the C++ code with the current flag
    g++ -Wno-enum-compare $flag -o $cpp_file_name $cpp_file_name.cpp

    # Execute the C++ code with the current flag as a command-line argument
    ./$cpp_file_name $flag

    # Remove the executable
    rm $cpp_file_name
done

# Determine the CSV file name based on the CPP file name
csv_file_name="${cpp_file_name}.csv"

# Execute the Python script for plotting, passing the CSV file name
python3 plotresults.py "$csv_file_name"


# Ask user for confirmation to delete all lines except the first one from the CSV file
read -p "Do you want to delete all lines except the first one from the CSV file? [Y/n] " answer

case $answer in
    [Yy]* )
        sed -i '2,$d' $csv_file_name
        echo "All lines except the first one from the CSV file have been deleted."
        ;;
    * )
        echo "Keeping all lines in the CSV file."
        ;;
esac
