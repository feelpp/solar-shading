#!/bin/bash

# Array of compiler flags
declare -a flagsArr=("-O1" "-O2" "-O3" "-Ofast" "-O1 -march=native" "-O2 -march=native" "-O3 -march=native" "-Ofast -march=native" "-O1 -march=native -funroll-loops" "-O2 -march=native -funroll-loops" "-O3 -march=native -funroll-loops" "-Ofast -march=native -funroll-loops" "-O1 -march=native -fdisable-tree-cunrolli" "-O2 -march=native -fdisable-tree-cunrolli" "-O3 -march=native -fdisable-tree-cunrolli" "-Ofast -march=native -fdisable-tree-cunrolli")

dir="src/reframe"

declare -A bestTimes
declare -A bestFlags

read -p "Enter the type of tests (real/int): " test_type

# Determine the appropriate file extension based on the test type
if [ "$test_type" = "int" ]; then
  file_extension="_int.cpp"
else
  file_extension="_real.cpp"
fi

# For each .cpp file in the directory
for file in $dir/*$file_extension
do
  # Extract the base name (without extension) of the file
  base_name=$(basename "$file" .cpp)
  
  bestTimes[$base_name]=INF
  bestFlags[$base_name]=""

  # For each flag combination
  for flags in "${flagsArr[@]}"
  do
    # Compile the file with the current flags
    g++ -DFROM_SHELL_SCRIPT -Wno-enum-compare $flags -o "${dir}/${base_name}_${flags// /_}" $file

    # If compilation was successful
    if [ $? -eq 0 ]
    then
      echo "Compilation of '${base_name}.cpp' with flags '${flags}' succeeded."
      
      # Run the compiled program and capture the output
      result=$(./"${dir}/${base_name}_${flags// /_}" $flags) 
      
      echo "$result"

      # Parse the execution time from the result
      exec_time=$(echo $result | awk -F ',' '{print $3}')
      
      if (( $(echo "$exec_time < ${bestTimes[$base_name]}" | bc -l) ))
      then
        bestTimes[$base_name]=$exec_time
        bestFlags[$base_name]=$flags
      fi
      
      rm "${dir}/${base_name}_${flags// /_}"
    else
      echo "Compilation with flags '${flags}' failed."
    fi
  done
done

# Display the summary
printf "\n\n"
printf "==============================================================================\n"
printf "%-25s %-10s %s\n" "FILE" "TIME (microseconds)" "FLAGS"
printf "==============================================================================\n"

for base_name in "${!bestTimes[@]}"
do
  printf "%-25s %-10s %s\n" "$base_name.cpp" "${bestTimes[$base_name]}" "'${bestFlags[$base_name]}'"
done

printf "==============================================================================\n"

SCRIPT_PATH="$(dirname "$(realpath "$0")")"
# Obtain the root directory path by removing /src from the script path
ROOT_PATH="${SCRIPT_PATH%/src}"

PYTHON_SCRIPT="$ROOT_PATH/plotresults.py"
FILES_PATH="$ROOT_PATH/../../results/csv/uniform_${test_type}.csv"
DESTINATION="$ROOT_PATH/../../results/images"

python3 $PYTHON_SCRIPT --csv_file $FILES_PATH --output_dir $DESTINATION