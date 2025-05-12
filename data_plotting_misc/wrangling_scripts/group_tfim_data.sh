#!/bin/bash

# Ensure the target directory exists
mkdir -p ../completed_fig3_data

# Loop through directories matching the pattern
for dir in prb_fig3_beta_*; do
    # Check if the directory exists and is a directory
    if [ -d "$dir" ]; then
        # Extract the value of xxx from the directory name
        beta_value=$(echo "$dir" | awk -F'_beta_' '{print $2}')

        # Check if xxx is smaller than 1.0
        if (( $(echo "$beta_value < 3.0" | bc -l) )); then
            # Navigate into the directory
            cd "$dir" || exit

            # Find all files matching the pattern and extract the largest yyy
            largest_file=$(ls PMRQMC_mpi-*.out 2>/dev/null | awk -F'-' '{print $2}' | awk -F'.' '{print $1}' | sort -n | tail -1)

            # If a largest file is found, copy it to the new location
            if [ -n "$largest_file" ]; then
                src_file="PMRQMC_mpi-${largest_file}.out"
                dest_file="../completed_fig3_data/beta_${beta_value}_data.out"
                cp "$src_file" "$dest_file"
                echo "Copied $src_file to $dest_file"
            fi

            # Return to the parent directory
            cd ..
        fi
    fi
done


