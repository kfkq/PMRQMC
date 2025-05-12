#!/bin/bash

# Ensure the target directory exists
mkdir -p completed_fig5_data

# Define ranges for L and beta
L_min=18
L_max=30
beta_min=1.17
beta_max=2.22

# Loop through directories matching the pattern
for dir in prx_fig5_L_*_beta_*; do
    # Check if the directory exists and is a directory
    if [ -d "$dir" ]; then
        # Extract the value of L from the directory name
        L_value=$(echo "$dir" | awk -F'_L_' '{print $2}' | awk -F'_beta_' '{print $1}')
        # Extract the value of beta from the directory name
        beta_value=$(echo "$dir" | awk -F'_beta_' '{print $2}')

        # Check if L and beta are within the specified ranges
        if (( L_value >= L_min && L_value <= L_max )) && (( $(echo "$beta_value >= $beta_min && $beta_value <= $beta_max" | bc -l) )); then
            # Navigate into the directory
            cd "$dir" || exit

            # Find all files matching the pattern and extract the largest yyy
            largest_file=$(ls PMRQMC_mpi-*.out 2>/dev/null | awk -F'-' '{print $2}' | awk -F'.' '{print $1}' | sort -n | tail -1)

            # If a largest file is found, copy it to the new location
            if [ -n "$largest_file" ]; then
                src_file="PMRQMC_mpi-${largest_file}.out"
                dest_file="../completed_fig5_data/L_${L_value}_beta_${beta_value}_data.out"
                cp "$src_file" "$dest_file"
                echo "Copied $src_file to $dest_file"
            fi

            # Return to the parent directory
            cd ..
        fi
    fi
done
