#!/bin/bash

# Ensure the target directory exists
mkdir -p completed_advmea_tfim_data

# Define ranges for L and beta
L_min=4
L_max=5
lam_min=1.0
lam_max=5.0

# Loop through directories matching the pattern
for dir in square_tfim_L_*_lam_*; do
    # Check if the directory exists and is a directory
    if [ -d "$dir" ]; then
        # Extract the value of L from the directory name
        L_value=$(echo "$dir" | awk -F'_L_' '{print $2}' | awk -F'_lam_' '{print $1}')
        # Extract the value of beta from the directory name
        lam_value=$(echo "$dir" | awk -F'_lam_' '{print $2}')

        # Check if L and beta are within the specified ranges
        if (( L_value >= L_min && L_value <= L_max )) && (( $(echo "$lam_value >= $lam_min && $lam_value <= $lam_max" | bc -l) )); then
            # Navigate into the directory
            cd "$dir" || exit

            # Find all files matching the pattern and extract the largest yyy
            largest_file=$(ls PMRQMC_mpi-*.out 2>/dev/null | awk -F'-' '{print $2}' | awk -F'.' '{print $1}' | sort -n | tail -1)

            # If a largest file is found, copy it to the new location
            if [ -n "$largest_file" ]; then
                src_file="PMRQMC_mpi-${largest_file}.out"
                dest_file="../completed_advmea_tfim_data/L_${L_value}_lam_${lam_value}_data.out"
                cp "$src_file" "$dest_file"
                echo "Copied $src_file to $dest_file"
            fi

            # Return to the parent directory
            cd ..
        fi
    fi
done
