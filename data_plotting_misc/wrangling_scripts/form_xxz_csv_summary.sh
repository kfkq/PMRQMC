#!/bin/bash

# Output CSV file
output_file="fig5_xxz_jackknife_data.csv"

# Write the CSV header
echo "L,beta,Hoffdiag_value,Hoffdiag_std,Hoffdiag_Eint_value,Hoffdiag_Eint_std,Hoffdiag_Fint_value,Hoffdiag_Fint_std,mc_steps,avg_q,max_q,avg_sgn,sgn_std,tot_cpu_time (s),wall_time (s),Tsteps,steps,Hoffdiag_test,Hoffdiag_Eint_test,Hoffdiag_Fint_test,offdiag_ES_value,offdiag_ES_std,offidag_FS_value,offdiag_FS_std" > "$output_file"

# Loop through all beta_xxx_data.out files in the current directory
for file in L_*_beta_*.out; do
    # Check if the file exists
    if [ -f "$file" ]; then
	# Extract L value from the filename
        L=$(echo "$file" | sed -n 's/^L_\([0-9]*\)_beta_[0-9.]*_data.out$/\1/p')
        # Extract beta value from the filename
	beta=$(echo "$file" | sed -n 's/^L_[0-9]*_beta_\([0-9.]*\)_data.out$/\1/p')

        # Initialize variables for each observable, new data fields, and tests
        hoffdiag_value="" hoffdiag_std=""
        hoffdiag_eint_value="" hoffdiag_eint_std=""
        hoffdiag_fint_value="" hoffdiag_fint_std=""
        mc_steps="" avg_q="" max_q="" avg_sgn="" sgn_std="" tot_cpu_time="" wall_time=""
        tsteps="" steps=""
        hoffdiag_test="False"
        hoffdiag_eint_test="False" hoffdiag_fint_test="False"
	offdiag_es_value="" offdiag_es_std=""
	offdiag_fs_value="" offdiag_fs_std=""

        # Process the file to extract relevant data
        while IFS= read -r line; do
            if [[ "$line" == *"Total of observable #1: H_{offdiag}"* ]]; then
                read -r line; hoffdiag_value=$(echo "$line" | awk '{print $4}')
                read -r line; hoffdiag_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #2: measure_Hoffdiag_Eint"* ]]; then
                read -r line; hoffdiag_eint_value=$(echo "$line" | awk '{print $4}')
                read -r line; hoffdiag_eint_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #3: measure_Hoffdiag_Fint"* ]]; then
                read -r line; hoffdiag_fint_value=$(echo "$line" | awk '{print $4}')
                read -r line; hoffdiag_fint_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total number of MC updates = "* ]]; then
                mc_steps=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total mean(q) = "* ]]; then
                avg_q=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total max(q) = "* ]]; then
                max_q=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total mean(sgn(W)) = "* ]]; then
                avg_sgn=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total std.dev.(sgn(W)) = "* ]]; then
                sgn_std=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total elapsed cpu time = "* ]]; then
                tot_cpu_time=$(echo "$line" | awk -F'= ' '{print $2}' | awk '{print $1}')
            elif [[ "$line" == *"Wall-clock time = "* ]]; then
                wall_time=$(echo "$line" | awk -F'= ' '{print $2}' | awk '{print $1}')
            elif [[ "$line" == *"Parameters: beta = "* ]]; then
                tsteps=$(echo "$line" | awk -F', ' '{print $2}' | awk -F'= ' '{print $2}')
                steps=$(echo "$line" | awk -F', ' '{print $3}' | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Observable #1: H_{offdiag},"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_test="True"
            elif [[ "$line" == *"Observable #2: measure_Hoffdiag_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_eint_test="True"
            elif [[ "$line" == *"Observable #3: measure_Hoffdiag_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_fint_test="True"
	    elif [[ "$line" == *"Total of derived observable: offdiagonal energy susceptibility"* ]]; then
                read -r line; offdiag_ES_value=$(echo "$line" | awk '{print $4}')
                read -r line; offdiag_ES_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: offdiagonal fidelity susceptibility"* ]]; then
                read -r line; offdiag_FS_value=$(echo "$line" | awk '{print $4}')
                read -r line; offdiag_FS_std=$(echo "$line" | awk '{print $4}')
            fi
        done < "$file"

        # Append the extracted values to the CSV file
        echo "$L,$beta,$hoffdiag_value,$hoffdiag_std,$hoffdiag_eint_value,$hoffdiag_eint_std,$hoffdiag_fint_value,$hoffdiag_fint_std,$mc_steps,$avg_q,$max_q,$avg_sgn,$sgn_std,$tot_cpu_time,$wall_time,$tsteps,$steps,$hoffdiag_test,$hoffdiag_eint_test,$hoffdiag_fint_test,$offdiag_ES_value,$offdiag_ES_std,$offdiag_FS_value,$offdiag_FS_std" >> "$output_file"
    else
        echo "File $file does not exist."
    fi
done

echo "CSV generation completed: $output_file"

