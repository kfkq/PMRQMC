#!/bin/bash

# Output CSV file
output_file="fig3_tfim_jackknife_data.csv"

# Write the CSV header
echo "beta,Hdiag_value,Hdiag_std,Hoffdiag_value,Hoffdiag_std,Hdiag_Eint_value,Hdiag_Eint_std,Hdiag_Fint_value,Hdiag_Fint_std,Hoffdiag_Eint_value,Hoffdiag_Eint_std,Hoffdiag_Fint_value,Hoffdiag_Fint_std,mc_steps,avg_q,max_q,avg_sgn,sgn_std,tot_cpu_time (s),wall_time (s),Tsteps,steps,Hdiag_test,Hoffdiag_test,Hdiag_Eint_test,Hdiag_Fint_test,Hoffdiag_Eint_test,Hoffdiag_Fint_test,diag_ES_value,diag_ES_std,diag_FS_value,diag_FS_std,offdiag_ES_value,offdiag_ES_std,offdiag_FS_value,offdiag_FS_std" > "$output_file"

# Loop through all beta_xxx_data.out files in the current directory
for file in beta_*.out; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract beta value from the filename
        beta=$(echo "$file" | sed -n 's/^beta_\([0-9.]*\)_data.out$/\1/p')

        # Initialize variables for each observable, new data fields, and tests
        hdiag_value="" hdiag_std=""
        hoffdiag_value="" hoffdiag_std=""
        hdiag_eint_value="" hdiag_eint_std=""
        hdiag_fint_value="" hdiag_fint_std=""
        hoffdiag_eint_value="" hoffdiag_eint_std=""
        hoffdiag_fint_value="" hoffdiag_fint_std=""
        mc_steps="" avg_q="" max_q="" avg_sgn="" sgn_std="" tot_cpu_time="" wall_time=""
        tsteps="" steps=""
        hdiag_test="False" hoffdiag_test="False"
        hdiag_eint_test="False" hdiag_fint_test="False"
        hoffdiag_eint_test="False" hoffdiag_fint_test="False"
	diag_es_value="" diag_es_std=""
	diag_fs_value="" diag_fs_std=""
	offdiag_es_value="" offdiag_es_std=""
	offdiag_fs_value="" offdiag_fs_std=""

        # Process the file to extract relevant data
        while IFS= read -r line; do
            if [[ "$line" == *"Total of observable #1: H_{diag}"* ]]; then
                read -r line; hdiag_value=$(echo "$line" | awk '{print $4}')
                read -r line; hdiag_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #2: H_{offdiag}"* ]]; then
                read -r line; hoffdiag_value=$(echo "$line" | awk '{print $4}')
                read -r line; hoffdiag_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #3: measure_Hdiag_Eint"* ]]; then
                read -r line; hdiag_eint_value=$(echo "$line" | awk '{print $4}')
                read -r line; hdiag_eint_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #4: measure_Hdiag_Fint"* ]]; then
                read -r line; hdiag_fint_value=$(echo "$line" | awk '{print $4}')
                read -r line; hdiag_fint_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #5: measure_Hoffdiag_Eint"* ]]; then
                read -r line; hoffdiag_eint_value=$(echo "$line" | awk '{print $4}')
                read -r line; hoffdiag_eint_std=$(echo "$line" | awk '{print $4}')
            elif [[ "$line" == *"Total of observable #6: measure_Hoffdiag_Fint"* ]]; then
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
            elif [[ "$line" == *"Observable #1: H_{diag},"* ]]; then
                [[ "$line" == *"test passed"* ]] && hdiag_test="True"
            elif [[ "$line" == *"Observable #2: H_{offdiag},"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_test="True"
            elif [[ "$line" == *"Observable #3: measure_Hdiag_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hdiag_eint_test="True"
            elif [[ "$line" == *"Observable #4: measure_Hdiag_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hdiag_fint_test="True"
            elif [[ "$line" == *"Observable #5: measure_Hoffdiag_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_eint_test="True"
            elif [[ "$line" == *"Observable #6: measure_Hoffdiag_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hoffdiag_fint_test="True"
	    elif [[ "$line" == *"Total of derived observable: diagonal energy susceptibility"* ]]; then
                read -r line; diag_ES_value=$(echo "$line" | awk '{print $4}')
                read -r line; diag_ES_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: diagonal fidelity susceptibility"* ]]; then
                read -r line; diag_FS_value=$(echo "$line" | awk '{print $4}')
                read -r line; diag_FS_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: offdiagonal energy susceptibility"* ]]; then
                read -r line; offdiag_ES_value=$(echo "$line" | awk '{print $4}')
                read -r line; offdiag_ES_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: offdiagonal fidelity susceptibility"* ]]; then
                read -r line; offdiag_FS_value=$(echo "$line" | awk '{print $4}')
                read -r line; offdiag_FS_std=$(echo "$line" | awk '{print $4}')
            fi
        done < "$file"

        # Append the extracted values to the CSV file
        echo "$beta,$hdiag_value,$hdiag_std,$hoffdiag_value,$hoffdiag_std,$hdiag_eint_value,$hdiag_eint_std,$hdiag_fint_value,$hdiag_fint_std,$hoffdiag_eint_value,$hoffdiag_eint_std,$hoffdiag_fint_value,$hoffdiag_fint_std,$mc_steps,$avg_q,$max_q,$avg_sgn,$sgn_std,$tot_cpu_time,$wall_time,$tsteps,$steps,$hdiag_test,$hoffdiag_test,$hdiag_eint_test,$hdiag_fint_test,$hoffdiag_eint_test,$hoffdiag_fint_test,$diag_ES_value,$diag_ES_std,$diag_FS_value,$diag_FS_std,$offdiag_ES_value,$offdiag_ES_std,$offdiag_FS_value,$offdiag_FS_std" >> "$output_file"
    else
        echo "File $file does not exist."
    fi
done

echo "CSV generation completed: $output_file"

