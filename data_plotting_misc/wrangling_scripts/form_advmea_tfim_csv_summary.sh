#!/bin/bash

# Output CSV file
output_file="advmeas_tfim_data.csv"

# Write the CSV header
echo "L,lam,Tsteps,steps,tot_steps,avg_sgn,avg_sgn_std,avg_q,max_q,a,a_std,a2,a2_std,aCorr,aCorr_std,aEint,aEint_std,aFint,aFint_std,b,b_std,b2,b2_std,bCorr,bCorr_std,bEint,bEint_std,bFint,bFint_std,ab,ab_std,abCorr,abCorr_std,abEint,abEint_std,abFint,abFint_std,h,h_std,h2,h2_std,hDiag,hDiag_std,hDiag2,hDiag2_std,hOffdiag,hOffdiag_std,hOffdiag2,hOffdiag2_std,zmag,zmag_std,hDiagCorr,hDiagCorr_std,hDiagEint,hDiagEint_std,hDiagFint,hDiagFint_std,hOffdiagCorr,hOffdiagCorr_std,hOffdiagEint,hOffdiagEint_std,hOffdiagFint,hOffdiagFint_std,diagES,diagES_std,diagFS,diagFS_std,offDiagES,offDiagES_std,offDiagFS,offDiagFS_std,Cv,Cv_std,a_therm,a2_therm,aCorr_therm,aEint_therm,aFint_therm,b_therm,b2_therm,bCorr_therm,bEint_therm,bFint_therm,ab_therm,abCorr_therm,abEint_therm,abFint_therm,h_therm,h2_therm,hDiag_therm,hDiag2_therm,hOffdiag_therm,hOffdiag2_therm,zmag_therm,hDiagCorr_therm,hDiagEint_therm,hDiagFint_therm,hOffdiagCorr_therm,hOffdiagEint_therm,hOffdiagFint_therm,diagES_therm,diagFS_therm,offDiagES_therm,offDiagFS_therm,Cv_therm,tot_cpu_time (s),wall_time (s)" >> "$output_file"

# Loop through all lam_xxx_data.out files in the current directory
for file in L_*_lam_*.out; do
    # Check if the file exists
    if [ -f "$file" ]; then
	# Extract L value from the filename
        L=$(echo "$file" | sed -n 's/^L_\([0-9]*\)_lam_[0-9.]*_data.out$/\1/p')
        # Extract lam value from the filename
	lam=$(echo "$file" | sed -n 's/^L_[0-9]*_lam_\([0-9.]*\)_data.out$/\1/p')

        # Initialize variables for each observable, new data fields, and tests
	# basic meta-parameters
	Tsteps="" steps="" tot_steps=""
	avg_sgn="" avg_sgn_std=""
	avg_q="" max_q=""
	# A outcomes 
	a="" a_std="" a2="" a2_std="" aCorr="" aCorr_std=""
	aEint="" aEint_std="" aFint="" aFint_std=""
	# B outcomes
	b="" b_std="" b2="" b2_std="" bCorr="" bCorr_std=""
	bEint="" bEint_std="" bFint="" bFint_std=""
	# AB outomes
	ab="" ab_std="" abCorr="" abCorr_std=""
	abEint="" abEint_std="" abFint="" abFint_std=""
	# Basic simple observables like H, Hdiag
	h="" h_std="" h2="" h2_std="" hDiag="" hDiag_std=""
	hDiag2="" hDiag2_std="" hOffDiag="" hOffDiag_std=""
	hOffDiag2="" hOffDiag2_std="" zmag="" zmag_std=""
	# Correlators, Eints and Fints of Hdiag,HOffidag
	hDiagCorr="" hDiagCorr_std="" hDiagEint="" hDiagEint_std=""
	hDiagFint="" hDiagFint_std="" hOffdiagCorr="" hOffdiagCorr_std=""
	hOffdiagEint="" hOffdiagEint_std="" hOffdiagFint="" hOffdiagFint_std=""
	# Derived observables
	diagES="" diagES_std="" diagFS="" diagFS_std=""
	offDiagES="" offDiagES_std="" offDiagFS="" offDiagFS_std=""
	Cv="" Cv_std=""
	# time meta-data
	tot_cpu_time="" wall_time=""
	# A therm tests
	a_therm="False" a2_therm="False" aCorr_therm="False"
	aEint_therm="False" aFint_therm="False"
	# B therm tests
	b_therm="False" b2_therm="False" bCorr_therm="False"
	bEint_therm="False" bFint_therm="False"
	# AB therm tests
	ab_therm="False" abCorr_therm="False"
	abEint_therm="False" abFint_therm="False"
	# Basic H therm tests
	h_therm="False" h2_therm="False"
	hDiag_therm="False" hDiag2_therm="False"
	hOffdiag_therm="False" hOffdiag2_therm="False"
	zmag_therm="False"
	# Complicated Hdiag, Hoffdiag therm tests
	hDiagCorr_therm="False"
	hDiagEint_therm="False"
	hDiagFint_therm="False"
	hOffdiagCorr_therm="False"
	hOffdiagEint_therm="False"
	hOffdiagFint_therm="False"
	# Derived observable therm tests
	diagES_therm="False"
	diagFS_therm="False"
	offDiagES_therm="False"
	offDiagFS_therm="False"
	Cv_therm="False"

        # Process the file to extract relevant data
        while IFS= read -r line; do
	    # meta-data (not directly an observable, basically)
            if [[ "$line" == *"Parameters: beta = "* ]]; then
                Tsteps=$(echo "$line" | awk -F', ' '{print $2}' | awk -F'= ' '{print $2}')
                steps=$(echo "$line" | awk -F', ' '{print $3}' | awk -F'= ' '{print $2}')
	    elif [[ "$line" == *"Total number of MC updates = "* ]]; then
                tot_steps=$(echo "$line" | awk -F'= ' '{print $2}')
	    elif [[ "$line" == *"Total mean(sgn(W)) = "* ]]; then
                avg_sgn=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total std.dev.(sgn(W)) = "* ]]; then
                avg_sgn_std=$(echo "$line" | awk -F'= ' '{print $2}')
	    elif [[ "$line" == *"Total mean(q) = "* ]]; then
                avg_q=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total max(q) = "* ]]; then
                max_q=$(echo "$line" | awk -F'= ' '{print $2}')
            elif [[ "$line" == *"Total elapsed cpu time = "* ]]; then
                tot_cpu_time=$(echo "$line" | awk -F'= ' '{print $2}' | awk '{print $1}')
            elif [[ "$line" == *"Wall-clock time = "* ]]; then
                wall_time=$(echo "$line" | awk -F'= ' '{print $2}' | awk '{print $1}')
	    # A-based observables
	    elif [[ "$line" == *"Total of observable #1: A"* ]]; then
                read -r line; a=$(echo "$line" | awk '{print $4}')
                read -r line; a_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #2: A^2"* ]]; then
                read -r line; a2=$(echo "$line" | awk '{print $4}')
                read -r line; a2_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #3: A(tau)A"* ]]; then
                read -r line; aCorr=$(echo "$line" | awk '{print $4}')
                read -r line; aCorr_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #4: A_Eint"* ]]; then
                read -r line; aEint=$(echo "$line" | awk '{print $4}')
                read -r line; aEint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #5: A_Fint"* ]]; then
                read -r line; aFint=$(echo "$line" | awk '{print $4}')
                read -r line; aFint_std=$(echo "$line" | awk '{print $4}')
	    # B-based observables
	    elif [[ "$line" == *"Total of observable #6: B"* ]]; then
                read -r line; b=$(echo "$line" | awk '{print $4}')
                read -r line; b_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #7: B^2"* ]]; then
                read -r line; b2=$(echo "$line" | awk '{print $4}')
                read -r line; b2_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #8: B(tau)B"* ]]; then
                read -r line; bCorr=$(echo "$line" | awk '{print $4}')
                read -r line; bCorr_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #9: B_Eint"* ]]; then
                read -r line; bEint=$(echo "$line" | awk '{print $4}')
                read -r line; bEint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #10: B_Fint"* ]]; then
                read -r line; bFint=$(echo "$line" | awk '{print $4}')
                read -r line; bFint_std=$(echo "$line" | awk '{print $4}')
	    # AB-based observables
	    elif [[ "$line" == *"Total of observable #11: AB"* ]]; then
                read -r line; ab=$(echo "$line" | awk '{print $4}')
                read -r line; ab_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #12: A(tau)B"* ]]; then
                read -r line; abCorr=$(echo "$line" | awk '{print $4}')
                read -r line; abCorr_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #13: AB_Eint"* ]]; then
                read -r line; abEint=$(echo "$line" | awk '{print $4}')
                read -r line; abEint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #14: AB_Fint"* ]]; then
                read -r line; abFint=$(echo "$line" | awk '{print $4}')
                read -r line; abFint_std=$(echo "$line" | awk '{print $4}')
	    # H-based observables
	    elif [[ "$line" == *"Total of observable #15: H"* ]]; then
                read -r line; h=$(echo "$line" | awk '{print $4}')
                read -r line; h_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #16: H^2"* ]]; then
                read -r line; h2=$(echo "$line" | awk '{print $4}')
                read -r line; h2_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #17: H_{diag}"* ]]; then
                read -r line; hDiag=$(echo "$line" | awk '{print $4}')
                read -r line; hDiag_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #18: H_{diag}^2"* ]]; then
                read -r line; hDiag2=$(echo "$line" | awk '{print $4}')
                read -r line; hDiag2_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #19: H_{offdiag}"* ]]; then
                read -r line; hOffdiag=$(echo "$line" | awk '{print $4}')
                read -r line; hOffdiag_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #20: H_{offdiag}^2"* ]]; then
                read -r line; hOffdiag2=$(echo "$line" | awk '{print $4}')
                read -r line; hOffdiag2_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #21: Z_magnetization"* ]]; then
                read -r line; zmag=$(echo "$line" | awk '{print $4}')
                read -r line; zmag_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #22: measure_Hdiag_corr"* ]]; then
                read -r line; hDiagCorr=$(echo "$line" | awk '{print $4}')
                read -r line; hDiagCorr_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #23: measure_Hdiag_Eint"* ]]; then
                read -r line; hDiagEint=$(echo "$line" | awk '{print $4}')
                read -r line; hDiagEint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #24: measure_Hdiag_Fint"* ]]; then
                read -r line; hDiagFint=$(echo "$line" | awk '{print $4}')
                read -r line; hDiagFint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #25: measure_Hoffdiag_corr"* ]]; then
                read -r line; hOffdiagCorr=$(echo "$line" | awk '{print $4}')
                read -r line; hOffdiagCorr_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #26: measure_Hoffdiag_Eint"* ]]; then
                read -r line; hOffdiagEint=$(echo "$line" | awk '{print $4}')
                read -r line; hOffdiagEint_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of observable #27: measure_Hoffdiag_Fint"* ]]; then
                read -r line; hOffdiagFint=$(echo "$line" | awk '{print $4}')
                read -r line; hOffdiagFint_std=$(echo "$line" | awk '{print $4}')
	    # Derived observables
	    elif [[ "$line" == *"Total of derived observable: diagonal energy susceptibility"* ]]; then
                read -r line; diagES=$(echo "$line" | awk '{print $4}')
                read -r line; diagES_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: diagonal fidelity susceptibility"* ]]; then
                read -r line; diagFS=$(echo "$line" | awk '{print $4}')
                read -r line; diagFS_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: offdiagonal energy susceptibility"* ]]; then
                read -r line; offDiagES=$(echo "$line" | awk '{print $4}')
                read -r line; offDiagES_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: offdiagonal fidelity susceptibility"* ]]; then
                read -r line; offDiagFS=$(echo "$line" | awk '{print $4}')
                read -r line; offDiagFS_std=$(echo "$line" | awk '{print $4}')
	    elif [[ "$line" == *"Total of derived observable: specific heat"* ]]; then
                read -r line; Cv=$(echo "$line" | awk '{print $4}')
                read -r line; Cv_std=$(echo "$line" | awk '{print $4}')
	    # thermalization tests
	    elif [[ "$line" == *"Observable #1: A,"* ]]; then
                [[ "$line" == *"test passed"* ]] && a_therm="True"
	    elif [[ "$line" == *"Observable #2: A^2,"* ]]; then
                [[ "$line" == *"test passed"* ]] && a2_therm="True"
	    elif [[ "$line" == *"Observable #3: A(tau)A,"* ]]; then
                [[ "$line" == *"test passed"* ]] && aCorr_therm="True"
	    elif [[ "$line" == *"Observable #4: A_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && aEint_therm="True"
	    elif [[ "$line" == *"Observable #5: A_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && aFint_therm="True"
	    elif [[ "$line" == *"Observable #6: B,"* ]]; then
                [[ "$line" == *"test passed"* ]] && b_therm="True"
	    elif [[ "$line" == *"Observable #7: B^2,"* ]]; then
                [[ "$line" == *"test passed"* ]] && b2_therm="True"
	    elif [[ "$line" == *"Observable #8: B(tau)B,"* ]]; then
                [[ "$line" == *"test passed"* ]] && bCorr_therm="True"
	    elif [[ "$line" == *"Observable #9: B_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && bEint_therm="True"
	    elif [[ "$line" == *"Observable #10: B_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && bFint_therm="True"
	    elif [[ "$line" == *"Observable #11: AB,"* ]]; then
                [[ "$line" == *"test passed"* ]] && ab_therm="True"
	    elif [[ "$line" == *"Observable #12: A(tau)B,"* ]]; then
                [[ "$line" == *"test passed"* ]] && abCorr_therm="True"
	    elif [[ "$line" == *"Observable #13: AB_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && abEint_therm="True"
	    elif [[ "$line" == *"Observable #14: AB_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && abFint_therm="True"
	    elif [[ "$line" == *"Observable #15: H,"* ]]; then
                [[ "$line" == *"test passed"* ]] && h_therm="True"
	    elif [[ "$line" == *"Observable #16: H^2,"* ]]; then
                [[ "$line" == *"test passed"* ]] && h2_therm="True"
	    elif [[ "$line" == *"Observable #17: H_{diag},"* ]]; then
                [[ "$line" == *"test passed"* ]] && hDiag_therm="True"
	    elif [[ "$line" == *"Observable #18: H_{diag}^2,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hDiag2_therm="True"
	    elif [[ "$line" == *"Observable #19: H_{offdiag},"* ]]; then
                [[ "$line" == *"test passed"* ]] && hOffdiag_therm="True"
	    elif [[ "$line" == *"Observable #20: H_{offdiag}^2,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hOffdiag2_therm="True"
	    elif [[ "$line" == *"Observable #21: Z_magnetization,"* ]]; then
                [[ "$line" == *"test passed"* ]] && zmag_therm="True"
	    elif [[ "$line" == *"Observable #22: measure_Hdiag_corr,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hDiagCorr_therm="True"
	    elif [[ "$line" == *"Observable #23: measure_Hdiag_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hDiagEint_therm="True"
	    elif [[ "$line" == *"Observable #24: measure_Hdiag_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hDiagFint_therm="True"
	    elif [[ "$line" == *"Observable #25: measure_Hoffdiag_corr,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hOffdiagCorr_therm="True"
	    elif [[ "$line" == *"Observable #26: measure_Hoffdiag_Eint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hOffdiagEint_therm="True"
	    elif [[ "$line" == *"Observable #27: measure_Hoffdiag_Fint,"* ]]; then
                [[ "$line" == *"test passed"* ]] && hOffdiagFint_therm="True"
	    elif [[ "$line" == *"Derived observable: diagonal energy susceptibility,"* ]]; then
                [[ "$line" == *"test passed"* ]] && diagES_therm="True"
	    elif [[ "$line" == *"Derived observable: diagonal fidelity susceptibility,"* ]]; then
                [[ "$line" == *"test passed"* ]] && diagFS_therm="True"
	    elif [[ "$line" == *"Derived observable: offdiagonal energy susceptibility,"* ]]; then
                [[ "$line" == *"test passed"* ]] && offDiagES_therm="True"
	    elif [[ "$line" == *"Derived observable: offdiagonal fidelity susceptibility,"* ]]; then
                [[ "$line" == *"test passed"* ]] && offDiagFS_therm="True"
	    elif [[ "$line" == *"Derived observable: specific heat,"* ]]; then
                [[ "$line" == *"test passed"* ]] && Cv_therm="True"


            fi
        done < "$file"
        # Append the extracted values to the CSV file
        echo "$L,$lam,$Tsteps,$steps,$tot_steps,$avg_sgn,$avg_sgn_std,$avg_q,$max_q,$a,$a_std,$a2,$a2_std,$aCorr,$aCorr_std,$aEint,$aEint_std,$aFint,$aFint_std,$b,$b_std,$b2,$b2_std,$bCorr,$bCorr_std,$bEint,$bEint_std,$bFint,$bFint_std,$ab,$ab_std,$abCorr,$abCorr_std,$abEint,$abEint_std,$abFint,$abFint_std,$h,$h_std,$h2,$h2_std,$hDiag,$hDiag_std,$hDiag2,$hDiag2_std,$hOffdiag,$hOffdiag_std,$hOffdiag2,$hOffdiag2_std,$zmag,$zmag_std,$hDiagCorr,$hDiagCorr_std,$hDiagEint,$hDiagEint_std,$hDiagFint,$hDiagFint_std,$hOffdiagCorr,$hOffdiagCorr_std,$hOffdiagEint,$hOffdiagEint_std,$hOffdiagFint,$hOffdiagFint_std,$diagES,$diagES_std,$diagFS,$diagFS_std,$offDiagES,$offDiagES_std,$offDiagFS,$offDiagFS_std,$Cv,$Cv_std,$a_therm,$a2_therm,$aCorr_therm,$aEint_therm,$aFint_therm,$b_therm,$b2_therm,$bCorr_therm,$bEint_therm,$bFint_therm,$ab_therm,$abCorr_therm,$abEint_therm,$abFint_therm,$h_therm,$h2_therm,$hDiag_therm,$hDiag2_therm,$hOffdiag_therm,$hOffdiag2_therm,$zmag_therm,$hDiagCorr_therm,$hDiagEint_therm,$hDiagFint_therm,$hOffdiagCorr_therm,$hOffdiagEint_therm,$hOffdiagFint_therm,$diagES_therm,$diagFS_therm,$offDiagES_therm,$offDiagFS_therm,$Cv_therm,$tot_cpu_time,$wall_time" >> "$output_file"
    else
        echo "File $file does not exist."
    fi
done

echo "CSV generation completed: $output_file"

