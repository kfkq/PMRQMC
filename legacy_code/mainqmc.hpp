//
// This program implements estimators for the titular quantities in the preprint:
// Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility.
// 
// This work augments the base code whose information is below.
//
//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians, Physical Review Research 6, 013281 (2024).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include<iostream>
#include<iomanip>
#include<complex>
#include<random>
#include<cstdlib>
#include<algorithm>
#include<bitset>
#include"divdiff.hpp"
#include"hamiltonian.hpp" // use a header file, which defines the Hamiltonian and the custom observables
#include"parameters.hpp"  // parameters of the simulation such as the number of Monte-Carlo updates

#define measurements (steps/stepsPerMeasurement)

#ifdef EXHAUSTIVE_CYCLE_SEARCH
	#define rmin 0                            // length r of sub-sequence is chosen randomly between rmin and rmax
	#define rmax cycle_max_len
	#define lmin r                            // cycle lengths must be between lmin and lmax
	#define lmax cycle_max_len
#else
	#define rmin (cycle_min_len-1)/2
	#define rmax (cycle_max_len+1)/2
	#define lmin 2*r-1
	#define lmax 2*r+1
#endif

static std::random_device rd;
static std::mt19937 rng;
static std::uniform_int_distribution<> dice2(0,1);
static std::uniform_int_distribution<> diceN(0,N-1);
static std::uniform_int_distribution<> diceNop(0,Nop-1);
static std::uniform_real_distribution<> val(0.0,1.0);
static std::geometric_distribution<> geometric_int(0.8);

ExExFloat beta_pow_factorial[qmax]; // contains the values (-beta)^q / q!
ExExFloat beta_div2_pow_factorial[qmax]; // contains the values (-beta/2)^q / q!
ExExFloat tau_pow_factorial[qmax]; // contains the values (-tau)^q / q!
ExExFloat tau_minus_beta_pow_factorial[qmax]; // contains the values (tau-beta)^q / q!
double factorial[qmax]; // contains the values q!
int cycle_len[Ncycles];
int cycles_used[Ncycles];
int n_cycles[Nop+3]; // numbers of cycles of lengths 0,1,2,...,Nop+2, the last three values are always zeros.
int cycle_min_len, cycle_max_len, found_cycles, min_index, max_index;

#ifndef MEASURE_CUSTOM_OBSERVABLES
#define Nobservables 0
#endif

const int N_all_observables = Nobservables + 15;
int valid_observable[N_all_observables];

int bin_length = measurements / Nbins;
double in_bin_sum[N_all_observables];
double bin_mean[N_all_observables][Nbins];
double in_bin_sum_sgn;
double bin_mean_sgn[Nbins];

int q;
int qmax_achieved=0;

divdiff* d; // primary QMC divided difference of [-\beta Ez_0,...,-\beta Ez_q]
divdiff* dfs; // QMC divided difference for fidsus [(-\beta/2)Ez_0,...,(-\beta/2)Ez_q]
divdiff* ds1; // a scratch divided diff for measurements
divdiff* ds2; // a scratch divided diff for measurements

std::bitset<N> lattice;
std::bitset<N> z;
std::bitset<Nop> P;
std::bitset<Nop> Pk, Pl;

int Sq[qmax];	// list of q operator indices
int Sq_backup[qmax];
int Sq_subseq[qmax];
int Sq_gaps[qmax];
double Energies[qmax+1];
double Energies_backup[qmax+1];
int eoccupied[qmax+1];
double currEnergy;
std::complex<double> old_currD, currD, currMDk, currMDl, currMD0;
std::complex<double> currD_partial[qmax], currMDk_trace[qmax], currMDl_trace[qmax], currMD0_trace[qmax];

ExExFloat one, currWeight;
unsigned long long step;
unsigned long long measurement_step;

double CalcEnergy(){ // calculate the energy <z | D_0 | z> of a given configuration of spins
	std::complex<double> sum = 0;
	for(int i=0;i<D0_size;i++) sum -= double(2*(int((D0_product[i] & (~lattice)).count())%2)-1) * D0_coeff[i];
	return sum.real();
}

std::complex<double> calc_d(int k){ // calculate d_k = <z | D_k | z> for the current configuration of spins
	std::complex<double> sum = 0;
	for(int i=0;i<D_size[k];i++) sum -= double(2*(int((D_product[k][i] & (~lattice)).count())%2)-1) * D_coeff[k][i];
	return sum;
}

void ApplyOperator(int k){
	lattice ^= P_matrix[k];
}

void GetEnergies(){
	currD = currD_partial[0] = 1;
	for(int i=0;i<q;i++){
		Energies[i] = CalcEnergy();
		ApplyOperator(Sq[i]);
		currD *= calc_d(Sq[i]);
		currD_partial[i+1] = currD; // all intermediate products as in eq. 32 from 0 to i
	}
	currEnergy = Energies[q] = CalcEnergy();
}

ExExFloat GetWeight(){
	d->CurrentLength=0; GetEnergies();
	for(int i=0;i<=q;i++) d->AddElement(-beta*Energies[i]);
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();
}

ExExFloat UpdateWeight(){
	int i, j, notfound, n=d->CurrentLength; double value;
	GetEnergies(); memset(eoccupied,0,(q+1)*sizeof(int));
	for(i=0;i<n;i++){
	        notfound = 1; value = d->z[i];
		for(j=0;j<=q;j++) if(eoccupied[j]==0 && value == -beta*Energies[j]){ notfound = 0; break; }
		if(notfound) break; eoccupied[j] = 1;
	}
	if(i==0) d->CurrentLength=0; else while(n>i){ d->RemoveElement(); n--; }
	j=0; while(i<=q){ while(eoccupied[j]) j++; d->AddElement(-beta*Energies[j++]); i++; }
        return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();
}

ExExFloat UpdateWeightReplace(double removeEnergy, double addEnergy){
	if(removeEnergy != addEnergy) if(d->RemoveValue(-beta*removeEnergy)) d->AddElement(-beta*addEnergy); else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real(); // use this value only when the values of q and currD are correct
}

ExExFloat UpdateWeightDel(double removeEnergy1, double removeEnergy2){
	if(d->RemoveValue(-beta*removeEnergy1) && d->RemoveValue(-beta*removeEnergy2))
		return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();  // use this value only when the values of q and currD are correct
	else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
}

ExExFloat UpdateWeightIns(double addEnergy1, double addEnergy2){
	d->AddElement(-beta*addEnergy1); d->AddElement(-beta*addEnergy2);
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();    // use this value only when the values of q and currD are correct
}

int NoRepetitionCheck(int* sequence, int r){ // check for absence of repetitions in a sequence of length r
	int i,j,rep = 1;
	for(i=0;i<r && rep;i++) for(j=0;j<i;j++) if(sequence[j]==sequence[i]){ rep = 0; break;}
	return rep;
}

void PickSubsequence(int r){ // randomly picks a sequential sub-sequence of length r from Sq
	int i,m; m = int(val(rng)*(q-r+1)); // m is random integer between 0 and q-r
	for(i=0;i<r;i++) Sq_subseq[i] = Sq[i+m];
	min_index = m; max_index = m+r-1;
}

int FindCycles(int r){  // find all cycles of length between lmin and lmax, each containing all operators of the array Sq_subseq of length r.
	int i,j,not_contained;
	int found_cycle_list[Ncycles];
	found_cycles = 0;                                  // found_cycles contains the number of cycles found. it is global variable.
	for(i=0;i<Ncycles;i++){
		if(cycle_len[i]<lmin || cycle_len[i]>lmax) continue;
		not_contained = 0; for(j=0;j<r;j++) if(!cycles[i].test(Sq_subseq[j])){ not_contained = 1; break;}
		if(not_contained) continue;
		found_cycle_list[found_cycles++] = i;
	}
	return found_cycles>0 ? found_cycle_list[int(val(rng)*found_cycles)] : -1; // returns one of the found cycles chosen randomly.
}

//unsigned int rng_seed;

void init(){
	int i; double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	ExExFloat curr_tau; ExExFloat curr_bmt; ExExFloat curr_beta2;
	for(q=1;q<qmax;q++){ 
		curr1*=(-beta)/q; curr2*=q; beta_pow_factorial[q] = curr1; factorial[q] = curr2;
		curr_tau*=(-tau)/q; tau_pow_factorial[q] = curr_tau;
		curr_bmt*=(tau-beta)/q; tau_minus_beta_pow_factorial[q] = curr_bmt;
		curr_beta2*=(-beta/(2*q)); beta_div2_pow_factorial[q] = curr_beta2;
	}
	rng.seed(rng_seed);
	//rng_seed = rd(); rng.seed(rng_seed);
	//rng.seed(100);
	lattice = 0; for(i=N-1;i>=0;i--) if(dice2(rng)) lattice.set(i); z = lattice; q=0;
	currWeight = GetWeight();
	for(i=0;i<Ncycles;i++) cycle_len[i] = cycles[i].count();
	cycle_min_len = 64; for(i=0;i<Ncycles;i++) cycle_min_len = min(cycle_min_len,cycle_len[i]);
	cycle_max_len = 0; for(i=0;i<Ncycles;i++) cycle_max_len = max(cycle_max_len,cycle_len[i]);
	for(i=0;i<Ncycles;i++) cycles_used[i] = 0;
	for(i=0;i<Nop+3;i++) n_cycles[i] = 0; for(i=0;i<Ncycles;i++) n_cycles[cycle_len[i]]++;
	for(i=0;i<N_all_observables;i++) in_bin_sum[i] = 0; in_bin_sum_sgn = 0;
	for(i=0;i<N_all_observables;i++) valid_observable[i] = 0;
#ifdef MEASURE_CUSTOM_OBSERVABLES
	for(i=0;i<Nobservables;i++) valid_observable[i] = 1;
#endif
#ifdef MEASURE_H
	valid_observable[Nobservables] = 1;
#endif
#ifdef MEASURE_H2
	valid_observable[Nobservables + 1] = 1;
#endif
#ifdef MEASURE_HDIAG
	valid_observable[Nobservables + 2] = 1;
#endif
#ifdef MEASURE_HDIAG2
	valid_observable[Nobservables + 3] = 1;
#endif
#ifdef MEASURE_HOFFDIAG
	valid_observable[Nobservables + 4] = 1;
#endif
#ifdef MEASURE_HOFFDIAG2
	valid_observable[Nobservables + 5] = 1;
#endif
#ifdef MEASURE_Z_MAGNETIZATION
	valid_observable[Nobservables + 6] = 1;
#endif
#ifdef MEASURE_HDIAG_CORR
	valid_observable[Nobservables + 7] = 1;
#endif
#ifdef MEASURE_HDIAG_INT1
	valid_observable[Nobservables + 8] = 1;
#endif
#ifdef MEASURE_HDIAG_INT2
	valid_observable[Nobservables + 9] = 1;
#endif
#ifdef MEASURE_HDIAG_INT3
	valid_observable[Nobservables + 10] = 1;
#endif
#ifdef MEASURE_HOFFDIAG_CORR
	valid_observable[Nobservables + 11] = 1;
#endif
#ifdef MEASURE_HOFFDIAG_INT1
	valid_observable[Nobservables + 12] = 1;
#endif
#ifdef MEASURE_HOFFDIAG_INT2
	valid_observable[Nobservables + 13] = 1;
#endif
#ifdef MEASURE_HOFFDIAG_INT3
	valid_observable[Nobservables + 14] = 1;
#endif
}

double Metropolis(ExExFloat newWeight){
	return min(1.0,fabs((newWeight/currWeight).get_double()));
}

void update(){
	int i,m,p,r,u; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight;
	if(v < 0.3){ //local move
		if(v<0.1 && q>=2){ // attempt to swap Sq[m] and Sq[m+1]
			m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
			if(Sq[m]!=Sq[m+1]){
				oldE = Energies[m+1]; old_currD = currD;
				p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p; GetEnergies(); newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
				if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else {
					Sq[m+1] = Sq[m]; Sq[m] = p; currD = old_currD;
					UpdateWeightReplace(Energies[m+1],oldE); Energies[m+1] = oldE;
				}
			}
		} else if(v<0.2){ // attempt to delete Sq[m] and Sq[m+1]
			if(q>=2){
				m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
				if(Sq[m]==Sq[m+1]){
					oldE = Energies[m]; oldE2 = Energies[m+1]; old_currD = currD;
					memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
					for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
					GetEnergies(); newWeight = UpdateWeightDel(oldE,oldE2);
					if(val(rng) < Metropolis(newWeight)/Nop) currWeight = newWeight; else{
						q+=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
						currD = old_currD; UpdateWeightIns(oldE,oldE2);
					}
				}
			}
		} else if(q+2<qmax){ // attempt to insert Sq[m] and Sq[m+1]
			m = int(val(rng)*(q+1)); // m is between 0 and q
			memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
			old_currD = currD; p = diceNop(rng);
			for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = Sq[m+1] = p;
			GetEnergies(); newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
			if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else{
				q-=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
				currD = old_currD; d->RemoveElement(); d->RemoveElement();
			}
		} else qmax_achieved = 1;
	} else if(v < 0.5){ // attempting a fundamental cycle completion
		int oldq, j = 0, inv_pr; double wfactor;
		u = geometric_int(rng); // a random integer u is picked according to geometric distribution
		if(q >= u+rmin){
			inv_pr = min(rmax,q-u)-(rmin)+1;
			r = int(val(rng)*inv_pr) + (rmin);  // r is random integer between rmin and min(rmax,q-u)
			PickSubsequence(r+u); // indexes of the subsequence are min_index, min_index+1,..., max_index=min_index+r+u-1
			std::shuffle(Sq_subseq,Sq_subseq+r+u,rng);
			if(NoRepetitionCheck(Sq_subseq,r)){
				for(i=0;i<u;i++) Sq_gaps[i] = Sq_subseq[i+r];
				m = FindCycles(r);
				if(found_cycles > 0){ // cycles[m] is one of the found cycles, containing all the operators of Sq_subseq
					P = cycles[m]; for(i=0;i<r;i++) P.reset(Sq_subseq[i]);
					p = P.count(); // here, p is length of the complement sequence S'
					if(q+p-r < qmax){
						memcpy(Sq_backup,Sq,q*sizeof(int)); oldq = q;
						if(r<p)	     for(i=q-1;i>max_index;i--) Sq[i+p-r] = Sq[i]; // shift the values to the right
						else if(r>p) for(i=max_index+1;i<q;i++) Sq[i+p-r] = Sq[i]; // shift the values to the left
						for(i=0;i<p;i++){ while(!P.test(j)) j++; Sq_subseq[i] = Sq[min_index+i] = j++;}
						for(i=0;i<u;i++) Sq[min_index+p+i] = Sq_gaps[i]; // S' contains the remaining operators
						std::shuffle(Sq+min_index,Sq+min_index+p+u,rng);
						q += p-r; // the length q may have changed
						newWeight = UpdateWeight();
						wfactor = found_cycles; FindCycles(p); wfactor /= found_cycles;
						wfactor *= factorial[p]/factorial[r];
						wfactor *= inv_pr; inv_pr = min(rmax,q-u)-(rmin)+1; wfactor /= inv_pr;
						if(val(rng) < Metropolis(newWeight*wfactor)){
							currWeight = newWeight; cycles_used[m] = 1;
						} else{
							q = oldq; memcpy(Sq,Sq_backup,q*sizeof(int)); currWeight = UpdateWeight();
						}
					} else qmax_achieved = 1;
				}
			}
		}
	} else if(v < 0.9 && q>=2){ // attempting a block swap
		m = q==2 ? 0 : int(val(rng)*(q-1)); // m is between 0 and (q-2)
		oldE = currEnergy; for(i=0;i<=m;i++) ApplyOperator(Sq[i]);
		for(i=0;i<=m;i++)  Sq_backup[q-1-m+i] = Sq[i];
		for(i=m+1;i<q;i++) Sq_backup[i-m-1] = Sq[i];
		for(i=0;i<q;i++) { p = Sq[i]; Sq[i] = Sq_backup[i]; Sq_backup[i] = p;}
		memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
		GetEnergies(); newWeight = UpdateWeightReplace(oldE,currEnergy);
		if(val(rng) < Metropolis(newWeight)){
			z = lattice; currWeight = newWeight;
		} else{ UpdateWeightReplace(currEnergy,oldE); currEnergy = oldE;
			lattice = z; memcpy(Sq,Sq_backup,q*sizeof(int));
			memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
		}
	} else{ // flip of a random spin
		p = diceN(rng);	lattice.flip(p); newWeight = UpdateWeight();
		if(val(rng) < Metropolis(newWeight)) { z = lattice; currWeight = newWeight;}
			else { lattice.flip(p); currWeight = UpdateWeight();}
	}
}

double meanq = 0;
double maxq = 0;

#ifdef MEASURE_CUSTOM_OBSERVABLES

std::complex<double> calc_MD0(int n){ // calculate <z | MD_0 | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD0_size[n];i++) sum -= double(2*(int((MD0_product[n][i] & (~lattice)).count())%2)-1) * MD0_coeff[n][i];
	return sum;
}

void calc_MD0_trace(int n){
	std::bitset<N> og_lattice = lattice;
	currMD0 = currMD0_trace[0] = calc_MD0(n);
	for(int i=0;i<q;i++){
		ApplyOperator(Sq[i]);
		currMD0_trace[i+1] = calc_MD0(n);
		currMD0 *= currMD0_trace[i+1];
	}
	lattice = og_lattice;
}

std::complex<double> calc_MD(int n, int k){ // calculate d_k = <z | MD_k | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD_size[n][k];i++) sum -= double(2*(int((MD_product[n][k][i] & (~lattice)).count())%2)-1) * MD_coeff[n][k][i];
	return sum;
}

void calc_MD_trace_k(int n, int k){
	std::bitset<N> og_lattice = lattice;
	currMDk = currMDk_trace[0] = 1;
	for(int i=0;i<q;i++){
		ApplyOperator(Sq[i]);
		currMDk_trace[i+1] = calc_MD(n, k);
		currMDk *= currMDk_trace[i+1];
	}
	lattice = og_lattice;
}

void calc_MD_trace_l(int n, int l){
	currMDl = currMDl_trace[0] = 1;
	for(int i=0;i<q;i++){
		ApplyOperator(Sq[i]);
		currMDl_trace[i+1] = calc_MD(n, l);
		currMDl *= currMDl_trace[i+1];
	}
}

#endif

double measure_H(){
	double R = d->z[q]/(-beta); //z[q] = (-beta) * E_{z_q}
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_H2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta));
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Hdiag(){
	return currEnergy;
}

double measure_Hdiag2(){
	return currEnergy*currEnergy;
}

double measure_Hoffdiag(){
	double R = 0;
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_Hoffdiag2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta)) + currEnergy*(currEnergy - 2*measure_H());
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Hdiag_corr(){
  double R = (d->z[0]/(-beta))/(d->divdiffs[q]*beta_pow_factorial[q]).get_double();
  ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
  // init ds2 with full divided difference
  double tot_num = 0;
  for(int j = q; j >= 0; j--) ds2->AddElement((tau - beta)*(d->z[j]/(-beta)));
  for(int j = 0; j <= q; j++) {
	  ds1->AddElement(-tau*(d->z[j]/(-beta)));
	  tot_num += (d->z[j]/(-beta))*((ds1->divdiffs[j]*tau_pow_factorial[j])*(ds2->divdiffs[q-j]*tau_minus_beta_pow_factorial[q-j])).get_double();
	  ds2->RemoveElement();
  }
  R *= tot_num;
  return R;
}

double measure_Hoffdiag_corr(){
  double R = measure_H2();
  R += measure_Hdiag_corr();
  R -= 2.0 * (d->z[0]/(-beta)) * measure_H();
  return R;
}

double measure_Hdiag_int1(){
  double R = -(d->z[0]/(-beta))/(d->divdiffs[q]*beta_pow_factorial[q]).get_double();
  double tot_numerator = 0;
  for(int j = 0; j <= q; j++) {
    d->AddElement(-beta*(d->z[j]/(-beta)));
    tot_numerator += (d->z[j]/(-beta))*(d->divdiffs[q+1]*beta_pow_factorial[q+1]).get_double();
    d->RemoveElement();
  }
  R *= tot_numerator;
  return R;
}

double measure_Hoffdiag_int1(){
  double R = beta * measure_H2();
  R += measure_Hdiag_int1();
  R -= 2.0 * beta * (d->z[0]/(-beta)) * measure_H();
  return R;
}

double measure_Hdiag_int2(){
	double R = 0.0;
  	double prefac = -(d->z[0]/(-beta))/(d->divdiffs[q]*beta_pow_factorial[q]).get_double();
  	double Ezj, Ezu;
 	for(int j = 0; j <= q; j++) {
    	Ezj = d->z[j]/(-beta);
    	d->AddElement(-beta*Ezj);
	 	R += prefac*Ezj*beta*(d->divdiffs[q+1]*beta_pow_factorial[q+1]).get_double();
	  	d->AddElement(0);
	  	R += prefac*Ezj*(j+1.0)*(d->divdiffs[q+2]*beta_pow_factorial[q+2]).get_double();
    	for(int u = 0; u <= j; u++){
      		Ezu = d->z[u]/(-beta);
      		d->AddElement(-beta*Ezu);
      		R += prefac*Ezj*Ezu*(d->divdiffs[q+3]*beta_pow_factorial[q+3]).get_double();
      		d->RemoveElement();
    	}
    	d->RemoveElement();
	  	d->RemoveElement();
	}
  	return R;
}

double measure_Hoffdiag_int2(){
  double R = (beta * beta * measure_H2()) / 2.0;
  R += measure_Hdiag_int2();
  R -= beta * beta * (d->z[0]/(-beta)) * measure_H();
  return R;
}

double measure_Hdiag_int3(){
  double R = -(d->z[0]/(-beta))/(d->divdiffs[q]*beta_pow_factorial[q]).get_double();
  double tot_numerator = 0;
  double factor = 0;
  for(int j = 0; j <= q; j++) {
	  // reset and init divdiffs
	  ds1->CurrentLength=0; ds2->CurrentLength=0;
	  for(int k = q; k >= j; k--) ds2->AddElement((-beta/2)*(d->z[k]/(-beta)));
	  ds2->AddElement((-beta/2)*(d->z[j]/(-beta)));
	  for(int k = j-1; k >= 0; k--) ds2->AddElement((-beta/2)*(d->z[k]/(-beta)));
	  // sum over r to accumulate contribution for fixed j
	  for(int r = 0; r <= j; r++) {
		  ds1->AddElement((-beta/2)*(d->z[r]/(-beta)));
		  factor = (d->z[j]/(-beta))*(ds1->divdiffs[r]*beta_div2_pow_factorial[r]).get_double();
		  tot_numerator += (factor*(beta/2)*(ds2->divdiffs[q+1-r]*beta_div2_pow_factorial[q+1-r]).get_double());
		  ds2->AddElement(0);
		  tot_numerator += (factor*(j-r+1)*(ds2->divdiffs[q+2-r]*beta_div2_pow_factorial[q+2-r]).get_double());
		  for(int u = r; u <= j; u++){
			  ds2->AddElement((-beta/2)*(d->z[u]/(-beta)));
			  tot_numerator += (factor*(d->z[u]/(-beta))*(ds2->divdiffs[q+3-r]*beta_div2_pow_factorial[q+3-r]).get_double());
			  ds2->RemoveElement();
		  }
		  ds2->RemoveElement();
		  ds2->RemoveElement();
	  }
  }
  R *= tot_numerator;
  return R;
}

double measure_Hoffdiag_int3(){
  double R = (beta * beta * measure_H2()) / 8.0;
  R += measure_Hdiag_int3();
  R -= (beta * beta * (d->z[0]/(-beta)) * measure_H()) / 4.0;
  return R;
}

double measure_Z_magnetization(){
	return (2.0*lattice.count() - N)/N;
}

std::string name_of_observable(int n){
	std::string s;
	if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
		s = Mnames[n];
#endif
	} else switch(n-Nobservables){
			case 0: s = "H";             break;
			case 1: s = "H^2";           break;
			case 2: s = "H_{diag}";      break;
			case 3: s = "H_{diag}^2";    break;
			case 4: s = "H_{offdiag}";   break;
			case 5: s = "H_{offdiag}^2"; break;
			case 6: s = "Z_magnetization"; break;
			case 7: s = "measure_Hdiag_corr"; break;
      		case 8: s = "measure_Hdiag_int1"; break;
			case 9: s = "measure_Hdiag_int2"; break;
			case 10: s = "measure_Hdiag_int3"; break; 
			case 11: s = "measure_Hoffdiag_corr"; break;
			case 12: s = "measure_Hoffdiag_int1"; break;
			case 13: s = "measure_Hoffdiag_int2"; break; 
			case 14: s = "measure_Hoffdiag_int3"; break;  
	}
	return s;
}

double measure_observable(int n){
	double R = 0;
	if(valid_observable[n]) if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
		int i,k,len,cont,l,lenk,lenl,j,t,r;
		std::complex<double> T = 0;
		std::complex<double> A0z = calc_MD0(n);
		std::complex<double> Akz = 0.0;
		std::complex<double> prefac = 0.0;
		std::complex<double> inner_prefac = 0.0;
		switch(n){
			// measures <O>
			case 0:
				T += A0z; // add contribution of diagonal of custom O
				for(k=0;k<MNop[n];k++){ // sum over D_j P_j of O (inside loop is eq. 31)
					// P says which P_j in H are involved or not (as bitset)
					// len is k above eq. 30 or eq. 31
					// if k < q, cannot have \delta_tilde{P} as below eq. 31
					P = MP[n][k]; len = P.count(); if(len>q) continue; 
					// cannot have repeats in Sq[q-k+1],Sq[q-k+2],...,Sq[q] or cannot be making marked P
					// because indexing from 0, actually checking Sq[q-k],Sq[q-k+1],...,Sq[q-1]
					if(!NoRepetitionCheck(Sq+(q-len),len)) continue;
					// Sq[q-k] * Sq[q-k+1] * ... * Sq[q-1] == P (i.e. each Sq[q-1-i] is in P is what test is doing)
					cont = 0; for(i=0;i<len;i++) if(!P.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					// additional factorial[len] corrects that we include all orderings (formula is for specific ordering)
					T +=	(d->divdiffs[q-len]/d->divdiffs[q]).get_double() *
							(beta_pow_factorial[q-len]/beta_pow_factorial[q]).get_double()/factorial[len] *
						(currD_partial[q-len]/currD) * calc_MD(n,k); // calc_MD(n,k) gives \tilde{A}(z)
						// currD_parital[q-len]/currD is be ratio in eq.32
				}
				R = (currD*T).real()/currD.real(); // we importance-sample Re(W_C A_C)/Re(W_C)
				break;
			// measures <O^2>
			case 4:
				T += A0z*A0z; // diag/diag component
				// compute (k=0, l) contribution
				for(l=0;l<MNop[n];l++){
					calc_MD_trace_l(n ,l);
					Pl = MP[n][l]; lenl = Pl.count(); if(lenl>q) continue; // need lenl > q
					if(!NoRepetitionCheck(Sq+(q-lenl),lenl)) continue; // check {Sq[q-lenl], ..., Sq[q-1]} no repeats
					// checks that Sq[q-1-i] \in P_l for i=0, i=lenl-1.
					cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					T += (A0z * currMDl_trace[q]) * (d->divdiffs[q-lenl]/d->divdiffs[q]).get_double() *
						(beta_pow_factorial[q-lenl]/beta_pow_factorial[q]).get_double()/(factorial[lenl]) *
						(currD_partial[q-lenl]/currD);
				}
				for(k=0;k<MNop[n];k++){ 
					// compute (k, l=0) contribution
					Akz = calc_MD(n,k);
					Pk = MP[n][k]; lenk = Pk.count(); if(lenk>q) continue; 
					if(!NoRepetitionCheck(Sq+(q-lenk),lenk)) continue;
					cont = 0; for(i=0;i<lenk;i++) if(!Pk.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					T += (A0z * Akz) * (d->divdiffs[q-lenk]/d->divdiffs[q]).get_double() *
						(beta_pow_factorial[q-lenk]/beta_pow_factorial[q]).get_double()/factorial[lenk] *
						(currD_partial[q-lenk]/currD);
					// compute (k, l) contributions
					for(l=0;l<MNop[n];l++){
						calc_MD_trace_l(n ,l);
						Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue; // need lenk + lenl > q
						if(!NoRepetitionCheck(Sq+(q-lenk-lenl),lenl)) continue; // check {Sq[q-lenk-lenl], ..., Sq[q-1-lenk]} no repeats
						// checks that Sq[q-1-lenk-i] \in P_l for i=0, i=lenl-1.
						cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[q-1-lenk-i])){ cont = 1; break;} if(cont) continue;
						T += (Akz * currMDl_trace[q-lenk]) * (d->divdiffs[q-lenk-lenl]/d->divdiffs[q]).get_double() *
							(beta_pow_factorial[q-lenk-lenl]/beta_pow_factorial[q]).get_double()/(factorial[lenk]*factorial[lenl]) *
							(currD_partial[q-lenk-lenl]/currD);
					}
				}
				R = (currD*T).real()/currD.real();
				break;
			// measures <O(\tau)O>
			case 5:
				// add (k=0, l=0) pure diagonal contribution with Leibniz rule (like Hdiag_corr)
				calc_MD0_trace(n);
				prefac = A0z / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
				ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
				for(int j = q; j >= 0; j--) ds2->AddElement((tau - beta)*(d->z[j]/(-beta))); // init ds2
				for(int j = 0; j <= q; j++) {
	  				ds1->AddElement(-tau*(d->z[j]/(-beta)));
	  				T += prefac*currMD0_trace[j]*((ds1->divdiffs[j]*tau_pow_factorial[j])*(ds2->divdiffs[q-j]*tau_minus_beta_pow_factorial[q-j])).get_double();
	  				ds2->RemoveElement();
  				}
				// add (k=0, l) cross term contribution
				prefac = A0z / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
				for(l=0;l<MNop[n];l++){
					// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
					Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
					calc_MD_trace_l(n ,l);
					ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
					for(int j = q; j >= lenl; j--) ds2->AddElement((-tau)*(d->z[j]/(-beta))); // init ds2
					ds2->AddElement(0); // add dummy value for below continue logic to hold
					for (int j = lenl; j <= q; j++){
						ds2->RemoveElement();
						ds1->AddElement((-beta+tau)*(d->z[j-lenl]/(-beta)));
						// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
						if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
						cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
						// add contribution if relevant
						T += prefac*currMDl_trace[j]*(ds1->divdiffs[j-lenl]*tau_minus_beta_pow_factorial[j-lenl]).get_double() * 
						(ds2->divdiffs[q-j]*tau_pow_factorial[q-j]).get_double()/(factorial[lenl]) * 
						(currD_partial[j-lenl]/currD_partial[j]);
					}
				}
				// add (k, l=0) and (k, l) contributions
				for(k=0;k<MNop[n];k++){ 
					Akz = calc_MD(n,k);
					// \delta(\tilde{P}_k) conditions
					Pk = MP[n][k]; lenk = Pk.count(); if(lenk>q) continue; 
					if(!NoRepetitionCheck(Sq+(q-lenk),lenk)) continue;
					cont = 0; for(i=0;i<lenk;i++) if(!Pk.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					// compute (k, l=0) contribution
					prefac = Akz / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
					ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
					for(int j = q-lenk; j >= 0; j--) ds2->AddElement((-tau)*(d->z[j]/(-beta))); // init ds2
					for(int j = 0; j <= q-lenk; j++){
						ds1->AddElement((-beta+tau)*(d->z[j]/(-beta)));
						T += prefac*currMD0_trace[j]*(ds1->divdiffs[j]*tau_minus_beta_pow_factorial[j]).get_double() * 
						(ds2->divdiffs[q-lenk-j]*tau_pow_factorial[q-lenk-j]).get_double()/(factorial[lenk]) * 
						(currD_partial[q-lenk]/currD);
						ds2->RemoveElement();
					}
					for(l=0;l<MNop[n];l++){
						// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
						Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
						calc_MD_trace_l(n ,l);
						// compute (k, l) Lebniz rule contribution
						std::complex<double> prefac = Akz / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
						ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
						for(int j = q-lenk; j >= lenl; j--) ds2->AddElement((-tau)*(d->z[j]/(-beta))); // init ds2
						ds2->AddElement(0);
						for (int j = lenl; j <= q-lenk; j++){
							ds2->RemoveElement();
							ds1->AddElement((-beta+tau)*(d->z[j-lenl]/(-beta)));
							// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
							if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
							cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
							// add contribution if relevant
							T += prefac*currMDl_trace[j]*(ds1->divdiffs[j-lenl]*tau_minus_beta_pow_factorial[j-lenl]).get_double() * 
							(ds2->divdiffs[q-lenk-j]*tau_pow_factorial[q-lenk-j]).get_double()/(factorial[lenk]*factorial[lenl]) * 
							(currD_partial[j-lenl]/currD) * (currD_partial[q-lenk]/currD_partial[j]);
						}
					}

				}
				R = (currD*T).real()/currD.real();
				break;
			// measures \int_0^\beta \dtau <O(\tau)O>
			case 1:
				// add (k=0, l=0) pure diagonal contribution with Leibniz rule (like Hdiag_int1)
				calc_MD0_trace(n);
				prefac = -A0z / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
				for(int j = 0; j <= q; j++) {
	  				d->AddElement(-beta*(d->z[j]/(-beta)));
	  				T += prefac*currMD0_trace[j]*(d->divdiffs[q+1]*beta_pow_factorial[q+1]).get_double();
	  				d->RemoveElement();
  				}
				// add (k=0, l) cross term contribution
				prefac = -A0z / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
				for(l=0;l<MNop[n];l++){
					// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
					Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
					calc_MD_trace_l(n ,l);
					ds1->CurrentLength=0; // reset scratch divdiffs
					for(int j = q; j >= lenl; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
					ds1->AddElement(0); // add dummy value for continue logic
					for(int j = lenl; j <= q; j++){
						ds1->RemoveElement(); // remove before check so things work with continue
						// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
						if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
						cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
						// add \Lambda_j elements to \Theta_j multiset
						for(int i = 0; i <= j-lenl; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
						// add contribution if relevant
						T += prefac*currMDl_trace[j]*(ds1->divdiffs[q-lenl+1]*beta_pow_factorial[q-lenl+1]).get_double() * 
						(currD_partial[j-lenl]/currD_partial[j])/(factorial[lenl]);
						// remove \Lambda_j elements from \Theta_j multiset
						for(int i = 0; i <= j-lenl; i++) ds1->RemoveElement();
						ds1->RemoveElement();
					}
				}
				// add (k, l=0) and (k, l) contributions
				for(k=0;k<MNop[n];k++){ 
					Akz = calc_MD(n,k);
					// \delta(\tilde{P}_k) conditions
					Pk = MP[n][k]; lenk = Pk.count(); if(lenk>q) continue; 
					if(!NoRepetitionCheck(Sq+(q-lenk),lenk)) continue;
					cont = 0; for(i=0;i<lenk;i++) if(!Pk.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					// compute (k, l=0) contribution
					prefac = -Akz / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
					ds1->CurrentLength=0; // reset scratch divdiffs
					for(int j = q-lenk; j >= 0; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
					for(int j = 0; j <= q-lenk; j++){
						// add \Lambda_j elements to \Theta_j multiset
						for(int i = 0; i <= j; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
						T += prefac*currMD0_trace[j]*(ds1->divdiffs[q-lenk+1]*beta_pow_factorial[q-lenk+1]).get_double() * 
						(currD_partial[q-lenk]/currD)/(factorial[lenk]);
						// remove \Lambda_j elements from \Theta_j multiset
						for(int i = 0; i <= j; i++) ds1->RemoveElement();
						ds1->RemoveElement();
					}
					for(l=0;l<MNop[n];l++){
						// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
						Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
						calc_MD_trace_l(n ,l);
						// compute (k, l) Lebniz rule contribution
						prefac = -Akz / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
						ds1->CurrentLength=0; // reset scratch divdiffs
						for(int j = q-lenk; j >= lenl; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
						ds1->AddElement(0); // add dummy value for continue logic
						for (int j = lenl; j <= q-lenk; j++){
							ds1->RemoveElement(); // remove before check so things work with continue
							// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
							if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
							cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
							// add \Lambda_j elements to \Theta_j multiset
							for(int i = 0; i <= j-lenl; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
							// add contribution if relevant
							T += prefac*currMDl_trace[j]*(ds1->divdiffs[q-lenk-lenl+1]*beta_pow_factorial[q-lenk-lenl+1]).get_double() * 
							(currD_partial[j-lenl]/currD) * (currD_partial[q-lenk]/currD_partial[j])/(factorial[lenk]*factorial[lenl]);
							// remove \Lambda_j elements from \Theta_j multiset
							for(int i = 0; i <= j-lenl; i++) ds1->RemoveElement();
						}
					}

				}
				R = (currD*T).real()/currD.real();
				break;
			// measures \int_0^\beta \dtau \tau <O(\tau)O>
			case 2: 
				// add (k=0, l=0) pure diagonal contribution with Leibniz rule (like Hdiag_int2)
				calc_MD0_trace(n);
				prefac = -A0z / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
				for(int j = 0; j <= q; j++) {
	  				d->AddElement(-beta*(d->z[j]/(-beta)));
	  				T += prefac*currMD0_trace[j]*beta*(d->divdiffs[q+1]*beta_pow_factorial[q+1]).get_double();
					d->AddElement(0);
					T += prefac*currMD0_trace[j]*(j+1.0)*(d->divdiffs[q+2]*beta_pow_factorial[q+2]).get_double();
					for(int u = 0; u <= j; u++){
						d->AddElement(-beta*(d->z[u]/(-beta)));
						T += prefac*currMD0_trace[j]*(d->z[u]/(-beta))*(d->divdiffs[q+3]*beta_pow_factorial[q+3]).get_double();
						d->RemoveElement();
					}
	  				d->RemoveElement();
					d->RemoveElement();
  				}
				// add (k=0, l) cross term contribution
				prefac = -A0z / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
				for(l=0;l<MNop[n];l++){
					// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
					Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
					calc_MD_trace_l(n ,l);
					ds1->CurrentLength=0; // reset scratch divdiffs
					for(int j = q; j >= lenl; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
					ds1->AddElement(0); // add dummy value for continue logic
					for(int j = lenl; j <= q; j++){
						ds1->RemoveElement(); // remove before check so things work with continue
						// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
						if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
						cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
						// add \Lambda_j elements to \Theta_j multiset
						for(int i = 0; i <= j-lenl; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
						// add contribution if relevant
						T += prefac*currMDl_trace[j]*beta*(ds1->divdiffs[q-lenl+1]*beta_pow_factorial[q-lenl+1]).get_double() * 
						(currD_partial[j-lenl]/currD_partial[j])/(factorial[lenl]);
						ds1->AddElement(0);
						T += prefac*currMDl_trace[j]*(j-lenl+1.0)*(ds1->divdiffs[q-lenl+2]*beta_pow_factorial[q-lenl+2]).get_double() * 
						(currD_partial[j-lenl]/currD_partial[j])/(factorial[lenl]);
						for(int u = 0; u <= j-lenl; u++){
							ds1->AddElement(-beta*(d->z[u]/(-beta)));
							T += prefac*currMDl_trace[j]*(d->z[u]/(-beta))*(ds1->divdiffs[q-lenl+3]*beta_pow_factorial[q-lenl+3]).get_double() * 
							(currD_partial[j-lenl]/currD_partial[j])/(factorial[lenl]);
							ds1->RemoveElement();
						}
						// remove \Lambda_j elements from \Theta_j multiset
						for(int i = 0; i <= j-lenl; i++) ds1->RemoveElement();
						ds1->RemoveElement();
					}
				}
				// add (k, l=0) and (k, l) contributions
				for(k=0;k<MNop[n];k++){ 
					Akz = calc_MD(n,k);
					// \delta(\tilde{P}_k) conditions
					Pk = MP[n][k]; lenk = Pk.count(); if(lenk>q) continue; 
					if(!NoRepetitionCheck(Sq+(q-lenk),lenk)) continue;
					cont = 0; for(i=0;i<lenk;i++) if(!Pk.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					// compute (k, l=0) contribution
					prefac = -Akz / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
					ds1->CurrentLength=0; // reset scratch divdiffs
					for(int j = q-lenk; j >= 0; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
					for(int j = 0; j <= q-lenk; j++){
						// add \Lambda_j elements to \Theta_j multiset
						for(int i = 0; i <= j; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
						T += prefac*currMD0_trace[j]*beta*(ds1->divdiffs[q-lenk+1]*beta_pow_factorial[q-lenk+1]).get_double() * 
						(currD_partial[q-lenk]/currD)/(factorial[lenk]);
						ds1->AddElement(0);
						T += prefac*currMD0_trace[j]*(j+1.0)*(ds1->divdiffs[q-lenk+2]*beta_pow_factorial[q-lenk+2]).get_double() * 
						(currD_partial[q-lenk]/currD)/(factorial[lenk]);
						for(int u = 0; u <= j; u++){
							ds1->AddElement(-beta*(d->z[u]/(-beta)));
							T += prefac*currMD0_trace[j]*(d->z[u]/(-beta))*(ds1->divdiffs[q-lenk+3]*beta_pow_factorial[q-lenk+3]).get_double() * 
							(currD_partial[q-lenk]/currD)/(factorial[lenk]);
							ds1->RemoveElement();
						}
						// remove \Lambda_j elements from \Theta_j multiset
						ds1->RemoveElement(); // Remove 0
						for(int i = 0; i <= j; i++) ds1->RemoveElement();
						ds1->RemoveElement(); // Remove from e{-tau[]} multiset
					}
					for(l=0;l<MNop[n];l++){
						// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
						Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
						calc_MD_trace_l(n ,l);
						// compute (k, l) Lebniz rule contribution
						prefac = -Akz / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
						ds1->CurrentLength=0; // reset scratch divdiffs
						for(int j = q-lenk; j >= lenl; j--) ds1->AddElement((-beta)*(d->z[j]/(-beta))); // init ds1
						ds1->AddElement(0); // add dummy value for continue logic
						for (int j = lenl; j <= q-lenk; j++){
							ds1->RemoveElement(); // remove before check so things work with continue
							// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
							if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
							cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
							// add \Lambda_j elements to \Theta_j multiset
							for(int i = 0; i <= j-lenl; i++) ds1->AddElement((-beta)*(d->z[i]/(-beta)));
							// add contribution if relevant
							T += prefac*currMDl_trace[j]*beta*(ds1->divdiffs[q-lenk-lenl+1]*beta_pow_factorial[q-lenk-lenl+1]).get_double() * 
							(currD_partial[j-lenl]/currD) * (currD_partial[q-lenk]/currD_partial[j])/(factorial[lenk]*factorial[lenl]);
							ds1->AddElement(0);
							T += prefac*currMDl_trace[j]*(j-lenl+1.0)*(ds1->divdiffs[q-lenk-lenl+2]*beta_pow_factorial[q-lenk-lenl+2]).get_double() * 
							(currD_partial[j-lenl]/currD) * (currD_partial[q-lenk]/currD_partial[j])/(factorial[lenk]*factorial[lenl]);
							for(int u = 0; u <= j-lenl; u++){
								ds1->AddElement(-beta*(d->z[u]/(-beta)));
								T += prefac*currMDl_trace[j]*(d->z[u]/(-beta))*(ds1->divdiffs[q-lenk-lenl+3]*beta_pow_factorial[q-lenk-lenl+3]).get_double() * 
								(currD_partial[j-lenl]/currD) * (currD_partial[q-lenk]/currD_partial[j])/(factorial[lenk]*factorial[lenl]);
								ds1->RemoveElement();
							}
							ds1->RemoveElement(); // remove 0
							// remove \Lambda_j elements from \Theta_j multiset
							for(int i = 0; i <= j-lenl; i++) ds1->RemoveElement();
						}
					}

				}
				R = (currD*T).real()/currD.real();
				break;
			// measures \int_0^{\beta/2} \dtau \tau <O(\tau)O>
			case 3:
				// add (k=0, l=0) pure diagonal contribution with Leibniz rule (like Hdiag_int3)
				calc_MD0_trace(n);
				prefac = -A0z / (d->divdiffs[q]*beta_pow_factorial[q]).get_double();
				for(int j = 0; j <= q; j++) {
					// reset and init divdiffs
					ds1->CurrentLength=0; ds2->CurrentLength=0;
					for(int i = q; i >= j; i--) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
					for(int i = j; i >= 0; i--) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
					// sum over r to accumulate contribution for fixed j
					for(int r = 0; r <= j; r++) {
						ds1->AddElement((-beta/2.0)*(d->z[r]/(-beta)));
						inner_prefac = prefac * currMD0_trace[j] * (ds1->divdiffs[r]*beta_div2_pow_factorial[r]).get_double();
						T += (beta/2.0)*inner_prefac*(ds2->divdiffs[q+1-r]*beta_div2_pow_factorial[q+1-r]).get_double();
						ds2->AddElement(0);
						T += (inner_prefac*(j-r+1.0)*(ds2->divdiffs[q+2-r]*beta_div2_pow_factorial[q+2-r]).get_double());
						for(int u = r; u <= j; u++){
							ds2->AddElement((-beta/2.0)*(d->z[u]/(-beta)));
							T += (d->z[u]/(-beta))*inner_prefac*(ds2->divdiffs[q+3-r]*beta_div2_pow_factorial[q+3-r]).get_double();
							ds2->RemoveElement();
						}
	  					ds2->RemoveElement();
						ds2->RemoveElement();
					}
  				}
				// add (k=0, l) cross term contribution
				prefac = -A0z / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
				for(l=0;l<MNop[n];l++){
					// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
					Pl = MP[n][l]; lenl = Pl.count(); if(lenl>q) continue;
					calc_MD_trace_l(n ,l);
					for (int j = lenl; j <= q; j++){
						// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
						if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
						cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
						ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
						// build e^{-\beta/2 [\Lambda_j]}
						for(int i = j; i <= q; i++) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
						// build e^{-\beta/2 [\Omega_{j,r}]} (aka add [\overline{\Theta_{j,r}}] to [\Lambda_j]
						for(int i = j-lenl; i>=0; i--) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
						// sum over r
						for(int r = 0; r <= j-lenl; r++){
							// build e^{-\beta/2 [\Theta_{j,r}]}
							//for(int i = 0; i <= r; i++) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
							ds1->AddElement((-beta/2)*(d->z[r]/(-beta)));
							inner_prefac = prefac*(currMDl_trace[j]/factorial[lenl])*(currD_partial[j-lenl]/currD_partial[j])*(ds1->divdiffs[r]*beta_div2_pow_factorial[r]).get_double();
							// add (beta/2)*e^{-beta/2[\Omega_{j,r}]} contribution
							T += inner_prefac*(beta/2)*(ds2->divdiffs[q+1-lenl-r]*beta_div2_pow_factorial[q+1-lenl-r]).get_double();
							// add N([\overline{\theta}_{j,r}]) e^{-beta/2[\Omega_{j,r}] + {0}} contribution
							ds2->AddElement(0);
							T += inner_prefac*(j-lenl-r+1.0)*(ds2->divdiffs[q+2-lenl-r]*beta_div2_pow_factorial[q+2-lenl-r]).get_double();
							// add \sum_u contribution
							for(int u = r; u<= j-lenl; u++){
								ds2->AddElement((-beta/2)*(d->z[u]/(-beta)));
								T += inner_prefac*(d->z[u]/(-beta))*(ds2->divdiffs[q+3-lenl-r]*beta_div2_pow_factorial[q+3-lenl-r]).get_double();
								ds2->RemoveElement();
							}
							ds2->RemoveElement(); // removes AddElement(0)
							ds2->RemoveElement(); // removes AddElement((-beta/2)*(d->z[r]/(-beta)));
						}
					}
				}
				// add (k, l=0) and (k, l) contributions
				for(k=0;k<MNop[n];k++){ 
					Akz = calc_MD(n,k);
					// \delta(\tilde{P}_k) conditions
					Pk = MP[n][k]; lenk = Pk.count(); if(lenk>q) continue; 
					if(!NoRepetitionCheck(Sq+(q-lenk),lenk)) continue;
					cont = 0; for(i=0;i<lenk;i++) if(!Pk.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
					// compute (k, l=0) contribution
					prefac = -(Akz / factorial[lenk]) / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
					for (int j = 0; j <= q-lenk; j++){
						ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
						// build e^{-\beta/2 [\Lambda_j]}
						for(int i = j; i <= q-lenk; i++) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
						// build e^{-\beta/2 [\Omega_{j,r}]} (aka add [\overline{\Theta_{j,r}}] to [\Lambda_j]
						for(int i = j; i>=0; i--) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
						// sum over r
						for(int r = 0; r <= j; r++){
							// build e^{-\beta/2 [\Theta_{j,r}]}
							//for(int i = 0; i <= r; i++) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
							ds1->AddElement((-beta/2)*(d->z[r]/(-beta)));
							inner_prefac = prefac*currMD0_trace[j]*(currD_partial[q-lenk]/currD)*(ds1->divdiffs[r]*beta_div2_pow_factorial[r]).get_double();
							// add (beta/2)*e^{-beta/2[\Omega_{j,r}]} contribution
							T += inner_prefac*(beta/2)*(ds2->divdiffs[q+1-lenk-r]*beta_div2_pow_factorial[q+1-lenk-r]).get_double();
							// add N([\overline{\theta}_{j,r}]) e^{-beta/2[\Omega_{j,r}] + {0}} contribution
							ds2->AddElement(0);
							T += inner_prefac*(j-r+1.0)*(ds2->divdiffs[q+2-lenk-r]*beta_div2_pow_factorial[q+2-lenk-r]).get_double();
							// add \sum_u contribution
							for(int u = r; u<= j; u++){
								ds2->AddElement((-beta/2)*(d->z[u]/(-beta)));
								T += inner_prefac*(d->z[u]/(-beta))*(ds2->divdiffs[q+3-lenk-r]*beta_div2_pow_factorial[q+3-lenk-r]).get_double();
								ds2->RemoveElement();
							}
							ds2->RemoveElement(); // removes AddElement(0)
							ds2->RemoveElement(); // removes AddElement((-beta/2)*(d->z[r]/(-beta)));
						}
					}
					// compute \sum_l (k, l) Leibniz rule contribution
					for(l=0;l<MNop[n];l++){
						// initial 1^{(j)}_{\tilde{P}_l}(S_{i_q}) check
						Pl = MP[n][l]; lenl = Pl.count(); if(lenk+lenl>q) continue;
						calc_MD_trace_l(n ,l);
						prefac = -(Akz / factorial[lenk]) / ((d->divdiffs[q]*beta_pow_factorial[q]).get_double());
						for (int j = lenl; j <= q-lenk; j++){
							// remainder of 1^{(j)}_{\tilde{P}_l}(S_{i_q}) condition
							if(!NoRepetitionCheck(Sq+(j-lenl),lenl)) continue;
							cont = 0; for(i=0;i<lenl;i++) if(!Pl.test(Sq[j-lenl+i])){ cont = 1; break;} if(cont) continue;
							ds1->CurrentLength=0; ds2->CurrentLength=0; // reset scratch divdiffs
							// build e^{-\beta/2 [\Lambda_j]}
							for(int i = j; i <= q-lenk; i++) ds2->AddElement((-beta/2.0)*(d->z[i]/(-beta)));
							// build e^{-\beta/2 [\Omega_{j,r}]} (aka add [\overline{\Theta_{j,r}}] to [\Lambda_j]
							for(int i = j-lenl; i>=0; i--) ds2->AddElement((-beta/2.0)*(d->z[i]/(-beta)));
							// sum over r
							for(int r = 0; r <= j-lenl; r++){
								// build e^{-\beta/2 [\Theta_{j,r}]}
								//for(int i = 0; i <= r; i++) ds2->AddElement((-beta/2)*(d->z[i]/(-beta)));
								ds1->AddElement((-beta/2.0)*(d->z[r]/(-beta)));
								inner_prefac = prefac*(currMDl_trace[j]/factorial[lenl])*(currD_partial[j-lenl]/currD)*(currD_partial[q-lenk]/currD_partial[j]) *
								(ds1->divdiffs[r]*beta_div2_pow_factorial[r]).get_double();
								// add (beta/2)*e^{-beta/2[\Omega_{j,r}]} contribution
								T += inner_prefac*(beta/2.0)*(ds2->divdiffs[q+1-lenl-lenk-r]*beta_div2_pow_factorial[q+1-lenl-lenk-r]).get_double();
								// add N([\overline{\theta}_{j,r}]) e^{-beta/2[\Omega_{j,r}] + {0}} contribution
								ds2->AddElement(0);
								T += inner_prefac*(j-lenl-r+1.0)*(ds2->divdiffs[q+2-lenl-lenk-r]*beta_div2_pow_factorial[q+2-lenl-lenk-r]).get_double();
								// add \sum_u contribution
								for(int u = r; u<= j-lenl; u++){
									ds2->AddElement((-beta/2.0)*(d->z[u]/(-beta)));
									T += inner_prefac*(d->z[u]/(-beta))*(ds2->divdiffs[q+3-lenl-lenk-r]*beta_div2_pow_factorial[q+3-lenl-lenk-r]).get_double();
									ds2->RemoveElement();
								}
								ds2->RemoveElement(); // removes AddElement(0)
								ds2->RemoveElement(); // removes AddElement((-beta/2)*(d->z[r]/(-beta)));
							}
						}
					}

				}
				R = (currD*T).real()/currD.real();
				break;
		}
#endif
	} else  switch(n-Nobservables){
			case 0:	R = measure_H(); break;
			case 1:	R = measure_H2(); break;
			case 2:	R = measure_Hdiag(); break;
			case 3:	R = measure_Hdiag2(); break;
			case 4:	R = measure_Hoffdiag(); break;
			case 5:	R = measure_Hoffdiag2(); break;
			case 6: R = measure_Z_magnetization(); break;
			case 7: R = measure_Hdiag_corr(); break;
            case 8: R = measure_Hdiag_int1(); break;
			case 9: R = measure_Hdiag_int2(); break;
			case 10: R = measure_Hdiag_int3(); break;
			case 11: R = measure_Hoffdiag_corr(); break; 
			case 12: R = measure_Hoffdiag_int1(); break;
			case 13: R = measure_Hoffdiag_int2(); break;
			case 14: R = measure_Hoffdiag_int3(); break;
	}
	return R;
}

void measure(){
	double R, sgn; int i;
	currWeight = GetWeight(); sgn = currWeight.sgn();
	meanq += q; if(maxq < q) maxq = q; in_bin_sum_sgn += sgn;
	if((measurement_step+1) % bin_length == 0){
		in_bin_sum_sgn /= bin_length; bin_mean_sgn[measurement_step/bin_length] = in_bin_sum_sgn; in_bin_sum_sgn = 0;
	}
	for(i=0;i<N_all_observables;i++){
		R = measure_observable(i); in_bin_sum[i] += R*sgn;
		if((measurement_step+1) % bin_length == 0){
			in_bin_sum[i] /= bin_length; bin_mean[i][measurement_step/bin_length] = in_bin_sum[i]; in_bin_sum[i] = 0;
		}
	}
}
