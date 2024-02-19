# aggregate-falcon
This repository contains supplemantary material for the paper "Aggregating Falcon Signatures With LaBRADOR".

It contains the following files:

- `jl.py`: his python script computes the constant 'C_1' and slack 'sqrt(lambda/C_2)' of the Johnson-Lindenstrauss lemma, Lemma 2.2 in the paper. 
The constant 'c' (C_1 in the paper) is used in 'proof_size_estimate.py' as JL_const and the slack (sqrt(lambda/C_2 in paper)) is used as JL_slack
To run the program, simply run the python script, for example >>> python3 ./jl.py

- `compute_weight.sage`: This SageMath script computes the weight 'w' and infinity norm bound 'gamma' of the different challenge sets as presented in Table 3 in Section D.4 of the paper. The weight and the inifinity norm bound directly define the bound on the square of the l2-norm (tau in python script and T_2 in paper) and the bound on the operator norm (T in python script and T_op in paper) and are used in 'proof_size_estimate.py'.
To run the program, simply run the SageMath script, for example >>> load("compute_weight.sage")

- `SIS_hardness.sage`:  This SageMath script computes the estimated hardness of Module-SIS needed for our aggregate signature. 
It is derived from corresponding code provided to us by Gregor Seiler.
Its output is then used as 'kappa_lim' in the python script 'proof_size_estimate.py' 
Concretely, it computes the corresponding 'kappa_lim' lists for the classes
	1) FALCON_64_128 (used for two-splitting 128-bit security)
	2) FALCON_128_128 (used for almost-fully-splititing 128-bit security
	3) FALCON_128_256 (used for two-splitting 256-bit security)
	4) FALCON_256_256 (used for almost-fully-splitting 256-bit security)
To run the program, simply run the SageMath script, for example >>> load("SIS_hardness.sage")

- `proof_size_estimate.py`: This python script computes the estimated sizes of our aggregate signature as presented in Section 6 and Appendix E of our paper. It is derived from corresponding code provided to us by Gregor Seiler.
It first provides the numbers for the comparision with [JRS23], Squirrel and Chipmunk.
Then it provides the numbers stored in 'estimates-lin-.csv' to derive the plots through 'plot_paper.py' later.
To run the program, simply run the python script, for example >>> python3 ./proof_size_estimate.py

- `plot_paper.py`: After having run 'proof_size_estimate.py', this python script computes the corresponding plots for Section 6.
To obtain Figure 3 (left) in the paper, run the python script with the only function not commented out being 'plot_512_1024_AS_lin()'
To obtain Figure 3 (right) in the paper, run the python script with the only function not commented out being 'plot_512_AS_lin_no_salt()'
To obtain Figure 4 in the paper, first run the python script with the only function not commented out being 'plot_512_AS_2S_FS_lin()', then run it again with the only function not commented out being 'plot_1024_AS_2S_FS_lin()'.

- `poly_arith/*`: Experiments with polynomial arithmetic for various choices of parameters. Check the `README.md` inside for instructions.

