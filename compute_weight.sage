# MatRiCT+ invertibility bound heuristic model computation (https://gitlab.com/raykzhao/matrict_plus)
# Ron Steinfeld, 29 July 2021
# Adapted by authors of the submission to encompass larger infinity norm bounds (gamma), 7 February 2024
# Please cite the MatRiCT+ paper accepted at IEEE S&P 2022 conference (full version available at https://eprint.iacr.org/2021/545) when using scripts

''' 
SageMath script to compute the weight 'w' and infinity norm bound 'gamma' of the different challenge sets as presented in Table 3 in Section D.4 of the submission.
The weight and the inifinity norm bound directly define the bound on the square of the l2-norm (tau in python script and T_2 in paper) and the bound on the operator norm (T in python script and T_op in paper) and are used in 'proof_size_estimate.py'.
To run the program, simply run the SageMath script, for example >>> load("compute_weight.sage")
'''

import math
import numpy as np
import time
import csv
from sage.modules.free_module_integer import IntegerLattice

# maximal number of witnesses in each iteration of LaBRADOR (cf. Equation 9 in submission)
# max_r = 6 ceil(sqrt(N)) + 1 for N <= 10~000, r_max <= 601
max_r = 601

# simple function to compute a prime number with the right splitting behaviour 
def compute_prime(num_bit,ell):
	for i in range(2^(num_bit-1),2^(num_bit+3)):
		if i % (4*ell) == (2*ell+1) and is_prime(i):
			return(i)
	return 0

# compute gama and weight for the almost fully splitting case	
def compute_weight_almost_fully_splitting(d,ell,q,secparam):
	max_gam = 100 # maximal gamma we are interested in;
	gam = 1
	tilde_w = 1 
	delt = d//ell # degree of the irred factors of X^d+1 mod q (delta in paper)
	while True: 
		# FOLLOWING THE FORMULAS OF LEMMA D.1
		A = RDF(gamma((tilde_w+1)/2)/sqrt(pi()*(ell*gam)^tilde_w)) # tilde_w'th moment of Gaussian with var 1/(2*gam*ell) 
		eta = RDF(ell^tilde_w * factorial(ell-tilde_w) / factorial(ell)) # eta as in paper 
		expect_EM2 = eta * RDF(1/q + (1-1/q)*A) # Expected value of M2
		expect_B = expect_EM2^delt # expected bound for B (with delta independent coefficients per CRT slot)
		
		# FOLLOWING THE EQUATION 2 
		final_term = (5+2*ell)*max_r * expect_B # (5+2l)Br = one of the additive terms in the knowledge error of LaBRADOR 

		# logs of Results
		lgEM2 = round(RDF(log(expect_EM2)/log(2)),1)
		lgB = round(RDF(log(expect_B)/log(2)),1)
		lgfinal_term = round(RDF(log(final_term)/log(2)),1)
	
		# increase the weight if security level is not yet reached and weight is below the number of CRT slots	
		if lgfinal_term > -secparam and tilde_w < ell:
			tilde_w += 1
			continue
		# increase gamma, reset weight to 1 if security level is not yet reached but weight is reaching the number of CRT slots 
		if lgfinal_term > - secparam and tilde_w >= ell:
			if gam >= max_gam:
				print("ERROR: TOO LARGE GAMMA >= ",max_gam)
				break
			else: 
				tilde_w = 1
				gam += 1
				continue
		else:
			final_w = tilde_w*delt
			binom = binomial(d,final_w)*(2*gam)^(final_w)
			if log(binom,2) < secparam: # double check if size of challenge space is big enough
				print("ERROR, CHALLENGE SPACE TOO SMALL")
				break
			# stop once securtiy level is reached
			else:
				print("\n")
				print ("=== COMPUTE WEIGHT - ALMOST FULLY SPLITTING ===")
				print("\n")		
				print ("Input Pars:")
				print ("ring degree (d in paper) = ", d) 
				print ("no. of irred factors (r in paper) = ", ell)
				print ("bit length of modulus (q in paper) = ", math.ceil(log(q)/log(2)))
				print ("security paramater (lambda in paper)= ", secparam)
				print("\n")
				print ("Output Pars:")
				print ("log of M2 = ",lgEM2)
				print ("log of well-spreadness = ",lgB)
				print ("log of final term in knowledge error (5+2l)Br= ",lgfinal_term)
				print ("weight of each S_i set (tilde_w in paper)= ", tilde_w)
				print ("total weight of callenge elements (w in paper)= ", final_w)
				print ("infinity norm bound (gamma in paper) = ", gam)
				break	

# compute gama and weight for the two splitting case	
def compute_weight_two_splitting(d,secparam):
	w = 1
	gam = 1
	max_gam = 100 # maximal gamma we are interested in;
	while True:
		binom = binomial(d,w)*(2*gam)^w
		# in two-splitting: B = 1/|C| -> |C|/(5+2l)r has to be larger than 2^{secpar} 
		if log(binom/(9*max_r),2) < secparam: # only condition: size of challenge space is big enough; we divide by (5+2l)r where l=2 and r<= max_r
			if w < d-1:
				w += 1
				continue
			if w >= d-1:
				if gam >= max_gam:
					print("ERROR, MAX_GAM =", max_gam, " WAS REACHED")	
				else:
					w = 1
					gam += 1
				continue
		else:
			break
	print("\n")
	print ("=== COMPUTE WEIGHT - TWO SPLITTING ===")
	print("\n")		
	print ("Input Pars:")
	print ("ring degree (d in paper) = ", d) 
	print ("security paramater (lambda in paper)= ", secparam)
	print("\n")
	print ("Output Pars:")
	print("total weight of challenge elements (w in paper)= ",w)
	print ("infinity norm bound (gamma in paper) = ", gam)
	
	
#### Aggregating Falcon Signatures With LaBRADOR ####
#### Reproducing the numbers of Table 1 of the submission ######

### ALMOST FULLY SPLITTING

## aiming security level of 128-bits
secparam=128

# splitting up to level 4
d=128
ell=d//4
q=compute_prime(47,ell)
#q=compute_prime(63,ell) #gives slightly smaller parameters
compute_weight_almost_fully_splitting(d,ell,q,secparam)

## aiming security level of 256-bits
secparam=256

# splitting up to level 8
d =256
ell=d//8
q=compute_prime(47,ell)
#q=compute_prime(63,ell) #gives slightly smaller parameters
compute_weight_almost_fully_splitting(d,ell,q,secparam)

### TWO SPLITTING

## aiming security level of 128-bits
secparam=128
d=64
compute_weight_two_splitting(d,secparam)
	
## aiming security level of 256-bits
secparam=256
d=128
compute_weight_two_splitting(d,secparam)
