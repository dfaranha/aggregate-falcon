''' 
This python script computes the constant 'C_1' and slack 'sqrt(lambda/C_2)' of the Johnson-Lindenstrauss lemma, Lemma 2.2 in the submission. 
The constant 'c' (C_1 in the paper) is used in 'proof_size_estimate.py' as JL_const and the slack (sqrt(lambda/C_2 in paper)) is used as JL_slack
To run the program, simply run the python script, for example >>> python3 ./jl.py
'''

from scipy.stats import chi2,norm
from scipy.special import binom
from math import log2,log,sqrt,ceil,floor,pi

def ceil_decimal(x, numdecials):
    if numdecials == 0:
        return ceil(x)
    factor = 10**numdecials
    return ceil(x*factor)/factor

def floor_decimal(x, numdecimals):
    if numdecimals == 0:
        return floor(x)
    factor = 10**numdecimals
    return floor(x*factor)/factor

#Added margin to get more sound proof for Labrador strengthening
def tailbounds_chl21(secparam, round=False, numdecimals=2,margin=0):
    target = 2.0**(-secparam-margin) 
    projdim = 2*secparam

    #Pr[|<r, w>| > alpha norm(w)] <= target * 1/projdim (for union bound over coords of proj)
    inv = norm.isf(target*(1/projdim))
    alpha = inv / sqrt(2)

    #Pr[norm(Pw)^2 < beta norm(w)^2] <= target
    inv = chi2.ppf(target, projdim)
    beta = inv / 2.0

    #Pr[norm(Pw)^2 > gamma norm(w)^2] <= target
    inv = chi2.isf(target, projdim)
    gamma = inv / 2.0

    if(round):
        alpha = ceil_decimal(alpha, numdecimals)
        beta = floor_decimal(beta, numdecimals)
        gamma = ceil_decimal(gamma, numdecimals)

    return projdim, alpha, beta, gamma

projdim128, alpha128, beta128, gamma128 = tailbounds_chl21(128)
assert alpha128 < 9.75
assert beta128 > 30
assert gamma128 < 337
projdim256, alpha256, beta256, gamma256 = tailbounds_chl21(256, margin=1)

#Proof of Cor 3.3 in [GHL21] works as long as sqrt(beta)b < q/(8d)
#Output: c s.t. the lemma holds for b < q/(cd)
def jl_glh21_normreq(beta, round=False, numdecimals=2):
    c = 8 * sqrt(beta)
    if round:
        return ceil_decimal(c, numdecimals)
    return c
    
assert jl_glh21_normreq(beta128) < 45

#Output: c s.t. the lemma holds for b < q/c
def jl_labrador_normreq(secparam, projdim, alpha, beta, round=False, numdecimals=2):
    #Case 1:
    c1 = sqrt(beta)/(1-alpha/10)

    #Case 2
    #Finding number of coordinates of proj (not mod q) that can have magnitude 
    #greater than q/120 before prob is < 2^-secparam 
    union_bound = 0
    last_val = union_bound
    bound = 2**(-secparam)
    prob_other_coords = 3.0**(-projdim)
    i = 0
    while union_bound < bound and i <= projdim:
        last_val = union_bound
        i += 1
        prob_other_coords *= 3
        union_bound += binom(projdim, i)*prob_other_coords
    case2numdigits = i

    c2 = 120*sqrt(beta) / sqrt(case2numdigits)
    c = max(c1, c2)

    #Case 3
    prob_qhalf = 2 * norm.sf(sqrt(2)*alpha/2)
    c3 = 11*sqrt(2*beta)/(sqrt(pi)*((1/sqrt(2))-prob_qhalf-0.3))
    c = max(c, c3)

    if round:
        return ceil_decimal(c, numdecimals)
    return c

assert jl_labrador_normreq(128, 256, alpha128, beta128) <= 125

c128 = jl_labrador_normreq(128, 256, alpha128, beta128)
c256 = jl_labrador_normreq(256, 512, alpha256, beta256)

def print_jl_msg(secparam):
    projdim, alpha, beta, _ = tailbounds_chl21(secparam) 
    c = jl_labrador_normreq(secparam, projdim, alpha, beta, round=True, numdecimals=0) 
    beta_rounded = floor(beta)
    print(f"=== Modular Johnson-Lindenstrauss for secparam = {secparam} ===")
    print("Let w be a fixed vector in Z^d.")
    print(f"Let {projdim} be the dimension of our projections.")
    print(f"Let b <= q/{c}.") #c corresponds to C_1 in Lemma 2.2 in paper
    print()
    print(f"If norm(w)_2 > b, then the probability that norm(proj)_2 < sqrt({beta_rounded}) norm(w)_2") # beta_rounded corresponds to C_2 in Lemma 2.2 in paper
    print(f"is heuristically less than 2^-{secparam}.")
    print()
    print(f"This gives a slack of {ceil_decimal(sqrt(projdim/(2*beta)), 2)} for the approximate norm check.")

# Run function for security parameter lambda = 128 -> used for Falcon-512
print_jl_msg(128)
# Run function for security parameter lambda = 256 -> used for Falcon-1024
print_jl_msg(256)
