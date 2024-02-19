''' This SageMath script computes the estimated hardness of Module-SIS needed for our aggregate signature. 
It is derived from corresponding code provided to us by Gregor Seiler.
Its output is then used as 'kappa_lim' in the python script 'proof_size_estimate.py' 
Concretely, it computes the corresponding 'kappa_lim' lists for the classes
	1) FALCON_64_128 (used for two-splitting 128-bit security)
	2) FALCON_128_128 (used for almost-fully-splititing 128-bit security
	3) FALCON_128_256 (used for two-splitting 256-bit security)
	4) FALCON_256_256 (used for almost-fully-splitting 256-bit security)
To run the program, simply run the SageMath script, for example >>> load("SIS_hardness.sage")'''

def GSAhermite(bs):
    return (bs/(2*pi*e)*numerical_approx(pi*bs)^(1/bs))^(1/(2*bs-2))

def findblocksize(h):
    bs=9999
    for i in range(50, 1000):
        if GSAhermite(i) <= h:
            bs = i
            break
    return bs

def BDGLcost(bs):
    return bs*log(sqrt(3/2),2)

# m: rows
# q: modulus
# b: L2 norm bound
def SIShardness(m, q, b, verbose=False):
    if b >= q:
        print("Norm bound bigger than q!\n")
        return
    
    # When b = q^(m/d) h^d , h is maximal at d = 2 m log_b(q)
    d = round(2*m*log(q, b))
    h = numerical_approx((b*q^(-m/d)))^(1/d) # calc root HF for lattice dimension d
    #print(d)
    #print(h)
    #for i in range(m+1, n, 5):
    #    t = (b*q^(-m/i))^(1/i) # calc root HF for col i (lattice dim)
    #    if t > h: 
    #        h = t
    #        d = i
    
    bs = findblocksize(h)
    if bs > d:
        #print("Required block size bigger than lattice dimension!\n")
        return 292
    
    cost = BDGLcost(bs)
    if verbose:
        print("SIS hardness in l2 norm:\n");
        print("  SIS lattice dimension used: %d\n" % d);
        print("  Root Hermite factor: %.5f\n"% h);
        print("  BKZ block size: %d\n" % bs);
        print("  Classical Core-SVP bit-cost for BDGL16 sieve: %.2f\n" % cost)

    return cost

#Binary search for the largest beta such that SIS is hard
def binary_search(sec, d, q, kappa, max_beta):
    l = 1
    r = max_beta
    cost = 0
    while l < r:
        mid = (l+r)//2
        #print("mid", mid)
        cost = SIShardness(m=kappa*d, q=q, b=mid)
        if abs(sec-cost) < 1 and cost > sec: # stop if the cost is close enough the the desired security level
            break
        if cost < sec: 
            r = mid
        else: 
            l = mid + 1
        
    print("Classical Core-SVP bit-cost for BDGL16 sieve: %.2f\n" % cost)
    return r-1

def build_table(sec, d, q, max_kappa):
    l = [-1]
    for kappa in range(1, max_kappa+1): # change
        beta = binary_search(sec, d, q, kappa, q)
        l.append(beta)
        print("kappa", kappa, "beta", beta, "\n")
    return l

###### ====== SEC PARAMETER 128 ===== ######

print('Security Parameter 128')
sec = 128

# class FALCON_64_128
d = 64 # ring degree 
logq = 47 #Lower bound for q, for hardness estimation
q = 2**logq-1
max_kappa = 35 #max kappa chosen such that maximal beta value is reached

l = build_table(sec, d, q, max_kappa)
print('Ring Degree 64 - Security Parameter 128')
print(l)

# class FALCON_128_128
d = 128 # ring degree 
logq = 47 #Lower bound for q, for hardness estimation
q = 2**logq-1
max_kappa = 18 #max kappa chosen such that maximal beta value is reached

l = build_table(sec, d, q, max_kappa)
print('Ring Degree 128 - Security Parameter 128')
print(l)

###### ====== SEC PARAMETER 256 ===== ######

print('Security Parameter 256')
sec = 256

# class FALCON_128_256
d = 128 # ring degree 
logq = 47 #Lower bound for q, for hardness estimation
q = 2**logq-1
max_kappa = 29 #max kappa chosen such that maximal beta value is reached

l = build_table(sec, d, q, max_kappa)
print('Ring Degree 128 - Security Parameter 256')
print(l)

# class FALCON_256_256
d = 256 # ring degree 
logq = 47 #Lower bound for q, for hardness estimation
q = 2**logq-1
max_kappa = 15 #max kappa chosen such that maximal beta value is reached

l = build_table(sec, d, q, max_kappa)
print('Ring Degree 256 - Security Parameter 256')
print(l)

# For d = 64 
#   128-bit security l = [-1, 271, 2687, 16383, 73727, 524287, 917503, 2621439, 7340031, 18874367, 50331647, 117440511, 251658239, 1073741823, 1207959551, 2684354559, 5368709119, 10737418239, 21474836479, 68719476735, 137438953471, 137438953471, 240518168575, 549755813887, 755914244095, 1374389534719, 4398046511103, 4398046511103, 8796093022207, 13194139533311, 35184372088831, 35184372088831, 52776558133247, 87960930222079, 140737488355326, 140737488355326]
# after that just repetitions
#
# For d = 128:
#   128-bit security l = [-1, 2687, 73727, 917503, 7340031, 50331647, 251658239, 1207959551, 5368709119, 21474836479, 137438953471, 240518168575, 755914244095, 4398046511103, 8796093022207, 35184372088831, 52776558133247, 140737488355326, 140737488355326]
#   256-bit security l = [-1, 463, 5887, 49151, 229375, 917503, 3407871, 11534335, 34603007, 100663295, 268435455, 805306367, 1744830463, 4294967295, 9663676415, 21474836479, 47244640255, 103079215103, 206158430207, 549755813887, 1099511627775, 1649267441663, 3161095929855, 6047313952767, 13194139533311, 21990232555519, 39582418599935, 70368744177663, 140737488355326, 140737488355326]
# after that just repetitions
#
# For d = 256
#   256-bit security l = [-1, 5887, 229375, 3407871, 34603007, 268435455, 1744830463, 9663676415, 47244640255, 206158430207, 1099511627775, 3161095929855, 13194139533311, 39582418599935, 140737488355326, 140737488355326]
# after that just repetitions
