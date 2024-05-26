''' This python script computes the estimated sizes of our aggregate signature as presented in Section 6 and Appendix E of our submission.
It is derived from corresponding code provided to us by Gregor Seiler.
It first provides the numbers for the comparision with [JRS23], Squirrel and Chipmunk.
Then it provides the numbers stored in 'estimates-lin-.csv' to derive the plots through 'plot_paper.py' later.
To run the program, simply run the python script, for example >>> python3 ./proof_size_estimate.py '''

import math
import sys
import time
import csv
from enum import Enum

### FALCON PARAMETER SETS ###

''' Defining parameter sets for FALCON_DEG_SEC:
	Parameters depend on ring degree DEG and security level SEC, thus different classes for the relevant combinations

	d = ring degree
	q = Falcon modulus
	beta = norm bound of signatures
	bit_lin = bit length of a signature without nonce (from Falcon specifications, Table 3.3 and Algorithm 10)
	SECPARAM = aimed security level
	JL_const = the constant in the Johnson-Lindenstrauss Lemma (for the corresponding security level)
	JL_slack = the slack in the Johnson-Lidenstrauss Lemma (for the given security level)
	kappa_lim[kappa] = largest beta for which that kappa gives us SECPARAM-bit MSIS security (using SIS_hardness.py script)
'''

# all configurations we need        
class FALCON_64_128(): 
    def __init__(self) -> None:
        self.d = 64
        self.q = 12289
        self.beta = 5834
        self.bit_len = 5328-320 
        self.SECPARAM = 128
        self.JL_const = 120
        self.JL_slack = math.sqrt(128/30)
        self.kappa_lim = [-1, 271, 2687, 16383, 73727, 524287, 917503, 2621439, 7340031, 18874367, 50331647, 117440511, 251658239, 1073741823, 1207959551, 2684354559, 5368709119, 10737418239, 21474836479, 68719476735, 137438953471, 137438953471, 240518168575, 549755813887, 755914244095, 1374389534719, 4398046511103, 4398046511103, 8796093022207, 13194139533311, 35184372088831, 35184372088831, 52776558133247, 87960930222079, 140737488355326, 140737488355326]
      
class FALCON_128_128(): 
    def __init__(self) -> None:
        self.d = 128
        self.q = 12289
        self.beta = 5834
        self.bit_len = 5328-320 
        self.SECPARAM = 128
        self.JL_const = 120
        self.JL_slack = math.sqrt(128/30)
        self.kappa_lim = [-1, 2687, 73727, 917503, 7340031, 50331647, 251658239, 1207959551, 5368709119, 21474836479, 137438953471, 240518168575, 755914244095, 4398046511103, 8796093022207, 35184372088831, 52776558133247, 140737488355326, 140737488355326]

class FALCON_128_256(): 
    def __init__(self) -> None:
        self.d = 128
        self.q = 12289
        self.beta = 5834
        self.bit_len = 10240-320
        self.SECPARAM = 256
        self.JL_const = 168
        self.JL_slack = math.sqrt(128/30)
        self.kappa_lim = [-1, 463, 5887, 49151, 229375, 917503, 3407871, 11534335, 34603007, 100663295, 268435455, 805306367, 1744830463, 4294967295, 9663676415, 21474836479, 47244640255, 103079215103, 206158430207, 549755813887, 1099511627775, 1649267441663, 3161095929855, 6047313952767, 13194139533311, 21990232555519, 39582418599935, 70368744177663, 140737488355326, 140737488355326]
   
class FALCON_256_256(): 
    def __init__(self) -> None:
        self.d = 256
        self.q = 12289
        self.beta = 5834
        self.bit_len = 10240-320
        self.SECPARAM = 256
        self.JL_const = 168
        self.JL_slack = math.sqrt(128/30)
        self.kappa_lim = [-1, 5887, 229375, 3407871, 34603007, 268435455, 1744830463, 9663676415, 47244640255, 206158430207, 1099511627775, 3161095929855, 13194139533311, 39582418599935, 140737488355326, 140737488355326]


### LABRADOR CHALLENGE SETS ###

''' Defining challenge sets for FALCON_DEG_SEC in the case of
    1) two-splitting rings  (2_SPLIT)
    2) almost fully-splitting rings (ALMOST_FULL_SPLIT)
    
    (Parameters depend on ring degree DEG and security level SEC, thus different classes for the relevant combinations)

    rho 	= number of CRT slots (ell in the paper)
    tau 	= bound on the square of l_2 norm (given by the weight, T_2 in the paper)
    T 		= bound on operator norm (given by the weight, T_op in the paper)
    rep 	= boolian to indicate whether parallel repetition is necessary

    Parameters are set such that 

    i) size of challenge set > 2^SECPARAM
    ii) probability of non-invertibility < 2^(-SEPARAM)

    using the SageMath program 'compute_weight.sage' in this repository
'''

# best subring -> best parameter sets        
class CHAL_2_SPLIT_64_128():
    def __init__(self) -> None:
        self.rho = 2
        self.tau = math.ceil(172/2) #omega * gam^2 = 43*(2**2)=172
        self.T = math.ceil(86/2) #omega * gam = 43*2 =86
        self.rep = False
        
class CHAL_2_SPLIT_128_256():
    def __init__(self) -> None:
        self.rho = 2
        self.tau = math.ceil(296/2.5) #omega * gam^2 = 74*(2**2)=296
        self.T = math.ceil(148/2.5) #omega * gam = 74*2=148
        self.rep = False

# best subring -> best parameter sets                
class CHAL_ALMOST_FULL_SPLIT_128_128():
    def __init__(self) -> None:
        self.rho = 128//4
        self.tau = math.ceil(1024/2) # 64*(4**2)=1024
        self.T =  math.ceil(256/2) # 64*4=256
        self.rep = False

class CHAL_ALMOST_FULL_SPLIT_256_256():
    def __init__(self) -> None:
        self.rho = 256//8
        self.tau = math.ceil(1296/2.5) #144*(3**2)=1296
        self.T = math.ceil(432/2.5) #144*3=432
        self.rep = False   

#beta_fun is a lambda that takes kappa as input and outputs a norm. Needed for Thm 5.1.
def get_kappa(beta_fun, q_,kappa_lim):
    kappa = 0
    beta = 0
    for i in range(1, len(kappa_lim)):
        beta = beta_fun(i)
        if beta <= kappa_lim[i]:
            kappa = i 
            break

    if beta >= q_:
        raise Exception("Beta must be smaller than the LaBRADOR modulus.")
    if kappa == 0:
        raise Exception("The table kappa_lim has not found a kappa for such a large beta", beta)
    return kappa

##### Printing sizes #####

def format_size(nbits):
    # Should we print KiB instead?
    nbytes = math.ceil(nbits / 8)
    nkb = round(nbytes / 1000, 2)
    nmb = round(nbytes / 1000**2, 2)
    ngb = round(nbytes / 1000**3, 2)
    if nbytes < 1000:
        return "{v:10.2f} B".format(v=nbytes)
    if nkb < 1000:
        return "{v:10.2f} kB".format(v=nkb)
    if nmb < 1000:
        return "{v:10.2f} MB".format(v=nmb)
    return "{v:10.2f} GB".format(v=ngb) 

#### LaBRADOR script utils ####

def gaussianentropy(sig):
    a = 1
    if (sig >= 4):
        a = math.floor(sig/2)
        sig /= a;

    d = 1/(2*sig**2)
    n = 0
    for i in range(-math.ceil(15*sig), 0):
        n += math.exp(-i**2*d)
    n = 2*n + 1
    logn = math.log(n)
    e = 0
    for i in range(-math.ceil(15*sig), 0):
        f = math.exp(-i**2*d)
        e += f*(math.log(f) - logn)
    e = (-2*e + logn)/(n*math.log(2.0))

    return(e+math.log2(a))

def l2norm(list):
    return math.sqrt(sum(list[i]**2 for i in range(len(list))))

#####   Starting params  ######

def get_initial_params(num_sigs, falcon, chal, scal, verbose=False):
    ## signature l2-norm (contains both s_i's and epsilon_i's as well as their conjugates):
    beta_sig_ell2 = math.sqrt(num_sigs) * 2 * falcon.beta # called beta_2^(1) in paper
    ## quotient l2-norm (contains v_i's for lifting the modulus):
    beta_quotient_ell2 = math.sqrt(num_sigs) * (falcon.beta* falcon.d + math.sqrt(falcon.d) + 1) # called beta_2^(2) in paper
    ## sum of both l2-norm bounds gives final l2-norm bound
    beta_labrador = beta_sig_ell2 + beta_quotient_ell2
    
    # Condition 1
    # for beta_2^(1)
    bound_cond_1_sig = math.sqrt(falcon.SECPARAM) * beta_sig_ell2 * falcon.JL_const
    # for beta_2^(2)
    bound_cond_1_quo = math.sqrt(falcon.SECPARAM) * beta_quotient_ell2 * falcon.JL_const
        
    # Condition 2
    # for beta_inf^(1)
    bound_cond_2_sig = (falcon.JL_slack)**2 * (beta_sig_ell2)**2 * 4 * (falcon.d + 2)
    # for beta_inf^(2)
    bound_cond_2_quo = falcon.JL_slack * beta_quotient_ell2 * 6 * falcon.q
    
    
    q_bitlen = math.ceil( max(math.log2(bound_cond_1_sig), math.log2(bound_cond_1_quo), math.log2(bound_cond_2_sig), math.log2(bound_cond_2_quo)))
    
    '''    
    if verbose:
        print(num_sigs)
        print("beta_2^1: ", beta_sig_ell2, "beta_2^2: ", beta_quotient_ell2)
        print("bound through condition 1 for signature: ", bound_cond_1_sig)
        print("bound through condition 1 for quotient: ", bound_cond_1_quo)
        print("bound through condition 2 for signature: ", bound_cond_2_sig)
        print("bound through condition 2 for quotient: ", bound_cond_2_quo)
        print("bit length of Labrador modulus: ", q_bitlen)
    '''    
    if falcon.SECPARAM > q_bitlen * (falcon.d / chal.rho):
    	print("ERROR: q^{d/l} not big enough")
        
    q_ = 2**q_bitlen-1
    
    ## Starting Constraints
    n = scal * num_sigs
    r = 6 * math.ceil(math.sqrt(num_sigs)) + 1
    return q_, n, [r], [beta_labrador,0]

###############################
class Stage(Enum):
    FIRST = 1
    MID = 2 
    SECLAST = 3
    LAST = 4 #Meant for last round optimization, not used any more.

class Iteration():
    def __init__(self, q_, d, slack, n, r_list, beta_list, chal, secparam, kappa_lim, stage=Stage.MID, prevnu=1, prevmu=1) -> None:
        self.q_ = q_
        self.d = d
        self.slack = slack
        self.n = n 
        self.r_list = r_list
        self.beta_list = beta_list
        self.chal = chal
        self.stage = stage
        self.prevnu = prevnu
        self.prevmu = prevmu
        self.secparam = secparam
        self.kappa_lim = kappa_lim
    
        #Initializing internal variables
        if stage == Stage.LAST:
             self.init_last()
        else:
             self.init_normal()

    def init_normal(self):
        #Combined norm of the old witness vectors.
        self.beta = l2norm(self.beta_list)
        self.logq = math.ceil(math.log2(self.q_))

        #sig=sigma, the usual symbol for standard derivation.
        #They assume s_i are Gaussian with this standard deviation (p.16).
        self.sigs = [self.beta_list[i]/math.sqrt(self.r_list[i]*self.n*self.d) for i in range(0, len(self.r_list))] 

        #In the paper they model z as having SD sigs*sqrt(r*tau).
        #This seems to be the new SD of z_0 and z_1, when they keep track of the norm and SD of the 
        #subcomponents of the witness separately.
        self.sigz = math.sqrt(self.sigs[0]**2*(1.0+(self.r_list[0]-1)*self.chal.tau) + sum([self.sigs[i]**2*self.r_list[i]*self.chal.tau for i in range(1, len(self.r_list))]))

        #An upper bound on the standard deviation on the g_ij's (p.19, max not mentioned)
        self.sigh = math.sqrt(2*self.n*self.d)*max(self.sigs)**2

        if self.stage == Stage.SECLAST:
            self.t,self.b = 1, 1
        else:
            #z is decomposed into t=2 parts wrt to the basis b (p.16)
            self.t,self.b = 2,round(math.sqrt(math.sqrt(12.0)*self.sigz))

        #t_i and h_ij are split into t1 parts wrt the basis b1.
        #t1 is different compared to the paper. I assume the max is to make sure 
        #we don't get weird dividing by log(1) issues.
        self.t1 = round(self.logq/math.log2(math.sqrt(12.0)*self.sigz/self.b))
        self.t1 = max(2,self.t1)
        self.t1 = min(14,self.t1)
        self.b1 = math.ceil(2**(self.logq/self.t1))

        #g_ij is split into t2 parts wrt the basis b2.
        #t2 and b2 are both different from the paper.
        self.t2 = round(math.log(math.sqrt(12)*self.sigh)/math.log(math.sqrt(12)*self.sigz/self.b))
        self.t2 = max(1,self.t2)
        self.b2 = math.ceil((math.sqrt(12)*self.sigh)**(1/self.t2))

        #Total number of witness elements
        sumr = sum(self.r_list)
        #Computes the combined l2-norm upper bound for z_0 and z_1, sqrt(2/b**2 * gamma) (see p.15)
        #gamma = sigz *sqrt(nd), because the norm contribution of each Zq coefficient is on average the SD,
        #and there are nd coeffcients.
        self.nextbeta_list = [self.sigz/self.b*math.sqrt(float(self.t)*self.n*self.d), 0]
        
        #Rank of inner commitments
        newbeta1_fun = lambda kappa: math.sqrt(self.b1**2/12*self.t1*sumr*kappa*self.d + (self.b1**2*self.t1+self.b2**2*self.t2)/12*(sumr**2+sumr)/2*self.d)
        newbeta_fun = lambda kappa: math.sqrt(self.nextbeta_list[0]**2 + newbeta1_fun(kappa)**2)
        kappa_norm = lambda kappa: max(6*self.chal.T*self.b*self.slack*newbeta_fun(kappa),2*self.b*self.slack*newbeta_fun(kappa)+4*self.chal.T*self.slack*self.beta)
        self.kappa = get_kappa(kappa_norm, self.q_,self.kappa_lim)
        self.nextbeta_list[1] = newbeta1_fun(self.kappa)

        #The rank of the outer commitments. 
        kappa1_norm = lambda kappa: 2*self.slack*newbeta_fun(kappa)
        self.kappa1 = get_kappa(kappa1_norm, self.q_,self.kappa_lim)

        self.m = self.t1*sumr*self.kappa + (self.t1+self.t2)*(sumr**2+sumr)/2

    #Tail is the version of the last iteration that is in the paper.
    def init_last(self):
        #Same as main
        self.beta = l2norm(self.beta_list)
        self.logq = math.ceil(math.log2(self.q_))
        # size += 4*128;  # challenges

        #Same as main
        self.sigs = [self.beta_list[i]/math.sqrt(self.r_list[i]*self.n*self.d) for i in range(0, len(self.r_list))]
        self.sigz = math.sqrt(self.sigs[0]**2*(1.0+(self.r_list[0]-1)*self.chal.tau) + sum([self.sigs[i]**2*self.r_list[i]*self.chal.tau for i in range(1, len(self.r_list))]))
        self.sigh = math.sqrt(2*self.n*self.d)*max(self.sigs)**2

        #The norm of z when it is not decomposed into z_0 and z_1.
        self.beta = self.sigz*math.sqrt(self.n*self.d)
        #Picks the inner com rank according to Theorem 5.1, but with the new expression for the norm of z
        #We also don't need to multiply with the slack since we are not going to recurse further.
        kappa_norm = lambda kappa: max(6*self.chal.T*self.beta,2*self.beta+4*self.chal.T*self.beta)
        self.kappa = get_kappa(kappa_norm, self.q_,self.kappa_lim)
        
        #For printing purposes
        self.kappa1 = 0
        self.m = 0

    def size_step(self):
        numproj = 2 if self.stage == Stage.FIRST else 1
        jl_proj = numproj * 2 * self.secparam*gaussianentropy(self.beta/math.sqrt(2.0))

        jl_proof = math.ceil(self.secparam/self.logq)*self.d*self.logq

        numouter = 0 if self.stage == Stage.SECLAST else 2
        outer_com = numouter * self.kappa1 * self.d * self.logq

        return jl_proj + jl_proof + outer_com
    

    def size_ti(self):
        return sum(self.r_list)*self.kappa*self.d*self.logq
    
    def size_gij(self):
        if self.stage == Stage.LAST:
            #According to Section 5.6, there will be 2*nu+1 g_ij polynomials
            return (1+2*self.r_list[0])*self.d*gaussianentropy(self.sigh) # quadratic garbage polys, assumes t = 1
        
        return ((self.r_list[0]**2 + self.r_list[0])/2)*self.d*gaussianentropy(self.sigh)

    def size_hij(self):
        if self.stage == Stage.LAST:
            #According to Section 5.6, there will be 2r-1 h_ij polynomials.
            return (2*sum(self.r_list)-1)*self.d*self.logq  # linear garbage polys
        
        return ((self.r_list[0]**2 + self.r_list[0])/2)*self.d*self.logq
    
    def size_lastmsg(self):
        return self.size_ti() + self.size_gij() + self.size_hij() + self.size_z()
    
    def size_all(self):
        return self.size_step() + self.size_lastmsg()

    def parallel_reps(self):
        if not self.chal.rep:
            return 1
        for k in range(1,20):
            if self.secparam < math.log2( k / self.d) + k * math.log2(self.q_): # <=> 2^(-secparam) > d/k * (1/q)^k
                return k     
        # TODO this is super hacky, but we should never end in this case
        return 21
        
    def size_z(self):
        if self.stage == Stage.LAST:
            #TODO: times two to avoid exponential loss in weak opening norms (ask sebastian)
            return 2* self.parallel_reps()*self.n*self.d*gaussianentropy(self.sigz)
        return self.parallel_reps()*self.n*self.d*gaussianentropy(self.sigz)

    def next_it(self, recursion_strategy, nextstage=Stage.MID):
        nu, mu = recursion_strategy(self, nextstage)
        n_ = math.ceil(self.parallel_reps()*self.n/nu)
        m_ = math.ceil(self.m/mu) #Without the ceils and floors, this would make newm=n/nu=newn.
        n_ = max(n_,m_) #New rank
        r_ = [self.t*nu,mu]; #New multiplicity.
        it = Iteration(self.q_, self.d, self.slack, n_, r_, self.nextbeta_list, self.chal, self.secparam, self.kappa_lim, nextstage, nu, mu,)
        return it

    def __str__(self):
        return f"n: {self.n:10d}, r: {sum(self.r_list):10d}, m: {math.ceil(self.m):10d}, ν: {self.prevnu:4d}, μ: {self.prevmu:4d}, kappa: {self.kappa:4d}, kappa1: {self.kappa1:4d}, log(β) {math.log2(self.beta):5.2f} norm balance: {100*(self.beta_list[1]-self.beta_list[0])/max(self.beta_list):.2f} size: {format_size(self.size_step())}{', ' +format_size(self.size_lastmsg())}"

    def lastmsg_str(self):
        #  self.size_z() + self.size_t_i() + self.size_g_ij() + self.size_hij()
        return f"z: {format_size(self.size_z())}, ti: {format_size(self.size_ti())}, gij: {format_size(self.size_gij())}, hij: {format_size(self.size_hij())}"

#### Recursion ####

def recursion_strategy_4(it, nextstage):
    best_nu, best_mu = 1,1
    best = it.next_it((lambda a,b: (1,1)), nextstage)
    best_size = best.parallel_reps() * 2 * best.n +  best.m
    res = 50
    for nu in range(1,res):
        for mu in range(1,res):
            candidate = it.next_it((lambda a,b: (nu,mu)), nextstage)
            cand_size =  candidate.parallel_reps() * 2 * candidate.n +  candidate.m
            if  cand_size < best_size:
                best_size = cand_size
                best_nu, best_mu = nu, mu
    # print(best_nu,best_mu)
    return best_nu,best_mu

def recursion_strategy_final_4(it, nextstage):
    best_nu, best_mu = 1,1
    best = it.next_it((lambda a, b: (1,1)), Stage.LAST)
    best_size = best.size_all()
    res = 100
    for nu in range(1,res):
        for mu in range(1,res):
            candidate = it.next_it((lambda a,b: (nu,mu)), Stage.LAST)
            if candidate.size_all() < best_size:
                best_size = candidate.size_all()
                best_nu, best_mu = nu, mu
    # print(best_nu,best_mu)
    return best_nu,best_mu

def recursion_to_depth(initial_it, depth):
    iters = [initial_it]
    it = initial_it
    strategies = [] if depth <= 3 else [recursion_strategy_4] * (depth-3) 
    for strat in strategies:
        it = it.next_it(strat, Stage.MID)
        iters.append(it)

    seclast = it.next_it(recursion_strategy_4, Stage.SECLAST)
    last = seclast.next_it(recursion_strategy_final_4, Stage.LAST)
    iters.append(seclast)
    iters.append(last)

    total = 0
    for it in iters:
        total += it.size_step()
    
    total += iters[len(iters)-1].size_lastmsg()

    return iters, total

def recursion_to_depth_no_last_opt(initial_it, depth):
    iters = [initial_it]
    it = initial_it
    strategies = [] if depth <= 2 else [recursion_strategy_4] * (depth-2) 
    for strat in strategies:
        it = it.next_it(strat, Stage.MID)
        iters.append(it)

    seclast = it.next_it(recursion_strategy_4, Stage.SECLAST)
    iters.append(seclast)

    total = 0
    for it in iters:
        total += it.size_step()
    
    total += iters[len(iters)-1].size_lastmsg()

    return iters, total

def search(num_sigs, max_depth, falcon, chal, scal, verbose=False):
    q_, n, r_list, beta_labrador_list = get_initial_params(num_sigs, falcon, chal, scal, verbose=verbose)
    if verbose: print("Initial Modulus (in bits):", math.ceil(math.log2(q_)))

    it = Iteration(q_, falcon.d, falcon.JL_slack, n, r_list, beta_labrador_list, chal, falcon.SECPARAM, falcon.kappa_lim, Stage.FIRST)
    if verbose: print(f'Requires {it.parallel_reps()} parallel repetition(s)')

    best_iterations, best_size = recursion_to_depth_no_last_opt(it, 1)

    for d in range(2, max_depth + 1):
        iterations, size = recursion_to_depth_no_last_opt(it, d)
        if size < best_size:
            best_iterations, best_size = iterations, size
    if verbose:
        #for depth, it in enumerate(best_iterations):
        #    print(f"{depth:2d}", it)
        #print("Last message:",best_iterations[len(best_iterations)-1].lastmsg_str())
        print("Total (with salt):", format_size(best_size + 320*num_sigs))
        print("Total (without salt) ", format_size(best_size))
        print("Trivial solution (with salt): ", format_size(num_sigs*falcon.bit_len + 320*num_sigs))
        print("Trivial solution (without salt): ", format_size(num_sigs*falcon.bit_len))
        print("Compression rate (with salt): ", round((best_size + 320*num_sigs)/(num_sigs*falcon.bit_len + 320*num_sigs),2))

    return best_size

## FUNCTIONS TO COMPUTE NUMBERS AND STORE IN CSV FILES
    
def compute_sizes_best(file_name, sizes):
    f = open(file_name,"a")
    max_depth = 15
  
    f.write("Num-Sigs,Naive-512,Naive-1024,Falcon-512-2S,Falcon-512-AS,Falcon-1024-2S,Falcon-1024-AS\n")
    for num_sigs in sizes:
        sys.stdout.write('\r')
        sys.stdout.write(f'{num_sigs}')
        sys.stdout.flush()
        # Falcon-512 2S (two-splitting), AS (almost-fully-splitting)
        f512_2S = search(num_sigs, max_depth, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8)
        f512_AS = search(num_sigs, max_depth, FALCON_128_128(), CHAL_ALMOST_FULL_SPLIT_128_128(),4)
	# Falcon-1024 2S (two-splitting), AS (almost-fully-splitting)
        f1024_2S = search(num_sigs, max_depth, FALCON_128_256(), CHAL_2_SPLIT_128_256(),8)
        f1024_AS = search(num_sigs, max_depth, FALCON_256_256(), CHAL_ALMOST_FULL_SPLIT_256_256(),4)
        # Also include naive concatenative 
        f.write(f'{num_sigs}, {num_sigs * FALCON_128_128().bit_len}, {num_sigs * FALCON_256_256().bit_len}, {f512_2S}, {f512_AS}, {f1024_2S},{f1024_AS} \n')
        f.flush()
    f.close()
    print("Done!")

## NEW WITH BEST SUBRING 
# COMPARISON WITH [JRS23]

print("\n RESULTS FOR COMPARISON WITH [JRS23] (Table 4 in the submission)")
print("\n For Falcon-512 (degree 64, secpar 121) and N=500")
search(500, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)
print("\n For Falcon-512 (degree 64, secpar 121) and N=1000")
search(1000, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)
print("\n For Falcon-512 (degree 64, secpar 121) and N=2000")
search(2000, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)

# COMPARISON WITH SQUIRREL AND CHIPMUNK

print("\n RESULTS FOR COMPARISON WITH SQUIRREL & CHIPMUNK (Table 5 in the submission)")
print("\n For Falcon-512 (degree 64, secpar 121) and N=1024")
search(1024, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)
print("\n For Falcon-512 (degree 64, secpar 121) and N=4096")
search(4096, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)
print("\n For Falcon-512 (degree 64, secpar 121) and N=8192")
search(8192, 15, FALCON_64_128(), CHAL_2_SPLIT_64_128(),8, True)

# Computes all relevant aggregate signature sizes to run plot.py to obtain the different figures in Section 6.2
compute_sizes_best("estimates-lin-.csv",  range(100, 10101, 100)) # start, end, step size
