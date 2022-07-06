import pickle
from progress.bar import Bar
import logging
import math
import numpy as np
from Crypto.Util import number
import galois
from lattice import *
from codetiming import Timer
from datetime import datetime


def prime(i, primes):
    for prime in primes:
        if not (i == prime or i % prime):
            return False
    primes.add(i)
    return i


def n_primes(n):
    primes = set([2])
    i, p = 2, 0
    while True:
        if prime(i, primes):
            p += 1
            if p == n:
                return list(primes)
        i += 1


class LinearHomoF2:
    def __init__(self, equations):
        self.Fp = galois.GF(2)
        self.dim = len(equations)
        self.A = self.Fp(equations)
        self.t = self.A.left_null_space()
        self.max = 2**len(self.t)
        self.dp = -1
        self.has_been_zero = 0
        if(len(self.t > 0)):
            self.dp = 0
        else:
            logging.warning(
                "No non-trivial solutions exist to the system of Homogenoeus equations")
            self.has_been_zero = 1
        logging.info("Solution space is {}".format(self.t))

    def bin_array(self, num, m):
        """Convert a positive integer num into an m-bit bit vector"""
        return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8).tolist()

    def next(self):
        sol = self.Fp([0]*self.dim)
        if self.has_been_zero == 1:
            logging.warning("No more solutions exist, returning empty array")
            return sol

        bits = self.bin_array(self.dp, len(self.t))
        for i in range(len(bits)):
            sol += bits[i]*self.t[i]
        self.dp = (self.dp + 1) % self.max
        if(self.dp == 0):
            self.has_been_zero = 1
        return sol


def primesfrom2to(n):
    sieve = np.ones(n//3 + (n % 6 == 2), dtype=bool)
    sieve[0] = False
    for i in range(int(n**0.5)//3+1):
        if sieve[i]:
            k = 3*i+1 | 1
            sieve[((k*k)//3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i & 1))//3::2*k] = False
    return np.r_[2, 3, ((3*np.nonzero(sieve)[0]+1) | 1)].tolist()


def prime_base(N, alpha):
    '''
    list of all primes from 2 to floor(N^alpha)
    '''
    n = int(math.log(N)**alpha)
    return primesfrom2to(n)


def is_smooth(x, P):
    '''
    Checks if a given number is pt-smooth
    '''
    y = x
    for p in P:
        while p.divides(y):
            y /= p
    return abs(y) == 1


def factorize_smooth(n, primes):
    '''
    returns empty array if not pt-smooth, otherwise returns factorization
    same behavior for negative numbers
    '''
    if(n < 0):
        n = -n
    exponents = [0]*len(primes)
    for i in range(len(primes)):
        if(n == 1):
            return exponents
        if(n < primes[i]):
            return []
        while(n % primes[i] == 0):
            n = n//primes[i]
            exponents[i] += 1
    if(n == 1):
        return exponents
    return []


def sr(x, prec):
    '''
    rounding funtion
    '''
    return round(x * (10 ** prec))


def generate_basis(prime_base: list, multiplier: int, prec: int) -> IntegerMatrix:
    '''
    generates basis matrix as per schnorr's original algorithm
    [ln(2) x (10^prec) ...........  multiplier x ln(2) x (10^prec)]
    [0 ln(3) x (10^prec) ....                                     ]
    .
    .
    .
    [..............................multiplier x ln(pt) x (10^prec)]

    '''

    n = len(prime_base)

    B = [[0]*(n+1) for _ in range(n)]
    for i in range(n):
        B[i][n] = sr(multiplier*math.log(prime_base[i]), prec)
        B[i][i] = sr(math.log(prime_base[i]), prec)

    # B = IntegerMatrix.from_matrix(B)
    return B


def check_fac_relation(e: list, prime_base: list, N: int, s_array = None) -> list:
    '''
    given e vector, checks if s = u - v*N is pt-smooth, returns their factorization if true
    '''
    assert len(e) == len(prime_base)
    [u, v] = [1, 1]
    T = [0]*len(prime_base)
    for i in range(len(e)):
        if(e[i] < 0):
            exp = -1*e[i]
            v *= prime_base[i]**(exp)
        else:
            u *= prime_base[i]**e[i]
            T[i] = e[i]

    s = u - v*N
    if(s_array != None):
        s_array.append(s)
    S = factorize_smooth(s, prime_base)
    S.insert(0, int(s < 0))
    T.insert(0, 0)

    if(len(S) == 1):
        return []
    return [tuple([u, s]), tuple(T), tuple(S)]


def fac_relations_cvp(N, P, c, prec=10, independent=False):
    # 0. not short enough
    # 1. working
    # -1. repeated

    timer = Timer(logger=None)
    ret = []  # return time

    n = len(P)
    logging.info("dimension: {}".format(n))
    s_sizes = []
    # for i in P:
    #     if(N % i == 0):
    #         return [i, N/i]

    if(independent):
        bit_length = N.bit_length()
        multiplier = 2**(bit_length*c)
    else:
        multiplier = N**c

    # Basis matrix for the prime number lattice
    ret.append(n+1)
    timer.start()
    Basis = generate_basis(P, multiplier, prec)
    run_time = timer.stop()
    logging.warning("Basis matrix generated in {} seconds".format(run_time))
    ret.append(run_time)

    refs = [Basis[i][i] for i in range(len(P))]
    # Target vector for CVP [0,..., N^c*log(N)]
    target = [0]*(len(P)+1)
    target[-1] = sr(multiplier*math.log(N), prec)

    relations = {}
    # Solve CVP here
    trial = 0
    bar = Bar('Finding relations', max=n+2)

    timer.start()
    reductions = []
    succ = []
    while(len(relations) < n+2):
        # assuming the probablity of getting same permuation more then once is very low
        trial += 1
        np.random.shuffle(Basis)

        # B_reduced = bkz_reduce(Basis, 30)  # try tuning the block size
        temp_timer = Timer(logger=None)
        temp_timer.start()

        B_reduced = lll_reduce(Basis)
        e_reduced = cvp_babai(B_reduced, target)
        w = B_reduced.multiply_left(e_reduced)

        reduction_time = temp_timer.stop()
        reductions.append(reduction_time)
        logging.warning("* {} reduction runtime: {} seconds".format(
            trial, reduction_time))

        e = []
        for i in range(len(w)-1):
            assert w[i] % refs[i] == 0
            e.append(w[i]//refs[i])

        # implement a checker function without actually computing u and v to reject longer vectors, similar to SVP one in Ritter's paper.
        rel = check_fac_relation(e, P, N, s_sizes)

        if(len(rel) == 0):
            # amend it, do not continue reduce it strongly
            logging.info(">>>> not short enough")
            succ.append(1)
            continue

        key = rel[0]

        assert (key[0] - key[1]) % N == 0

        if key not in relations:
            logging.info(
                ">>>> found new fac-relation: {}".format(key))
            relations[key] = rel[1:]
            bar.next()
            succ.append(0)

        else:
            logging.info(">>>> relation already exists")
            succ.append(-1)

    run_time = timer.stop()
    logging.info("Found n+2 relations in {} seconds".format(run_time))
    ret.append(run_time)
    ret.append([reductions, succ])

    bar.finish()
    P.insert(0, -1)

    return [relations, ret]


def shuffle_basis(Basis, rank, dimension):
    Shuff = []

    # Copying the Basis matrix
    for x in range(rank):
        Shuff.append([])
        for y in range(dimension):
            Shuff[x].append(Basis[x][y])

    np.random.shuffle(Shuff)
    return IntegerMatrix.from_matrix(Shuff)


def ratios(succ):
    minus_one, zero, one = [0,0,0]
    for t in succ:
        if(t == -1):
            minus_one += 1
        if(t == 0):
            zero += 1
        if(t == 1):
            one += 1
    return [minus_one/len(succ),zero/len(succ),one/len(succ)]


def fac_relations_svp(N, P, c, prec=10, beta=30):
    '''
    If independent, similar to ritters version:
    c = c1 - c2
    prec = c2
    '''
    timer = Timer(logger=None)
    ret = []

    n = len(P)
    ret.append(n+2)
    multiplier = 10**c


    timer.start()
    # Ritters basis, size: n+1 x n+2
    Basis = generate_basis(P, multiplier, prec)
    refs = [Basis[i][i] for i in range(len(P))]

    for vecs in Basis:
        vecs.insert(0, 0)

    target = [0]*(len(P)+2)
    target[0] = 1
    target[-1] = sr(multiplier*math.log(N), prec)
    Basis.insert(0, target)

    run_time = timer.stop()
    logging.warning("Basis matrix generated in {} seconds".format(run_time))
    ret.append(run_time)

    rank = len(Basis)
    dimension = len(Basis[0])

    Basis = lll_reduce(Basis, 0.99)

    relations = {}
    bar = Bar('Finding relations', max=n+2)
    succ = []
    reductions = []

    Rounds = 0
    timer.start()

    while(len(relations) < n+2):
        temp_timer = Timer(logger=None)
        temp_timer.start()

        Rounds += 1
        Basis = shuffle_basis(Basis, rank, dimension)  # shuffl
        enums = BKZs(Basis)(beta)

        reduction_time = temp_timer.stop()
        reductions.append(reduction_time)
        logging.warning("* {} reduction runtime: {} seconds".format(
            Rounds, reduction_time))


        for vec in enums:
            if(len(relations) == n+2):
                break
            if(abs(vec[0]) != 1):
                continue
            
            e = []
            for i in range(1, dimension-1):
                assert vec[i]%refs[i-1] == 0
                e.append(vec[i]//refs[i-1])
            
            rel = check_fac_relation(e,P,N)

            if(len(rel) == 0):
                logging.info(">>>> not short enough")
                succ.append(1)
                continue
            key = rel[0]

            assert (key[0] - key[1]) % N == 0

            if key not in relations:
                logging.info(">>>> found new fac-relation: {}".format(key))
                relations[key] = rel[1:]
                bar.next()
                succ.append(0)

            else:
                logging.info(">>>> relation already exists")
                succ.append(-1)
        for vec in Basis:
            if(len(relations) == n+2):
                break
            if(abs(vec[0]) != 1):
                continue
            
            e = []
            for i in range(1, dimension-1):
                assert vec[i]%refs[i-1] == 0
                e.append(vec[i]//refs[i-1])
            
            rel = check_fac_relation(e,P,N)

            if(len(rel) == 0):
                logging.info(">>>> not short enough")
                succ.append(1)
                continue
            key = rel[0]

            assert (key[0] - key[1]) % N == 0

            if key not in relations:
                logging.info(">>>> found new fac-relation: {}".format(key))
                relations[key] = rel[1:]
                bar.next()
                succ.append(0)

            else:
                succ.append(-1)
                logging.info(">>>> relation already exists")

    run_time = timer.stop()
    logging.info("Found n+2 relations in {} seconds".format(run_time))
    ret.append(run_time)
    ret.append([reductions, succ])


    bar.finish()
    P.insert(0, -1)

    return [relations,ret]


def solve_linear(N, a_b, primes):
    a_plus_b_mod2 = []
    for i in range(len(a_b)):
        temp = []
        for j in range(len(a_b[i][0])):
            temp.append((a_b[i][0][j] + a_b[i][1][j]) % 2)
        a_plus_b_mod2.append(temp)
    solver = LinearHomoF2(a_plus_b_mod2)
    while(True):
        if(solver.has_been_zero == 1):
            break
        sol = solver.next()
        logging.info("* C = {}".format(sol))
        A = np.array([0]*len(a_b[0][0]))
        B = np.array([0]*len(a_b[0][1]))
        for i in range(len(sol)):
            if(sol[i] == 1):
                A += np.array(a_b[i][0])
                B += np.array(a_b[i][1])
        A += B
        A = A//2
        a = 1
        b = 1
        for i in range(len(A)):
            a *= primes[i]**int(A[i])  # Fuckin hate numpy
            b *= primes[i]**int(B[i])
        assert (a**2 - b**2) % N == 0

        x = math.gcd(N, a+b)
        y = math.gcd(N, a-b)

        logging.info(">> gcd(N, a+b): {}, gcd(N, a-b): {}".format(x, y))

        if(1 < x < N):
            print("Factor found: {}".format(x))
            return x

        if(1 < y < N):
            print("Factor found: {}".format(y))
            return y

    print("Better luck next time:(")
    return 1


def schnorr(N, alpha, c, prec):
    ret = []

    timer = Timer(logger=None)

    timer.start()
    P = prime_base(N, alpha)  # 1
    run_time = timer.stop()
    logging.warning(
        "Time for generating prime basis: {} seconds".format(run_time))
    logging.info("prime base: {}".format(P))
    ret.append(run_time)

    [relations, timing] = fac_relations_cvp(N, P, c, prec, False)  # 2
    ret.append(timing)

    a_b = list(relations.values())  # 3

    timer.start()
    fac = solve_linear(N, a_b, P)  # 4
    run_time = timer.stop()
    ret.append(run_time)
    logging.warning(
        "Time find a non-trivial solution to linear equations: {} seconds".format(run_time))
    if(fac == 1):
        ret.append(0)
    else:
        ret.append(1)
    return ret


def ritter(N, t, c1, c2, beta):
    ret = []

    timer = Timer(logger=None)
    timer.start()

    P = n_primes(t) # Generate list of first t primes

    run_time = timer.stop()
    logging.warning(
        "Time for generating prime basis: {} seconds".format(run_time))
    logging.info("prime base: {}".format(P))
    ret.append(run_time)

    [relations, timing] = fac_relations_svp(N, P, c=c1 - c2, prec=c2, beta=beta) # Get >= n+2 fac-relations
    a_b = list(relations.values())
    ret.append(timing)

    timer.start()
    fac = solve_linear(N, a_b, P) # Use linear solver to get a non trivial factor
    run_time = timer.stop()
    ret.append(run_time)
    logging.warning(
        "Time find a non-trivial solution to linear equations: {} seconds".format(run_time))
    
    if(fac == 1):
        ret.append(0)
    else:
        ret.append(1)
    return ret


def main_schnorr():
    bits_low = 20
    bits_high = 24
    bits_step = 2

    c_low = 1.1
    c_high = 1.3
    c_step = 0.05

    alpha_low = 1.6
    alpha_high = 1.9
    alpha_step = 0.1

    timing = {}
    timer = Timer(logger=None)

    now = datetime.now()
    current_date = str(now.date())
    current_time = str(now.time())

    logging.basicConfig(level=logging.INFO,
                        filename="./logs/" + current_date + current_time + ".log", format='%(message)s')
    for bits in range(bits_low, bits_high+1, bits_step):
        p = number.getPrime(bits//2)
        q = number.getPrime(bits//2)
        N = p*q
        print("N = {}:".format(N))
        logging.warning(
            "------------------------------------------------------------------------\nN = {}".format(N))
        alpha = alpha_low
        alphas = {}
        while(alpha < alpha_high):
            c = c_low
            cs = {}
            while(c < c_high):
                print("alpha = {}, c = {}".format(alpha, c))
                logging.warning("alpha = {}, c = {}".format(alpha, c))
                timer.start()
                ret = schnorr(N, alpha, c, 5)
                overall = timer.stop()
                ret.append(overall)
                logging.warning(
                    "Overall runtime: {} seconds".format(overall))
                cs[c] = ret
                c += c_step
            alphas[alpha] = cs
            alpha += alpha_step
        timing[bits] = alphas
    with open("./timing/" + current_date+current_time + '.pkl', 'wb') as fp:
        pickle.dump(timing, fp)


def main_ritter():
    bits_low = 20
    bits_high = 20
    bits_step = 2

    c1_low = 10
    c1_high = 12
    c1_step = 2

    c2_low = 3
    c2_high = 3
    c2_step = 1

    t_low = 125
    t_high = 125
    t_step = 10

    beta_low = 30
    beta_high = 30
    beta_step = 10

    timing = {}
    timer = Timer(logger=None)

    now = datetime.now()
    current_date = str(now.date())
    current_time = str(now.time())

    logging.basicConfig(level=logging.INFO,
                        filename="./logs/ritter/" + current_date + current_time + ".log", format='%(message)s')
    for bits in range(bits_low, bits_high+1, bits_step):
        p = number.getPrime(bits//2)
        q = number.getPrime(bits//2)
        N = p*q
        print("N = {}:".format(N))
        logging.warning(
            "------------------------------------------------------------------------\nN = {}".format(N))
        c1 = c1_low
        c1s = {}
        while(c1 <= c1_high):
            c2 = c2_low
            c2s = {}
            while(c2 <= c2_high):
                t = t_low
                ts = {}
                while(t <= t_high):
                    beta = beta_low
                    betas = {}
                    while(beta <= beta_high):
                        print("bits = {}, c1 = {}, c2 = {}, t = {}, beta = {}".format(bits, c1,c2,t,beta))
                        logging.warning("bits = {}, c1 = {}, c2 = {}, t = {}, beta = {}".format(bits, c1,c2,t,beta))
                        timer.start()
                        ret = ritter(N, t,c1,c2,beta)
                        overall = timer.stop()
                        ret.append(overall)
                        logging.warning(
                            "Overall runtime: {} seconds".format(overall))
                        betas[beta] = ret
                        beta += beta_step
                    ts[t] = betas
                    t += t_step
                c2s[c2] = ts
                c2 += c2_step
            c1s[c1] = c2s
            c1 += c1_step
        timing[bits] = c1s
    with open("./timing/ritter/" + current_date+current_time + '.pkl', 'wb') as fp:
        pickle.dump(timing, fp)

def test():
    multiplier = 10**10
    prec = 5
    P = n_primes(30)
    Basis = IntegerMatrix.from_matrix(generate_basis(P,multiplier,prec))
    print("Before reduction:")
    print(Basis)
    BKZs(Basis)(20)
    print("After reduction:")
    print(Basis)


if __name__ == "__main__":
    # test()
    main_ritter()


"""
To do:
tune hyperparameters
try sieving techniques
get stats on: if reduction is not strong enough or problem with permutation
"""
