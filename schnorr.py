import pickle
from progress.bar import Bar
import logging
import math
import numpy as np
from Crypto.Util import number
import galois
from lattice import *


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
            logging.warn(
                "No non-trivial solutions exist to the system of Homogenoeus equations")
            self.has_been_zero = 1

    def bin_array(self, num, m):
        """Convert a positive integer num into an m-bit bit vector"""
        return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8).tolist()

    def next(self):
        sol = self.Fp([0]*self.dim)
        if self.has_been_zero == 1:
            logging.warn("No more solutions exist, returning empty array")
            return sol

        bits = self.bin_array(self.dp, len(self.t))
        for i in range(len(bits)):
            sol += bits[i]*self.t[i]
        self.dp = (self.dp + 1) % self.max
        if(self.dp == 0):
            self.has_been_zero = 1
        return sol


def get_primes(count):
    primes = []

    n = 2
    while len(primes) != count:
        for i in primes:
            if n % i == 0:
                break
        else:
            primes.append(n)
        n += 1

    return primes


def primesfrom2to(n):
    '''
    returns a list of primes from 2 to n-1
    '''
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
    n = int((math.log(N))**alpha)
    return primesfrom2to(n)


def is_smooth(x: int, P: list) -> bool:
    '''
    Checks if a given number is pt-smooth
    '''
    y = x
    for p in P:
        while y % p == 0:
            y //= p
    return abs(y) == 1


def factorize_smooth(n, primes):
    '''
    returns empty array if not pt-smooth
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
            n //= primes[i]
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
    '''

    n = len(prime_base)

    B = [[0]*(n+1) for _ in range(n)]
    for i in range(n):
        B[i][n] = sr(multiplier*math.log(prime_base[i]), prec)
        B[i][i] = sr(math.log(prime_base[i]), prec)

    # B = IntegerMatrix.from_matrix(B)
    return B


def check_fac_relation(e: list, prime_base: list, N: int) -> list:
    '''
    given e vector, checks if s = u - v*N is pt-smooth, returns their factorization if true
    '''
    assert len(e) == len(prime_base)
    [u, v] = [1, 1]
    T = [0]*len(prime_base)
    # print("Solution: {}".format(e))
    for i in range(len(e)):
        if(e[i] < 0):
            exp = -1*e[i]
            v *= prime_base[i]**(exp)
        else:
            u *= prime_base[i]**e[i]
            T[i] = e[i]

    s = u - v*N

    S = factorize_smooth(s, prime_base)
    S.insert(0, int(s < 0))
    T.insert(0, 0)

    if(len(S) == 1):
        return []
    return [tuple([u, s]), tuple(T), tuple(S)]


def fac_relations_cvp(N, P, c, prec=10, independent=False):

    n = len(P)
    logging.info("dimension: {}".format(n))

    if(independent):
        bit_length = N.bit_length()
        multiplier = 2**(bit_length*c)
    else:
        multiplier = N**c

    # Basis matrix for the prime number lattice
    Basis = generate_basis(P, multiplier, prec)
    refs = [Basis[i][i] for i in range(len(P))]

    # Target vector for CVP [0,..., N^c*log(N)]
    target = [0]*(len(P)+1)
    target[-1] = sr(multiplier*math.log(N), prec)

    relations = {}
    # Solve CVP here
    trial = 0
    bar = Bar('Finding relations', max=n+2)
    while(len(relations) < n+2):
        # assuming the probablity of getting same permuation more then once is very low
        trial += 1
        np.random.shuffle(Basis)

        # B_reduced = bkz_reduce(Basis, 6)  # try tuning the block size
        B_reduced = lll_reduce(Basis)
        e_reduced = cvp_babai(B_reduced, target)
        w = B_reduced.multiply_left(e_reduced)

        e = []
        for i in range(len(w)-1):
            assert w[i] % refs[i] == 0
            e.append(w[i]//refs[i])

        # implement a checker function without actually computing u and v to reject longer vectors, similar to SVP one in Ritter's paper.
        rel = check_fac_relation(e, P, N)

        if(len(rel) == 0):
            # amend it, do not continue reduce it strongly
            logging.info("Trial {}: not short enough".format(trial))
            continue

        key = rel[0]

        assert (key[0] - key[1]) % N == 0

        if key not in relations:
            logging.info(
                "Trial {}: found new fac-realtion: {}".format(trial, key))
            relations[key] = rel[1:]
            bar.next()

        else:
            logging.info("Trial {}: relation already exists".format(trial))

    bar.finish()

    P.insert(0, -1)

    return relations


def fac_relations_svp(N, P, c, prec=10, independent=False):
    n = len(P)
    logging.info("dimension: {}".format(n))
    if(independent):
        multiplier = 10**(c)
    else:
        multiplier = N**c

    # Basis matrix for the prime number lattice
    Basis = generate_basis(P, multiplier, prec)
    refs = [Basis[i][i] for i in range(len(P))]
    for i in range(len(Basis)):
        Basis[i].insert(0, 0)
    target = [0]*(len(P)+2)
    target[-1] = sr(multiplier*math.log(N), prec)
    target[0] = 1
    Basis.append(target)

    bar = Bar('Finding relations', max=n+3)
    trial = 0
    relations = {}
    while(len(relations) < n+3):
        trial += 1
        np.random.shuffle(Basis)

        B_reduced = bkz_reduce(Basis, 40)
        # B_reduced = lll_reduce(Basis)
        for vec in B_reduced:
            if(abs(vec[0]) != 1):
                continue
            e = []
            for i in range(1, n+1):
                assert vec[i] % refs[i-1] == 0
                e.append(vec[i]//refs[i-1])
            rel = check_fac_relation(e, P, N)
            if(len(rel) == 0):
                logging.info("Trial {}: not short enough".format(trial))
                continue
            key = rel[0]
            if key not in relations:
                logging.info(
                    "Trial {}: found new fac-realtion: {}".format(trial, key))
                relations[key] = rel[1:]
                bar.next()
            else:
                logging.info("Trial {}: relation already exists".format(trial))
            if(len(relations) >= n+3):
                break
    bar.finish()

    P.insert(0, -1)

    return relations


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
        logging.info("C = {}".format(sol))
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


def schnorr(N, alpha, c, prec, t=0):
    # should generate primes by fixing t(no. of primes)
    if(t > 0):
        P = get_primes(t)
    else:
        P = prime_base(N, alpha)

    logging.info("prime base: {}".format(P))

    relations = fac_relations_svp(N, P, c, prec, True)

    a_b = list(relations.values())
    fac = solve_linear(N, a_b, P)

# ========================================================================


def main():

    bits = 40
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    print("N: {} = {}*{}".format(N, p, q))

    alpha = 2.0
    c = 1.1  # C should be really small
    prec = 5
    open("./logs/" + str(N) + '.log', "w").close()
    logging.basicConfig(filename="./logs/" + str(N) + '.log',
                        encoding='utf-8', level=logging.INFO)

    logging.info(
        'N: {} = {} * {} \nalpha: {} \nc = {} \nprecision = {}'.format(N, p, q, alpha, c, prec))

    schnorr(N, alpha, c, prec)
    logging.basicConfig()


def ritter_test():
    N = 2131438662079
    t = 125
    c1 = 15  # multiplier
    c2 = 3  # prec
    open("./logs/" + str(N) + '.log', "w").close()
    logging.basicConfig(filename="./logs/" + str(N) + '.log',
                        encoding='utf-8', level=logging.INFO)

    logging.info(
        'N: {}\nt: {} \nc1 = {}\nc2 = {} '.format(N, t, c1, c2))

    schnorr(N, 0, c1, c2, t)
    logging.basicConfig()
    pass


def test():
    print(get_primes(100))


if __name__ == "__main__":
    ritter_test()
    # main()


"""
To do:
tune hyperparameters
try sieving techniques
get stats on: if reduction is not strong enough or problem with permutation
"""
