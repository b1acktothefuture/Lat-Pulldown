import pickle
from progress.bar import Bar
import math
import numpy as np
from Crypto.Util import number
import galois
from lattice import *


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
            n = n//primes[i]
            exponents[i] += 1
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


def fac_relation(e: list, prime_base: list, N: int) -> list:
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


def schnorr(N, alpha, c, prec=10, independent=False, save=False):
    P = prime_base(N, alpha)
    n = len(P)

    for i in P:
        if(N % i == 0):
            return [i, N/i]

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
    bar = Bar('Finding relations', max=n+2)
    while(len(relations) < n+2):
        # assuming the probablity of getting same permuation more then once is very low
        np.random.shuffle(Basis)

        B_reduced = bkz_reduce(Basis, 4)  # try tuning the block size
        e_reduced = cvp_babai(B_reduced, target)
        w = B_reduced.multiply_left(e_reduced)

        e = []
        for i in range(len(w)-1):
            assert w[i] % refs[i] == 0
            e.append(w[i]//refs[i])

        # implement a checker function without actually computing u and v to reject longer vectors, similar to SVP one in Ritter's paper.
        rel = fac_relation(e, P, N)

        if(len(rel) == 0):
            # amend it, do not continue reduce it strongly
            continue

        key = rel[0]

        assert (key[0] - key[1]) % N == 0

        if key not in relations:
            relations[key] = rel[1:]
            bar.next()

    bar.finish()

    P.insert(0, -1)
    if(save):
        with open(str(N) + '.pkl', 'wb') as f:
            pickle.dump(relations, f)
        with open(str(N) + '_primes.pkl', 'wb') as f:
            pickle.dump(P, f)

    return [relations,  P]


def solve_linear_mod2(a_b):
    '''
    solves set of homogeneous equations mod 2
    As of now it returns only one solution
    '''
    a_plus_b_mod2 = []
    for i in range(len(a_b)):
        temp = []
        for j in range(len(a_b[i][0])):
            temp.append((a_b[i][0][j] + a_b[i][1][j]) % 2)
        a_plus_b_mod2.append(temp)
    F2 = galois.GF(2)
    A = F2(a_plus_b_mod2)
    t = A.left_null_space()
    return t[0]


def main():

    bits = 30
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    print("N: {} = {}*{}".format(N, p, q))

# ========================================================================

    alpha = 1.5
    c = 1.1  # C should be really small
    [relations, primes] = schnorr(N, alpha, c, 5, False, True)

    a_b = list(relations.values())
    sol = solve_linear_mod2(a_b)

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
    y = math.gcd(N, a+b)

    if(1 < x < N):
        print("Factor found: {}".format(x))
        return

    if(1 < y < N):
        print("Factor found: {}".format(y))
        return

    print("Better luck next time:(")


def test():
    N = 707309
    with open(str(N) + '.pkl', 'rb') as f:
        relations = pickle.load(f)

    with open(str(N) + '_primes.pkl', 'rb') as f:
        primes = pickle.load(f)

    a_b = list(relations.values())
    sol = solve_linear_mod2(a_b)

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

    print(len(primes))
    for i in range(len(A)):
        a *= primes[i]**int(A[i])
        b *= primes[i]**int(B[i])

    print("Are squares congruent: ", (a**2 - b**2) % N)


if __name__ == "__main__":
    # test()
    main()
