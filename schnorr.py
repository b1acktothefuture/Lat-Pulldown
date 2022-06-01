from fpylll import IntegerMatrix, SVP
import math
import numpy as np
from Crypto.Util import number
import pickle
from progress.bar import Bar
from lattice import *


def primesfrom2to(n):
    sieve = np.ones(n//3 + (n % 6 == 2), dtype=bool)
    sieve[0] = False
    for i in range(int(n**0.5)//3+1):
        if sieve[i]:
            k = 3*i+1 | 1
            sieve[((k*k)//3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i & 1))//3::2*k] = False
    return np.r_[2, 3, ((3*np.nonzero(sieve)[0]+1) | 1)]


'''
list of all primes from 2 to floor(N^alpha)
'''


def prime_base(N, alpha):
    n = int(math.log(N)**alpha)
    return primesfrom2to(n)


'''
Checks if a given number is pt-smooth
'''


def is_smooth(x, P):
    y = x
    for p in P:
        while p.divides(y):
            y /= p
    return abs(y) == 1


'''
returns empty array if not pt-smooth
'''


def factorize_smooth(n, primes):
    if(n < 0):
        n = -n
    exponents = [0]*len(primes)
    for i in range(len(primes)):
        if(n == 1):
            return exponents
        if(n < primes[i]):
            return []
        while(n % primes[i] == 0):
            n /= primes[i]
            exponents[i] += 1
    return []


'''
rounding funtion
'''


def sr(x, prec):
    return round(x * (10 ** prec))


'''
generates basis matrix as per schnorr's original algorithm
'''


def generate_basis(prime_base: list, multiplier: int, prec: int) -> IntegerMatrix:

    n = len(prime_base)

    B = [[0]*(n+1) for _ in range(n)]
    for i in range(n):
        B[i][n] = sr(multiplier*math.log(prime_base[i]), prec)
        B[i][i] = sr(math.log(prime_base[i]), prec)

    # B = IntegerMatrix.from_matrix(B)
    return B


'''
given e vector, checks if s = u - v*N is pt-smooth, returns their factorization if true
'''


def relation(e: list, prime_base: list, N: int) -> list:
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

    S = [int(s < 0)] + factorize_smooth(s, prime_base)
    T = [0] + T

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

        B_reduced = bkz_reduce(Basis, 6)  # try tuning the block size
        e_reduced = cvp_babai(B_reduced, target)
        w = B_reduced.multiply_left(e_reduced)

        e = []
        for i in range(len(w)-1):
            assert w[i] % refs[i] == 0
            e.append(w[i]//refs[i])

        # implement a checker function without actually computing u and v to reject longer vectors, similar to SVP one in Ritter's paper.
        rel = relation(e, P, N)

        if(len(rel) == 0):
            continue

        key = rel[0]
        if key not in relations:
            relations[key] = rel[1:]
            bar.next()

    bar.finish()

    if(save):
        with open(str(N) + '.pkl', 'wb') as f:
            pickle.dump(relations, f)

    return relations


'''
transforms matrix into reduced row echelon form
source: https://rosettacode.org/wiki/Reduced_row_echelon_form#Python
'''


def ToReducedRowEchelonForm(M):
    if not M:
        return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i = r
        while M[i][lead] == 0:
            i ^= 1
            if i == rowCount:
                i = r
                lead ^= 1
                if columnCount == lead:
                    return
        M[i], M[r] = M[r], M[i]
        lv = M[r][lead]
        M[r] = [mrx * lv for mrx in M[r]]
        for i in range(rowCount):
            if i != r:
                lv = M[i][lead]
                M[i] = [iv - lv*rv for rv, iv in zip(M[r], M[i])]
        lead += 1


'''
solves set of equations mod 2
'''


def solve_f2(equations):
    # reduce to row echelon form
    # find if solutions exit
    # if infinite, use back tracking to find a solution
    ToReducedRowEchelonForm(equations)
    pass


def main():
    alpha = 1.6
    c = 1.1  # C should be really small

    bits = 20
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    print("N: {} = {}*{}".format(N, p, q))
    schnorr(N, alpha, c, 5, False, True)


def test():
    M = [[0, 1, 0], [1, 1, 1], [0, 0, 1]]
    ToReducedRowEchelonForm(M)
    print(M)
    pass


if __name__ == "__main__":
    # test()
    main()
