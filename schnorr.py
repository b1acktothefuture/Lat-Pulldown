from fpylll import IntegerMatrix, SVP
import math
import numpy as np
from Crypto.Util import number
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


def prime_base(N, alpha):
    n = int(math.log(N)**alpha)
    return primesfrom2to(n)


def is_smooth(x, P):
    y = x
    for p in P:
        while p.divides(y):
            y /= p
    return abs(y) == 1


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


def generate_basis(prime_base: list, multiplier: int, c: float, prec: int) -> IntegerMatrix:
    def sr(x):
        return round(x * (10 ** prec))
    n = len(prime_base)

    B = [[0]*(n+1) for _ in range(n)]
    for i in range(n):
        B[i][n] = sr(multiplier*math.log(prime_base[i]))
        B[i][i] = sr(math.log(prime_base[i]))

    B = IntegerMatrix.from_matrix(B)
    return B


def fac_relation(N, P, c, prec=1000):
    n = len(P)
    multiplier = N**c

    B = generate_basis(P, multiplier, c, prec)

    s = list(lll_reduce(B)[0])
    # Use approximate svp, stronger basis reduction, babai's nearest plane

    print(s)
    assert abs(s[0]) == 1

    if(s[0] == 1):
        for i in range(len(s)):
            s[i] = -1*s[i]

    e = [s[i+1] / (round((math.log(P[i])) * (10 ** prec))) for i in range(n)]

    u = 1
    v = 1
    a = [0]*(len(P)+1)

    for i in range(n):
        # assert e[i] in ZZ
        if e[i] > 0:
            a[i+1] = e[i]
            u *= P[i]**e[i]
        if e[i] < 0:
            v *= P[i]**(-1*e[i])

    b = [0] + factorize_smooth(u - v*N, P)
    if(len(b) == 1):
        return []

    if(u < v*N):
        b[0] = 1

    print("Success: ", u, u-v*N)
    return [a, b]


def factor(N, alpha, c, prec=1000):
    P = prime_base(N, alpha)
    n = len(P)

    for i in P:
        if(N % i == 0):
            return [i, N/i]

    print("Dimension: {}".format(n))

    return fac_relation(N, P, c, 3)


def main():
    alpha = 1.8
    c = 1.6

    bits = 40
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    print("N: {} = {}*{}".format(N, p, q))
    success = factor(N, alpha, c, 3)


def test():
    alpha = 1.8
    c = 1.6

    bits = 20
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    P = prime_base(N, alpha)
    multiplier = N**c

    print(generate_basis(P, multiplier, c, 5))


if __name__ == "__main__":
    test()
