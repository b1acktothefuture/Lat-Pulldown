import pickle
import json
from schnorr import *
import random


def primesfrom2to_test():
    '''
    Test the primesfrom2to function
    '''
    def SieveOfEratosthenes(n):

        # Create a boolean array
        # "prime[0..n]" and initialize
        #  all entries it as true.
        # A value in prime[i] will
        # finally be false if i is
        # Not a prime, else true.
        prime = [True for i in range(n+1)]
        p = 2
        while (p * p <= n):

            # If prime[p] is not
            # changed, then it is a prime
            if (prime[p] == True):

                # Update all multiples of p
                for i in range(p * p, n+1, p):
                    prime[i] = False
            p += 1
        ret = []
        # Print all prime numbers
        for p in range(2, n+1):
            if prime[p]:
                ret.append(p)
        return ret

    for i in range(1000):
        n = random.randint(5, 100000)
        ret_test = primesfrom2to(n)
        print("testing primesfrom2to(%d)..." % n)
        assert(ret_test == SieveOfEratosthenes(n-1))


def factorize_smooth_test(limit=1000, trials=100, exp_limit=10):
    p = primesfrom2to(limit)
    for i in range(trials):
        exps = []
        num = 1
        for p_i in p:
            exps.append(random.randint(0, exp_limit))
            num *= p_i ** exps[-1]
        print("testing factorize_smooth({}...)".format(str(num)[:10]))
        assert(factorize_smooth(num, p) == exps)


def ritter_basis(N, P, c, prec=10, independent=False):
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
    Basis.insert(0, target)
    B_t = IntegerMatrix.from_matrix(Basis)
    print(B_t)
    pass


def main():
    bits = 20
    p = number.getPrime(bits//2)
    q = number.getPrime(bits//2)
    N = p*q

    print("N: {} = {}*{}".format(N, p, q))

    alpha = 1.7
    c = 1.2  # C should be really small
    prec = 10
    P = prime_base(N, alpha)

    ritter_basis(N, P, c, prec)


if __name__ == "__main__":
    main()
