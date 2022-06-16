from fpylll import *
from sympy import Integer
from numpy import linalg as LA
import time


def svp():
    pass


def lll_reduce(B: list, delta: float = 0.99) -> IntegerMatrix:
    B_t = IntegerMatrix.from_matrix(B)
    LLL.reduction(B_t, delta)
    assert LLL.is_reduced(B_t, delta)
    return B_t


def bkz_reduce(A: list, block_size: int, prune: bool = True) -> IntegerMatrix:
    B = IntegerMatrix.from_matrix(A)
    if(prune):
        param = BKZ.Param(block_size=block_size,
                          strategies=BKZ.DEFAULT_STRATEGY)
    else:
        param = BKZ.Param(block_size=block_size)

    bkz_reduced = BKZ.reduction(B, param)
    return bkz_reduced


def cvp_babai(B: IntegerMatrix, t: list) -> list:
    M = GSO.Mat(B)
    _ = M.update_gso()
    w = M.babai(t)
    return list(w)


def l2_norm(x: list, y: list) -> float:
    assert len(x) == len(y)
    ret = 0
    for i in range(len(x)):
        ret += (x[i] - y[i])**2
    return ret**0.5


def test():
    '''
    Calling Babai's CVP on different basis of same Lattice:
    1. Default
    2. LLL reduced with delta = 0.99
    3. BKZ reduced

    // Sometimes BKZ gives worse norms
    '''

    FPLLL.set_random_seed(time.time())

    dim = 50
    mat = IntegerMatrix(dim, dim)
    mat.randomize("uniform", bits=20)

    A = mat[:-1]
    t = mat[-1]

    # print("Basis is: \n{}".format(A))
    print("Targer vector: \n{}".format(t))

    t1 = cvp_babai(A, t)
    print("Default: {}".format(l2_norm(A.multiply_left(t1), t)))

    B = lll_reduce(A)
    t2 = cvp_babai(B, t)
    print("LLL: {}".format(l2_norm(B.multiply_left(t2), t)))

    for i in range(2, 30):
        C = bkz_reduce(B, i)
        t3 = cvp_babai(C, t)
        print("BKZ (blocksize: {}): {}".format(
            i, l2_norm(C.multiply_left(t3), t)))

    t4 = CVP.closest_vector(A, t)
    print("CVP: {}".format(l2_norm(t4, t)))


def test_bkz():
    FPLLL.set_random_seed(time.time())

    dim = 20
    mat = IntegerMatrix(dim, dim)
    mat.randomize("uniform", bits=20)
    C = bkz_reduce(mat, 5)
    for i in C:
        print(i)
    pass


if __name__ == "__main__":
    test_bkz()
