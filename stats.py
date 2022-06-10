import pickle
from prettytable import PrettyTable

'''
structure

timing = {
bit_size:
    alpha:
        c:
            prime_basis generation time,
            [
                time to generate lattice basis,
                time to generate n+2 fac-relations,
                [
                    [time for each reduction],
                    [result of the reduction {1: not short enough, 0: found , -1: repeated}]
                ]

            ],
            time to solve linear eqautions,
            factor found {0 no / 1 yes},
            overall runtime
}
'''


def classify(data):
    l = len(data)
    minus_one = 0
    one = 0
    zero = 0
    for i in range(l):
        if data[i] == -1:
            minus_one += 1
        elif data[i] == 1:
            one += 1
        else:
            zero += 1
    minus_one /= l
    one /= l
    zero /= l
    return minus_one, zero, one, l


def print_stats(timing):
    t = PrettyTable(['bits', 'alpha', 'c', 'run time (seconds)', 'trials', 'not short',
                    '**found**', 'repeated'])
    for bits in timing:
        for alpha in timing[bits]:
            for c in timing[bits][alpha]:
                [repeat, found, not_short, l] = classify(
                    timing[bits][alpha][c][1][2][1])
                t.add_row([bits, round(alpha, 2), round(c, 2),
                           round(timing[bits][alpha][c][-1], 2), l, round(
                    not_short, 2), round(found, 2), round(repeat, 2)])
        t.add_row(['', '', '', '', '', '', '', ''])
    print(t)


def main():
    name = '2022-06-1012:42:42.549221'
    with open('./timing/' + name + '.pkl', 'rb') as fp:
        timing = pickle.load(fp)
    print_stats(timing)


if __name__ == "__main__":
    main()
