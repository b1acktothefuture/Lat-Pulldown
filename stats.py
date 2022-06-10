import pickle


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
    return minus_one, zero, one


def print_fac_times(timing):
    for bits in timing:
        print('** bit_size: {}'.format(bits))
        for alpha in timing[bits]:
            for c in timing[bits][alpha]:
                print("alpha: {}, c: {} = {} seconds".format(
                    round(alpha, 2), round(c, 2), round(timing[bits][alpha][c][-1], 2)))


def main():
    name = '2022-06-1012:42:42.549221'
    with open('./timing/' + name + '.pkl', 'rb') as fp:
        timing = pickle.load(fp)
    print_fac_times(timing)


if __name__ == "__main__":
    main()
