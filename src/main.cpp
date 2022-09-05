#include "Factoring.h"
#include "OtherTests.h"
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    cout << "You have entered " << argc
         << " arguments:" << "\n";
    long N_size = atol(argv[1]);
    cout << "\nN ~ 10^" << N_size;

    unsigned long n = stoul(argv[2],nullptr);
    cout << "\nn is the dimension of the lattice (default is n = 90) : " << n;

    double c1 = stod(argv[3]);
    RR c(c1);
    cout << "\nc is a parameter of the lattice (default c = 0.714): " << c;

    long s = stol(argv[4]);
    cout << "\ns is the pruning level (default s = 14): " << s;

    double reduce_ratio = stod(argv[5]);
    cout << "\nreduce_ratio: bypass not-reducing the minimal distance if it is lower than reduce_ratio*current minimal distance  (0 = never bypass, 1 = always bypass, default = 0.85): " <<reduce_ratio;

    long bkz_strong = stol(argv[6]);
    long bkz_slight = stol(argv[7]);
    cout << "\nBlock sizes of the strong (default 32) and slight (default 20) BKZ-reduce:\n\n";
    cout << "strong BKZ block size: " << bkz_strong;
    cout << "\nslight BKZ block size: " << bkz_slight;

    double A_factor = stod(argv[8]);
    cout << "\nA_factor controls the maximum start value of A = A_factor * 0.25 * sum(r_ii^2) (default alpha = 0.2): " << A_factor;

    unsigned long min_eqns = stoul(argv[9],nullptr);
    cout << "\nHow many equations should be found before stopping the program and calculating statistics\n\n"
            "Min. equations = " << min_eqns;
    
    long seed_type = stol(argv[10]);
    cout << "\nChoose a seed for the random number generator:\n"
            "-2 = timestamp (default)\n"
            "-1 = random device (sometimes buggy)\n"
            "any non negative number: this number is used as seed\n\n"
            "Seed type: " << seed_type;

    char cf = argv[11][0];
    cout << "\nUse continued fractions (y/n) (defaults to yes): " << cf;

    int scalingType = stoi(argv[12]);
    cout << "\nChoose the scaling type:\n"
            "0 = mixed (default)\n"
            "1 = every row with probability 1/2\n"
            "2 = every row with probability 1/4\n"
            "3 = every row with probability 3/4\n"
            "4 = first n/2 rows with prop. 1/4, other with prop. 1/2\n"
            "5 = first n/2 rows with prop. 1/2, other with prop. 1/4\n"
            "Scaling type: " << scalingType;

    cout << "\n\n";
    // N n c s reduce strong slight A_fac min seed cf scal
    // ./main 14 90 0.714 14 0.8 32 20 0.2 91 -1 y 0
    Timer timer = Timer();
    Factoring(FactoringSettings(getN(N_size), n, c, s, A_factor, reduce_ratio, 10000, bkz_strong, bkz_slight, min_eqns, seed_type,cf =='y',scalingType));
    int time = timer.stopTimer();

    cout << "Time taken: " << time << endl;
    
    std::ofstream outfile;
    outfile.open("timing.txt", std::ios_base::app); // append instead of overwrite
    cout <<"Here";
    outfile << time << "\n"; 
    // Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 91,1448726806,true,0));
    return 0;
}