#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "GSA.h"
#include "Print.h"
// g++ -o test test.cpp  -lntl -lgmp -lpthread
// ./test 100 1000 20 100000980001501 0.714
int main(int argc, char** argv){

    ofstream File("./GSA/GSA_d.txt");
    File.clear();
    long n = atol(argv[1]), accuracy_factor = atol(argv[2]), block_size = atol(argv[3]);

    // long n = 200, accuracy_factor = 1000, block_size = 20;
    // conv(N,"100000980001501");
    // conv(N, argv[4]);
    
    ZZ N(getN(atol(argv[4])));

    double _c = stod(argv[5]);
    RR c(_c);

    File << n << " " << accuracy_factor << " " << block_size <<" " << N << "\n" << c;

    GSA t(n,accuracy_factor,N,c);
    write(t.B,File);

    t.reduceBasis(block_size);
    t.GSO();

    write(t.B,File);
    write(t.U,File);
    write(t.sizes,File);

    File.close();
    
}