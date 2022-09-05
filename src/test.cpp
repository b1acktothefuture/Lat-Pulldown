#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "GSA.h"
#include "print.h"
// g++ -o test test.cpp  -lntl -lgmp -lpthread
int main(){
    ofstream File("./GSA/GSA_d.txt");
    File.clear();
    
    
    long n = 200, accuracy_factor = 1000, block_size = 20;
    ZZ N;
    conv(N,"100000980001501");
    RR c(0.714);

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