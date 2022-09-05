#ifndef PRINT_H
#define PRINT_H

#include <vector>
#include <NTL/matrix.h>
#include <iostream>
#include <fstream>

using namespace NTL;
using namespace std;


template <typename T>
void print(Mat<T> A){
    long n = A.NumRows(), m = A.NumCols();
    cout <<endl;
    for(long i = 0;i<n;i++){
        for(long j = 0;j<m;j++) 
            cout << A[i][j] << " ";
        cout <<endl;
    }
    cout << endl;
}

template <typename T>
void print(Vec<T> A){
    long n = A.length();
    cout << endl;
    for(long i = 0;i<n;i++) cout << A[i] << " ";
    cout << endl;
}

template <typename T>
void write(Mat<T> A,ofstream& f){
    long n = A.NumRows(), m = A.NumCols();
    f <<endl;
    for(long i = 0;i<n;i++){
        for(long j = 0;j<m;j++) 
            f << A[i][j] << " ";
        f <<endl;
    }
    f << endl;
}

template <typename T>
void write(Vec<T> A,ofstream& f){
    long n = A.length();
    f << endl;
    for(long i = 0;i<n;i++) f << A[i] << " ";
    f << endl;
}
#endif