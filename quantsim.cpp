#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <set>
#include <complex>
#include <assert.h>
#include <numbers>
#include <cmath>
#include <limits>
#include <random>
#include "QuantumState.cpp"
#include "quantumAlgorithms.cpp"

using namespace std;

typedef long double ld;
typedef long long ll;

int main(){
    srand(500);
    /*
    QuantumState q(2);
    q.printState();
    q.H(0);
    q.H(1);
    q.printState();
    q.CZ({0}, 1);
    q.printState();
    q.H(0);
    q.H(1);
    q.printState();*/

    
    vector<int> f(16, 0);
    f[9] = 1;
    cout << Grovers(4, f) << endl;

    cout << Shors(63) << endl;

    /*
    QuantumState q(2);
    q.probs[0] = sqrt((ld)1/3);
    q.probs[1] = sqrt((ld)1/3);
    q.probs[3] = sqrt((ld)1/3);
    q.printState();
    q.QFT();
    q.printState();

    QuantumState q2(2);
    q2.probs[0] = sqrt((ld)1/3);
    q2.probs[1] = sqrt((ld)1/3);
    q2.probs[3] = sqrt((ld)1/3);
    q2.QFTraw();
    q2.printState();
    */
    

    //cout << periodQuantum(8, 21) << endl;;
    
    return 0;
}