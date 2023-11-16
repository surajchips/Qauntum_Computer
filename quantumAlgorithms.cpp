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

using namespace std;

typedef long double ld;
typedef long long ll;
typedef complex<ld> comp;

//apply oracle on qubit n with inputs [0, n-1]
void oracle(QuantumState &q, int n, vector<int> &f){
    vector<int> cnt;
    for(int i = 0; i < n; i++) cnt.push_back(i);
    for(int i = 0; i < (1<<n); i++){
        if(f[i]){
            for(int e = 0; e < n; e++){
                if(!((i>>e)&1)){
                    q.X(e);
                }
            }
            q.CX(cnt, n);
            for(int e = 0; e < n; e++){
                if(!((i>>e)&1)){
                    q.X(e);
                }
            }
        }
    }
}

//Deutsch-Jozsa algorithm
//returns 0 if f(x) is constant, returns 1 if f(x) is balanced
bool DJ(int n, vector<int> &f){
    QuantumState q(n+1);
    q.X(n);
    for(int i = 0; i <= n; i++){
        q.H(i);
    }

    //q.printState();

    oracle(q, n, f);

    //q.printState();

    for(int i = 0; i <= n; i++){
        q.H(i);
    }
    q.X(n);

    //q.printState();
    int r = q.M();

    return (r!=0);
}

//Grovers algorithm
//returns x such that f(x)=1
int Grovers(int n, vector<int> f){
    QuantumState q(n);

    ld k = PI/4;
    k *= sqrt((ld)(1<<n));
    k -= 0.5;

    int rep = (int)(k+0.5);

    for(int i = 0; i < n; i++){
        q.H(i);
    }
    vector<int> cnt;
    for(int i = 1; i < n; i++) cnt.push_back(i);
    for(int i = 0; i < rep; i++){
        //oracle
        for(int i = 0; i < (1<< n); i++){
            if(!f[i]) continue;
            for(int e = 0; e < n; e++){
                if(!((i>>e)&1)){
                    q.X(e);
                }
            }
            q.CZ(cnt, 0);
            for(int e = 0; e < n; e++){
                if(!((i>>e)&1)){
                    q.X(e);
                }
            }
        }

        //amplification
        for(int e = 0; e < n; e++) q.H(e);
        for(int e = 0; e < n; e++) q.X(e);
        q.CZ(cnt, 0);
        for(int e = 0; e < n; e++) q.X(e);
        for(int e = 0; e < n; e++) q.H(e);
    }

    //q.printDistribution();
    //q.printState();

    return q.M();
}

//Brute force period finding algorithm
int period(int a, int N){
    int res = 1;
    int p = a;
    while(p != 1){
        p = (int)((ll)p*a%N);
        res++;
    }
    return res;
}

/// @brief find periods using shors algorithm
/// @param a base
/// @param N mod
/// @return period of function
int periodQuantum(int a, int N){
    int n = 0;
    int pw = 1;
    while(pw <= N){
        pw *= 2;
        n++;
    }
    int reg1 = 2*n;
    QuantumState qs(reg1+n);
    for(int i = 0; i < reg1; i++){
        qs.H(i);
    }
    qs.X(reg1);
    qs.ModExp(a, N, reg1);
    vector<int> v;
    for(int i = 0; i < n; i++) v.push_back(reg1+i);
    qs = qs.M(v).first;
    qs.QFT();

    //cheat: if 0 is found, just try again
    int res = qs.M();
    while(res == 0) res = qs.M();

    //continued fractions
    v.clear();
    int p = res;
    int q = 1<<(2*n);
    while(p != 0){
        v.push_back((int)((ld)q/p+0.5));
        swap(p, q);
        p -= v.back()*q;
        p = abs(p);
    }
    int dem = 0;
    for(int e = 0; e < v.size(); e++){
        p = 1;
        q = v[e];
        for(int i = e-1; i >= 0; i--){
            p += q*v[i];
            swap(p, q);
        }
        if(q >= N){
            break;
        }
        dem = q;
    }

    //test all multiples of r0
    pw = 1;
    for(int i = 0; i < dem; i++) pw = (pw*a)%N;
    int p2 = 1;
    for(int res = dem; res < N; res += dem){
        p2 = (p2*pw)%N;
        if(p2 == 1) return res;
    }
    exit(1);
    return -1;
}
int gcd(int a, int b){
    if(a == 0 || b == 0) return a+b;
    return gcd(b, a%b);
}
//Shors algorithm
//Returns non-trivial factor of N, or -1 if none found after 100 trials
int Shors(int N){
    cout << "Factorizing " << N << endl;
    for(int tr = 1; tr < 100; tr++){
        cout << "Beginning iteration " << tr << endl;
        int a = rand()%(N-2)+2;

        printf("Randomly picked a=%d\n", a);
        if(gcd(a, N) != 1) {
            int g = gcd(a, N);
            printf("%d is a factor of %d\n", g, N);
            return g;
        }
        int r = periodQuantum(a, N);
        printf("Period of a is %d mod N\n", r);
        if(r%2 == 1) {
            cout << "Period of a is odd, next pick\n";
            continue;
        }
        r /= 2;
        int f = 1;
        for(int i = 0; i < r; i++) {
            f *= a;
            f%=N;
        }
        f+=N-1;
        if(gcd(f, N) != 1){
            int g = gcd(f, N);
            printf("%d is a factor of %d\n", g, N);
            return g;
        }
        cout << "Bad luck. Next try\n";
    }
    return -1;
}