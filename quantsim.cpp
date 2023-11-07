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
const ld PI = 3.141592653589793268462643383279502884197169399375105820974944;
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

int period(int a, int N){
    int res = 1;
    int p = a;
    while(p != 1){
        //if(res%1000 == 0) cout << res << endl;
        p = (int)((ll)p*a%N);
        res++;
    }
    return res;
}

int periodQuantum(int a, int N){
    int n = 0;
    int p = 1;
    while(p <= N){
        p *= 2;
        n++;
    }
    QuantumState q(3*n);
}
int gcd(int a, int b){
    if(a == 0 || b == 0) return a+b;
    return gcd(b, a%b);
}
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
        int r = period(a, N);
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

    /*
    vector<int> f(16, 0);
    f[9] = 1;
    cout << Grovers(4, f) << endl;*/

    cout << Shors(176399) << endl;
    
    return 0;
}