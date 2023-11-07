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

using namespace std;

typedef long double ld;
typedef long long ll;
typedef complex<ld> comp;

struct QuantumState{
    int n;
    vector<comp> probs;
    default_random_engine generator;
    uniform_real_distribution<ld> distribution;

    QuantumState(int n): n(n){
        probs = vector<comp>(1<<n, comp(0.0, 0.0));
        probs[0] = comp(1.0, 0.0);
        distribution = uniform_real_distribution<ld>(0.0, 1.0);
    }

    void I(){ return; }

    void X(int q){
        assert(q<n);
        for(int i = 0; i < (1 << n); i++){
            if((i>>q)&1) continue;
            int alt = i+(1<<q);
            swap(probs[i], probs[alt]);
        }
    }

    void CX(vector<int> cnt, int q){
        assert(q<n);
        for(int i = 0; i < (1 << n); i++){
            bool pass = true;
            for(int qq : cnt){
                assert(qq<n);
                assert(qq!=q);
                pass &= (i>>qq)&1;
            }
            if(!pass) continue;
            if((i>>q)&1) continue;
            int alt = i+(1<<q);
            swap(probs[i], probs[alt]);
        }
    }

    void H(int q){
        assert(q<n);
        ld s2 = sqrt(2);
        vector<comp> probs2(1<<n);
        for(int i = 0; i < (1 << n); i++){
            if((i>>q)&1){
                int z = i^(1<<q);
                int o = i;
                probs2[z] += probs[i]/s2;
                probs2[o] -= probs[i]/s2;
            }
            else{
                int z = i;
                int o = i^(1<<q);
                probs2[z] += probs[i]/s2;
                probs2[o] += probs[i]/s2;
            }
        }
        for(int i = 0; i < (1<<n); i++) probs[i] = probs2[i];

        normalize();
    }

    void Z(int q){
        assert(q < n);
        for(int i = 0; i < (1<<n); i++){
            if((i>>q)&1){
                probs[i] = probs[i]*(ld)(-1.0);
            }
        }
        normalize();
    }

    void CZ(vector<int> cnt, int q){
        assert(q<n);
        for(int i = 0; i < (1 << n); i++){
            bool pass = true;
            for(int qq : cnt){
                assert(qq<n);
                assert(qq!=q);
                pass &= (i>>qq)&1;
            }
            if(!pass) continue;
            if((i>>q)&1){
                probs[i] = probs[i]*(ld)(-1.0);
            }
        }
        normalize();
    }

    void normalize(){
        ld mag = 0.0;
        for(int i = 0; i < (1 << n); i++){
            mag += norm(probs[i]);
        }
        mag = sqrt(mag);
        ld a = abs(probs[0]);
        comp c = conj(probs[0])/a;
        if(a < 1e-12) c = comp(1.0, 0.0);

        for(int i = 0; i < (1 << n); i++){
            probs[i] *= c;
        }
    }

    int M(){
        ld p = distribution(generator);
        int i = 0;
        while(i < probs.size() && p > norm(probs[i])) {
            i++;
            p -= norm(probs[i]);
        }
        if(i == probs.size()) return i-1;
        return i;
    }

    void printDistribution(){
        for(int i = 0; i < (1<<n); i++){
            printf("%d: %.6f\n", i, (double)norm(probs[i]));
        }
    }
    void printState(){
        for(int i = 0; i < (1<<n); i++){
            cout << i << ": " << probs[i] << endl;
        }
    }
};