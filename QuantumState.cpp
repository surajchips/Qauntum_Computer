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
#pragma once

using namespace std;

typedef long double ld;
typedef long long ll;
typedef complex<ld> comp;
const ld PI = 3.141592653589793268462643383279502884197169399375105820974944;

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

    void SWAP(int q1, int q2){
        assert(q1<n);
        assert(q2<n);
        for(int i = 0; i < (1 << n); i++){
            if((i>>q1)&1) continue;
            if(((i>>q2)&1) == 0) continue;
            int alt = i+(1<<q1)-(1<<q2);
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

    void Rk(int q, int k){
        assert(q < n);
        comp c = polar((ld)1.0, (ld)2*PI/(1<<k));
        for(int i = 0; i < (1<<n); i++){
            if((i>>q)&1){
                probs[i] *= c;
            }
        }
    }

    void CRk(vector<int> cnt, int q, int k){
        assert(q<n);
        comp c = polar((ld)1.0, (ld)2*PI/(1<<k));
        for(int i = 0; i < (1 << n); i++){
            bool pass = true;
            for(int qq : cnt){
                assert(qq<n);
                assert(qq!=q);
                pass &= (i>>qq)&1;
            }
            if(!pass) continue;
            if((i>>q)&1){
                probs[i] *= c;
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
        c /= mag;
        
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

    pair<QuantumState, int> M(vector<int> cnt){
        sort(cnt.begin(), cnt.end());
        reverse(cnt.begin(), cnt.end());
        int mbits = cnt.size();
        int qbits = n-cnt.size();
        vector<ld> mprobs(1<<mbits);
        for(int i = 0; i < (1<<n); i++){
            int mnum = 0;
            for(int c : cnt){
                mnum*=2;
                if((i>>c)&1){
                    mnum++;
                }
            }
            mprobs[mnum] += norm(probs[i]);
        }
        ld p = distribution(generator);
        int res = 0;
        while(res < mprobs.size() && p > mprobs[res]) {
            res++;
            p -= mprobs[res];
        }
        if(res == mprobs.size()) res--;

        QuantumState qs(qbits);
        qs.probs[0] = 0.0;

        int msk = 0;
        int mskr = 0;
        for(int i = 0; i < cnt.size(); i++){
            msk += (1<<cnt[i]);
            if((res>>i)&1) mskr += (1<<cnt[i]);
        }
        for(int i = 0; i < (1<<n); i++){
            if(((i^mskr)&msk) != 0) continue;
            int nmsk = i;
            for(int c : cnt){
                nmsk = ((nmsk>>(c+1))<<c)+(nmsk&((1<<c) - 1));
            }
            qs.probs[nmsk] = probs[i];
        }
        qs.normalize();
        return make_pair(qs, res);
    }

    //TO-DO: replace with gate call version
    void QFTraw(){
        vector<comp> probs2(1<<n);
        for(int x = 0; x < (1<<n); x++){
            for(int y = 0; y < (1<<n); y++){
                ld theta = 2*PI*x*y;
                theta /= (1<<n);
                probs2[y] += probs[x]*polar((ld)1.0, theta);
            }
        }
        for(int i = 0; i < (1<<n); i++) probs[i] = probs2[i];

        normalize();
    }

    //QFT implemented with quantum gates
    void QFT(){
        for(int i = 0; n-i-1 > i; i++){
            SWAP(i, n-i-1);
        }
        H(0);
        vector<int> v(1);
        for(int i = 1; i < n; i++){
            v[0] = i;
            CRk(v, 0, i+1);
        }
        for(int i = 1; i+1 < n; i++){
            H(i);
            v[0] = i+1;
            CRk(v, i, 2);
        }
        H(n-1);
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