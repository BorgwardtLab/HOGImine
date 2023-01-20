#ifndef bitsetll_h
#define bitsetll_h
#include <vector>

using namespace std;

class Bitset{
public:
    vector<long long> bits;
    int n = 0;

    Bitset(){}
    ~Bitset(){}
    Bitset(int n_){
        n = n_;
        bits = vector<long long>((n+63)>>6); //ceil(n/64)
    }


    inline void set(int i){
        bits[i>>6] |= 1ll<<(i & 63);
    }
    inline void reset(int i){
        bits[i>>6] &= ~(1ll<<(i & 63));
    }
    inline void flip(int i){
        bits[i>>6] ^= (1ll<<(i & 63));
    }

    inline bool get(int i){
        return bits[i>>6]&(1ll<<(i & 63));
    }


    inline Bitset bitwise_or(Bitset& s2){
        Bitset out(n);
        for(int i=0; i<bits.size(); i++)
            out.bits[i] = bits[i] | s2.bits[i];
        return out;
    }
    inline Bitset bitwise_and(Bitset& s2){
        Bitset out(n);
        for(int i=0; i<bits.size(); i++)
            out.bits[i] = bits[i] & s2.bits[i];
        return out;
    }

    inline int count(){
        int cnt = 0;
        for(int i=0; i<bits.size(); i++)
            cnt += __builtin_popcountll(bits[i]);
        return cnt;
    }

};

/*
inline Bitset bitwise_or(Bitset s1, Bitset s2){
    Bitset out(s1.n);
    for(int i=0; i<s1.bits.size(); i++)
        out.bits[i] = s1.bits[i] | s2.bits[i];
    return out;
}
*/

#endif
