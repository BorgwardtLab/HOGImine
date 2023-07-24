#ifndef UTILS_H
#define UTILS_H

#include <bits/stdc++.h>
#include "bitsetll.hpp"
#include "chi2.h"

#define hash_map std::unordered_map
#define hash_set std::unordered_set

using namespace std;

struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

#endif
