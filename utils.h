#ifndef UT
#define UT
#include <string>
#include <iostream>
#include <cstdint>
#include <cctype>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

using namespace std;

typedef uint64_t minimizer;
typedef unsigned int uint;


vector<string> getReads(string& ReadFile,uint n);
uint64_t nuc2int(char c);
string reversecomplement (const string& s);
void int2seq(minimizer min, uint n);
minimizer cat(uint32_t seed, uint32_t body, uint n);
minimizer getRepresent(minimizer min, uint n);
minimizer rc(minimizer min, uint n);

#endif
