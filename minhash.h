#ifndef MH
#define MH

#include <string>
#include <iostream>
#include <cstdint>
#include <cctype>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>


typedef uint32_t minimizer;
typedef unsigned int uint;


using namespace std;


minimizer seq2int(const string& seq);
minimizer seq2intStranded(const string& seq);
void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous);


vector<minimizer> minHash(uint16_t H, uint8_t k, const string& seq);
vector<minimizer> allHash(uint8_t k,const string& seq);
string randomSeq(uint32_t length);
uint32_t sketchUnorderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
uint32_t sketchOrderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
vector<minimizer> minHash(uint32_t H, uint8_t k, const string& seq);
vector<minimizer> minHashPart(uint32_t H, uint8_t k,const string& seq, uint8_t part);
double scoreFromAlignment(const string& seq1,const string& seq2);
unordered_set <minimizer> allKmerSet(uint8_t k,const string& seq);
vector<minimizer> minHashPart2(uint32_t H, uint8_t k, const string& seq, uint8_t part, const unordered_set<minimizer>& filter);
unordered_multimap<string,string> allKmerMapStranded(uint8_t k,const string& seq, uint8_t nuc);
unordered_multimap<string,string> minHashErrors(uint32_t H, uint8_t k, const string& seq, uint8_t nuc);
uint32_t sketchUnorderedComparisonError(const unordered_multimap<string, string>& map1, const unordered_multimap<string, string>& map2);
string mutate(string read,int n);

#endif
