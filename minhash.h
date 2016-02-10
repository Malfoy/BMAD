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
#include "utils.h"




using namespace std;


minimizer seq2int(const string& seq);
minimizer seq2intStranded(const string& seq);
void minHash2(uint64_t H, uint64_t k, const string& seq, vector<minimizer>& previous);
vector<minimizer> allHash(uint64_t k,const string& seq);
string randomSeq(uint64_t length);
uint64_t sketchUnorderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
uint64_t sketchOrderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
vector<minimizer> minHash(uint64_t H, uint64_t k, const string& seq);
vector<minimizer> minHashPart(uint64_t H, uint64_t k,const string& seq, uint64_t part);
double scoreFromAlignment(const string& seq1,const string& seq2);
unordered_set<minimizer> allKmerSet(uint64_t k,const string& seq);
vector<minimizer> minHashPart2(uint64_t H, uint64_t k, const string& seq, uint64_t part, const unordered_set<minimizer>& filter);
unordered_multimap<string,string> allKmerMapStranded(uint64_t k,const string& seq, uint64_t nuc);
unordered_multimap<string,string> minHashErrors(uint64_t H, uint64_t k, const string& seq, uint64_t nuc);
uint64_t sketchUnorderedComparisonError(const unordered_multimap<string, string>& map1, const unordered_multimap<string, string>& map2);
string mutate(string read,int n);
vector<minimizer> minHashGenomic(uint64_t H, uint64_t k, const string& seq, const unordered_set<minimizer>& filter);
vector<minimizer> allGenomicKmers(uint64_t k,const string& seq,const unordered_set <minimizer>& set);
void minHash2(uint64_t H, uint64_t k, const string& seq, vector<minimizer>& previous);
vector<minimizer> allQuasiGenomicKmers(uint64_t k,const string& seq,const unordered_multimap<minimizer,minimizer>& map,uint64_t nuc);
unordered_multimap<minimizer,minimizer> allKmerMap(const uint64_t k,const string& seq, const  uint64_t nuc);
unordered_map<minimizer,uint64_t > kmerCounting(const string& readFile, const uint64_t k);
unordered_map<minimizer,vector<readNumber>> indexReadSet(const string& readFile, const uint64_t k, const uint64_t seedSize,const unordered_multimap<minimizer,minimizer>& map);
unordered_multimap<minimizer,minimizer> getSolidMap(unordered_map<minimizer,uint64_t >& count, uint64_t T, uint64_t k, uint64_t nuc);
unordered_set<minimizer> getSolidSet(unordered_map<minimizer,uint64_t >& count, uint64_t T);





#endif
