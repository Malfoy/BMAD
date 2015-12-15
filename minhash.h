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
void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous);


vector<minimizer> allHash(uint8_t k,const string& seq);
string randomSeq(uint length);
uint sketchUnorderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
uint sketchOrderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2);
vector<minimizer> minHash(uint H, uint k, const string& seq);
vector<minimizer> minHashPart(uint H, uint k,const string& seq, uint part);
double scoreFromAlignment(const string& seq1,const string& seq2);
unordered_set<minimizer> allKmerSet(uint k,const string& seq);
vector<minimizer> minHashPart2(uint H, uint8_t k, const string& seq, uint8_t part, const unordered_set<minimizer>& filter);
unordered_multimap<string,string> allKmerMapStranded(uint8_t k,const string& seq, uint8_t nuc);
unordered_multimap<string,string> minHashErrors(uint H, uint8_t k, const string& seq, uint8_t nuc);
uint sketchUnorderedComparisonError(const unordered_multimap<string, string>& map1, const unordered_multimap<string, string>& map2);
string mutate(string read,int n);
vector<minimizer> minHashGenomic(uint H, uint k, const string& seq, const unordered_set<minimizer>& filter);
vector<minimizer> allGenomicKmers(uint k,const string& seq,unordered_set <minimizer> set);
void minHash2(uint H, uint k, const string& seq, vector<minimizer>& previous);
vector<minimizer> allQuasiGenomicKmers(uint k,const string& seq,unordered_multimap<minimizer,minimizer> map,uint nuc);
unordered_multimap<minimizer,minimizer> allKmerMap(const char k,const string& seq, const  char nuc);
vector<minimizer> allQuasiGenomicKmers(uint k,const string& seq,unordered_multimap<minimizer,minimizer> map,uint nuc);
unordered_map<minimizer,vector<readNumber>> indexReadSet(const string& readFile, const uint k, const uint seedSize,unordered_multimap<minimizer,minimizer> map);
#endif
