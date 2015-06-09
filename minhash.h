#include <string>
#include <iostream>
#include <cstdint>
#include <cctype>
#include <vector>
#include <unordered_set>

typedef uint32_t minimizer;
typedef unsigned int uint;


using namespace std;

minimizer seq2int(const string& seq);
minimizer seq2intStranded(const string& seq);
void updateMinimizer(minimizer&	min, char nuc,size_t k);
void updateMinimizerRC(minimizer&	min, char nuc,size_t k);
void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous);
