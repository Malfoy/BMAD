#include "minhash.h"
#include "utils.h"
#include <fstream>



using namespace std;



int main(int argc, char ** argv){
	srand (time(NULL));
	string ref,readFile;
	uint k(25),seedSize(10);
	unordered_multimap<minimizer,minimizer> map(allKmerMap(k,ref, seedSize));
	unordered_map<minimizer,vector<readNumber>> min2read(indexReadSet(readFile,k,seedSize,map));

	return 0;
}
