#include "minhash.h"



using namespace std;



int main(){
	srand (time(NULL));
	uint32_t H(1000);
	uint8_t k(11);
	cout<<"Test"<<endl;
	//random sequences generated
	string seq1(randomSeq(10000)),seq2(randomSeq(10000));
	//We compute the sketches
	vector<minimizer> sketch1(minHash(H,k,seq1)),sketch2(minHash(H,k,seq2));
	//we compare the sketches
	uint32_t res1(sketchHammingComparison(sketch1,sketch2));
	uint32_t res2(sketchComparison(sketch1,sketch2));
	cout<<"hamming comparison : "<<res1<<" set comparison : "<<res2<<endl;

	return 0;
}
