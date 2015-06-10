#include "minhash.h"
#include "nw.h"



using namespace std;



int main(){
	srand (time(NULL));
	uint32_t H(100);
	uint8_t k(11),part(10);
	//random sequences generated
	string seq1(randomSeq(10000)),seq2(randomSeq(10000)),seq3(randomSeq(100000)),align1,align2;
	//~ cout<<seq1<<endl;
	//~ cout<<seq2<<endl;
	nw(seq1,seq2,align1,align2,false);
	//~ cout<<align1<<endl;
	//~ cout<<align2<<endl;
	cout<<"alignment score : "<<scoreFromAlignment(align1,align2)<<endl;


	cout<<"Test minhash vanilla"<<endl;
	//We compute the sketches
	vector<minimizer> sketch1(minHash(H,k,seq1)),sketch2(minHash(H,k,seq2));
	//we compare the sketches
	uint32_t res1(sketchHammingComparison(sketch1,sketch2));
	uint32_t res2(sketchComparison(sketch1,sketch2));
	cout<<"hamming comparison : "<<res1<<" set comparison : "<<res2<<endl;

	cout<<"Test minhash with minimizer repartition"<<endl;
	//We compute the sketches but we force a reparition od the minimizers
	vector<minimizer> sketch3(minHashPart(H,k,seq1,part)),sketch4(minHashPart(H,k,seq2,part));
	//we compare the sketches
	uint32_t res3(sketchHammingComparison(sketch3,sketch4));
	uint32_t res4(sketchComparison(sketch3,sketch4));
	cout<<"hamming comparison : "<<res3<<" set comparison : "<<res4<<endl;

	unordered_set<minimizer> solidKmers(allKmerSet(k,seq3));
	cout<<"Test minhash with minimizer repartition and considering only the kmers that appear in seq3"<<endl;
	//We compute the sketches but we force a reparition od the minimizers
	vector<minimizer> sketch5(minHashPart2(H,k,seq1,part,solidKmers)),sketch6(minHashPart2(H,k,seq2,part,solidKmers));
	//we compare the sketches
	uint32_t res5(sketchHammingComparison(sketch5,sketch6));
	uint32_t res6(sketchComparison(sketch5,sketch6));
	cout<<"hamming comparison : "<<res5<<" set comparison : "<<res6<<endl;

	return 0;
}
