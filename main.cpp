#include "minhash.h"
#include "nw.h"



using namespace std;



int main(){
	srand (time(NULL));
	uint32_t H(100);
	uint8_t k(15),part(10);
	double errorRate(0.15);
	//random sequences generated
	string seq1(randomSeq(10000)),seq2(mutate(seq1,seq1.size()*errorRate)),seq3(mutate(seq1,seq1.size()*errorRate)),align1,align2;
	//~ cout<<seq1<<endl;
	//~ cout<<seq2<<endl;
	nw(seq1,seq2,align1,align2,false);
	//~ cout<<align1<<endl;
	//~ cout<<align2<<endl;
	cout<<"Alignment identity score : "<<scoreFromAlignment(align1,align2)<<endl;


	cout<<"Test minhash vanilla"<<endl;
	//We compute the sketches
	vector<minimizer> sketch1(minHash(H,k,seq1)),sketch2(minHash(H,k,seq2));
	//we compare the sketches
	uint32_t res1(sketchUnorderedComparison(sketch1,sketch2));
	uint32_t res2(sketchOrderedComparison(sketch1,sketch2));
	cout<<"Ordered comparison : "<<res1<<" unordered comparison : "<<res2<<endl;

	cout<<"Test minhash with minimizer repartition"<<endl;
	//We compute the sketches but we force a reparition of the minimizers
	vector<minimizer> sketch3(minHashPart(H,k,seq1,part)),sketch4(minHashPart(H,k,seq2,part));
	//we compare the sketches
	uint32_t res3(sketchUnorderedComparison(sketch3,sketch4));
	uint32_t res4(sketchOrderedComparison(sketch3,sketch4));
	cout<<"Ordered comparison : "<<res3<<" unordered comparison : "<<res4<<endl;

	unordered_set<minimizer> solidKmers(allKmerSet(k,seq3));
	cout<<"Test minhash with minimizer repartition and considering only the kmers that appear in seq3"<<endl;
	//We compute the sketches but we force a reparition of the minimizers
	vector<minimizer> sketch5(minHashPart2(H,k,seq1,part,solidKmers)),sketch6(minHashPart2(H,k,seq2,part,solidKmers));
	//we compare the sketches
	uint32_t res5(sketchUnorderedComparison(sketch5,sketch6));
	uint32_t res6(sketchOrderedComparison(sketch5,sketch6));
	cout<<"Ordered comparison : "<<res5<<" unordered comparison : "<<res6<<endl;


	cout<<"Test minhash with errors"<<endl;
	char nuc(3);

	unordered_multimap<string,string> map1 (minHashErrors(H,k,seq1,nuc));
	unordered_multimap<string,string> map2 (minHashErrors(H,k,seq2,nuc));

	uint32_t res7(sketchUnorderedComparisonError(map1,map2));
	cout<<"Unordered comparison  with error: "<<res7<<endl;

	return 0;
}
