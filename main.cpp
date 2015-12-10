#include "minhash.h"
#include "utils.h"
#include <fstream>



using namespace std;



int main(int argc, char ** argv){
	srand (time(NULL));
	uint H(1000);
	uint nuc(10);
	uint k(11),part(2);
	string readFile1("reads_virus_region_5000_0001.ref");
	string readFile2("readvirus.fa");
	vector<string> V1,Vref;
	Vref=getReads(readFile1, 1);
	string ref(Vref[0]);
	// cout<<"ref"<<ref.size()<<endl;
	V1=getReads(readFile2,100);
	for(int j(0);j<20;j+=4){
	uint i(0),res1(0),res2(0),res3(0),res4(0),res5(0),res6(0),res7(0),res8(0),res9(0),res10(0);
		k=11+j;
		cout<<"k:"<<k<<endl;
		for(;i+1<V1.size();i+=2){
			string seq1(V1[i]),seq2(V1[i+1]);
			if(seq1.size()<H or seq2.size()<H){continue;}

			// cout<<"Test minhash vanilla"<<endl;
			vector<minimizer> sketch1(minHash(H,k,seq1)),sketch2(minHash(H,k,seq2));
			res1+=(sketchUnorderedComparison(sketch1,sketch2));
			res2+=(sketchOrderedComparison(sketch1,sketch2));

			// cout<<"Test minhash with minimizer repartition"<<endl;
			vector<minimizer> sketch3(minHashPart(H,k,seq1,part)),sketch4(minHashPart(H,k,seq2,part));
			res3+=(sketchUnorderedComparison(sketch3,sketch4));
			res4+=(sketchOrderedComparison(sketch3,sketch4));

			unordered_set<minimizer> solidKmers(allKmerSet(k,ref));
			// cout<<"Test minhash with minimizer repartition and considering only the kmers that appear in ref"<<endl;
			vector<minimizer> sketch5(minHashGenomic(H,k,seq1,solidKmers)),sketch6(minHashGenomic(H,k,seq2,solidKmers));
			res5+=(sketchUnorderedComparison(sketch5,sketch6));
			res6+=(sketchOrderedComparison(sketch5,sketch6));

			// cout<<"Test the kmers that appear in seq3"<<endl;
			vector<minimizer> sketch7(allGenomicKmers(k,seq1,solidKmers)),sketch8(allGenomicKmers(k,seq2,solidKmers));
			res7+=(sketchUnorderedComparison(sketch7,sketch8));
			res8+=(sketchOrderedComparison(sketch7,sketch8));

			unordered_multimap<minimizer,minimizer> quasi(allKmerMap(k,ref,nuc));
			// cout<<"Test the kmers that quasi-appear in seq3"<<endl;
			vector<minimizer> sketch9(allQuasiGenomicKmers(k,seq1,quasi,nuc)),sketch10(allQuasiGenomicKmers(k,seq2,quasi,nuc));
			res9+=(sketchUnorderedComparison(sketch9,sketch10));
			res10+=(sketchOrderedComparison(sketch9,sketch10));

		}

		cout
		<<"minhash vanilla                  : "<<res1/i<<" "<<res2/i<<" "<<endl
		<<"minhash with minimizer reparition: "<<res3/i<<" "<<res4/i<<" "<<endl
		<<"minhash with only genomic kmers  : "<<res5/i<<" "<<res6/i<<" "<<endl
		<<"genomic kmers                    : "<<res7/i<<" "<<res8/i<<endl
		<<"quasi genomic kmer               : "<<res9/i<<" "<<res10/i<<endl<<endl;
	}




	return 0;
}
