#include "minhash.h"
#include "utils.h"
#include <fstream>



using namespace std;



int main(int argc, char ** argv){
	srand (time(NULL));
	ifstream refs("reads_virus_10k_0001.ref");
	ifstream readFile("reads.fa");
	string ref(getReads(refs,1)[0]);
	uint seedSize(8);
	vector<string> V1=getReads(readFile,10);
	// unordered_multimap<minimizer,minimizer> map(allKmerMap(k,ref, seedSize));
	// unordered_map<minimizer,vector<readNumber>> min2read(indexReadSet(readFile,k,seedSize,map));
	for(uint k(11);k<27;k+=2){
		uint res(0),res2(0),res3(0),res4(0);
		unordered_multimap<minimizer,minimizer> quasi(allKmerMap(k,ref,seedSize));
		unordered_set<minimizer> solidKmers(allKmerSet(k,ref));
		for(uint i(0);i+1<V1.size();i+=1){
			string seq1(V1[i]),seq2(V1[i+1]);
			vector<minimizer> sketch9(allQuasiGenomicKmers(k,seq1,quasi,seedSize)),sketch10(allQuasiGenomicKmers(k,seq2,quasi,seedSize));
			res+=(sketchUnorderedComparison(sketch9,sketch10));
			vector<minimizer> sketch7(allGenomicKmers(k,seq1,solidKmers)),sketch8(allGenomicKmers(k,seq2,solidKmers));
			res4+=(sketchUnorderedComparison(sketch7,sketch8));
		}
		// cout<<res<<endl;
		unordered_map<minimizer,uint8_t> count(kmerCounting("reads.fa", k));
		// cout<<"coutning"<<endl;
		quasi=getSolidMap(count,2,k, seedSize);
		unordered_set<minimizer>  solidKmers2=getSolidSet(count,2);
		for(uint i(0);i+1<V1.size();i+=1){
			string seq1(V1[i]),seq2(V1[i+1]);
			vector<minimizer> sketch9(allQuasiGenomicKmers(k,seq1,quasi,seedSize)),sketch10(allQuasiGenomicKmers(k,seq2,quasi,seedSize));
			vector<minimizer> sketch7(allGenomicKmers(k,seq1,solidKmers2)),sketch8(allGenomicKmers(k,seq2,solidKmers2));
			res3+=(sketchUnorderedComparison(sketch7,sketch8));
			res2+=(sketchUnorderedComparison(sketch9,sketch10));
		}
		cout<<"k: "<<k<<endl
		<<"intersection "<<interSet(solidKmers, solidKmers2)<<" missed genomic kmers "<<100*(double)inANotInB(solidKmers, solidKmers2)/solidKmers.size()
		<<"% kmer above threshold but not genomic "<<100*(double)inANotInB(solidKmers2, solidKmers)/solidKmers2.size()<<"%"<<endl
		<<"quasi genomic kmer with ref "<<res/V1.size()<<endl
		<<"genomic kmer with ref "<<res4/V1.size()<<endl
		<<"quasi genomic kmer without ref "<<res2/V1.size()<<endl
		<<"genomic kmer without ref "<<res3/V1.size()<<endl
		<<endl;
	}

	return 0;
}
