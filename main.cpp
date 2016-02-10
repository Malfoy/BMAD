#include "minhash.h"
#include "utils.h"
#include <fstream>


using namespace std;


int main(int argc, char ** argv){
	srand (time(NULL));
	string input("reads.fa"),seq;
	ifstream in(input);
	uint k(15),seedSize(5),threshold(3);
	cout<<"Counting kmers"<<endl;
	unordered_map<minimizer,uint64_t > count(kmerCounting(input, k));
	cout<<"Compute qmers from kmers"<<endl;
	auto quasi(getSolidMap(count,2,k, seedSize));
	cout<<"Indexing reads"<<endl;
	unordered_map<minimizer,vector<readNumber>> qmer2reads;
	vector<string> reads;
	vector<minimizer> mins;
	uint nRead(0);
	while(not in.eof()){
		seq=getRead(in);
		if(not seq.empty()){
			reads.push_back(seq);
			mins=(allQuasiGenomicKmers(k,seq,quasi,seedSize));
			for(uint i(0);i<mins.size();++i){
				qmer2reads[mins[i]].push_back(nRead);
			}
			++nRead;
		}
	}
	cout<<"Querying reads"<<endl;
	vector<readNumber> toAdd;
	unordered_map<readNumber,uint> connivence;
	for(uint i(0);i<reads.size();++i){
		connivence={};
		seq=reads[i];
		cout<<i<<endl;
		cout<<"has ";
		mins=(allQuasiGenomicKmers(k,seq,quasi,seedSize));
		for(uint ii(0);ii<mins.size();++ii){
			toAdd=qmer2reads[mins[ii]];
			for(uint iii(0);iii<toAdd.size();++iii){
				++connivence[toAdd[iii]];
			}
		}
		uint count(0);
		for ( auto it = connivence.cbegin(); it != connivence.cend(); ++it ){
			if(it->second>threshold){
				count++;
			}
		}
		cout<<count<<" friends"<<endl;
	}



	return 0;
}
