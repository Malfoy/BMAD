#include "minhash.h"
#include "xor.h"
#include "utils.h"


using namespace std;


void updateMinimizer(minimizer&	min, char nuc, uint k){
	// cout<<"go"<<endl;
	// int2seq(min, 17);
	minimizer offset(1);
	offset<<=(2*k);
	// cout<<offset<<endl;
	// int2seq(offset, 18);
	min<<=2;
	// int2seq(min, 17);
	min+=nuc2int(nuc);
	// int2seq(min, 17);
	min%=offset;
	// int2seq(min, 17);
	// cout<<"end"<<endl;
}


void updateMinimizerEnd(minimizer&	min, char nuc, uint k){
	min>>=2;
	min+=(nuc2int(nuc)<<(2*k-2));
}


void updateMinimizerRC(minimizer&	min, char nuc, uint k){
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
}


//return a sketch containing all kmers
vector<minimizer> allKmer(uint k,const string& seq){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint i(0);
	do{
		sketch.push_back(kmer);
		if(i+k<seq.size()){
			updateMinimizer(kmerS, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
			kmer=min(kmerRC,kmerS);
		}else{
			return sketch;
		}
		++i;
	}while(true);
	return sketch;
}


//return a sketch containing all kmers
vector<minimizer> allGenomicKmers(uint k,const string& seq,unordered_set <minimizer> set){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint i(0);
	do{
		if(set.unordered_set::count(kmer)!=0){
			sketch.push_back(kmer);
		}
		if(i+k<seq.size()){
			updateMinimizer(kmerS, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
			kmer=min(kmerRC,kmerS);
		}else{
			return sketch;
		}
		++i;
	}while(true);
	return sketch;
}


//return a set containing all kmers
unordered_set<minimizer> allKmerSet(uint k,const string& seq){
	unordered_set<minimizer> set;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint i(0);
	do{
		set.insert(kmer);
		if(i+k<seq.size()){
			updateMinimizer(kmerS, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
			kmer=min(kmerRC,kmerS);
		}else{
			return set;
		}
		++i;
	}while(true);
	return set;
}


//compute the sketch of SEQ with H minimizers of size K, the read is separated in PART parts with H/PART minimizers each
vector<minimizer> minHashPart(uint H, uint k,const string& seq, uint part){
	vector<minimizer> result;
	uint size(seq.size()/part);
	for(uint i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}


void minHash2(uint H, uint k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	//	hashValue=hash(kmer);
	hashValue=hash64(kmer);
	for(uint j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}
	for(uint i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		// minimizer kmerS(seq2intStranded(seq.substr(i,k)));
		// minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(i,k))));
		kmer=(min(kmerS,kmerRC));
		hashValue=hash64(kmer);
		//		hashValue=hash(kmer);
		for(uint j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=hash64(hashValue);
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}


//compute H minimizers of size k from seq, minhash vanilla
vector<minimizer> minHash(uint H, uint k, const string& seq){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	// cout<<kmer<<endl;
	//hashValue=hash(kmer);
	//here I use a bad hash function for the first hash computation, this COULD lead to bas results,we could put a state of the art hash function as Murmurhash3
	hashValue=hash64(kmer);
	// int2seq(kmer, 16);
	// cin.get();

	for(uint j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}
	for(uint i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		// int2seq(kmerS, 16);
		// int2seq(kmerRC, 16);
		// minimizer kmerS(seq2intStranded(seq.substr(i,k)));
		// minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(i,k))));
		kmer=(min(kmerS,kmerRC));
		// int2seq(kmer, 16);
		// cin.get();
		hashValue=hash64(kmer);
		//hashValue=hash(kmer);
		for(uint j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=hash64(hashValue);
		}
	}
	return sketchs;
}


void minHash3(uint H, uint k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	hashValue=hash64(kmer);
	//~ hashValue=hash(kmer);
	for(uint j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}

	for(uint i(1);i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
			hashValue=hash64(kmer);
			//~ hashValue=hash(kmer);
			for(uint j(0);j<H;++j){
				if(hashValue<sketch[j]){
					sketch[j]=hashValue;
					sketchs[j]=kmer;
				}
				hashValue=hash64(hashValue);
			}
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}


//Compute a sketch of SEQ, with H minimizer of size K , but only kmers in the set FILTER are considered
vector<minimizer> minHashGenomic(uint H, uint k, const string& seq, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	uint i(1);
	while(filter.unordered_set::count(kmer)==0){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);
		++i;
	}
	hashValue=hash64(kmer);
	for(uint j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}

	for(;i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
			hashValue=hash64(kmer);
			//~ hashValue=hash(kmer);
			for(uint j(0);j<H;++j){
				if(hashValue<sketch[j]){
					sketch[j]=hashValue;
					sketchs[j]=kmer;
				}
				hashValue=hash64(hashValue);
			}
		}
	}
	return sketchs;
}


//Index all kmers (in a stranded way) to allow 1 error, nuc is the siez of the 'seed'
unordered_multimap<string,string> allKmerMapStranded(uint k,const string& seq, uint nuc){
	unordered_multimap<string,string> sketch;
	for(size_t i(0); i+k<=seq.size(); ++i){
		string kmer(seq.substr(i,k));
		sketch.insert({kmer.substr(0,nuc),kmer.substr(nuc)});
	}
	return sketch;
}


//stranded function !!!
unordered_multimap<string,string> minHashErrors(uint H, uint k, const string& seq, uint nuc){
	unordered_multimap<string,string> map;
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);

	minimizer kmer(seq2intStranded(seq.substr(0,k)));
	// hashValue=hash(kmer);
	//here I use a bad hash function for the first hash computation, this COULD lead to bas results,we could put a state of the art hash function as Murmurhash3
	uint64_t hashValue=hash64(kmer);
	for(uint j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}
	for(uint i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmer, seq[i+k], k);
		hashValue=hash64(kmer);
		//hashValue=hash(kmer);
		for(uint j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=j;
			}
			hashValue=hash64(hashValue);
		}
	}

	for(uint i(0); i<H; ++i){
		string kmer(seq.substr(sketchs[i],k));
		map.insert({kmer.substr(0,nuc),kmer.substr(nuc)});
	}

	return map;
}
