#include "minhash.h"

using namespace std;


string randomSeq(uint32_t length){
	auto randchar=[]() -> char{
		const char charset[] ="ATCG";
		const uint32_t max_index = (sizeof(charset) - 1);
		return charset[ rand() % max_index ];
	};
	string str(length,0);
	generate_n( str.begin(), length, randchar );
	return str;
}


uint8_t nuc2int(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return 0;
}


char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
/*
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
*/
	return 'X';//error
}


string reversecomplement (const string& s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc += revcomp(s[i]);
	}
	return rc;
}


void updateMinimizer(minimizer&	min, char nuc, uint8_t k){
	minimizer offset(1<<(2*k));
	min<<=2;
	min+=nuc2int(nuc);
	min%=offset;
}


void updateMinimizerEnd(minimizer&	min, char nuc, uint8_t k){
	min>>=2;
	min+=(nuc2int(nuc)<<(2*k-2));
}


void updateMinimizerRC(minimizer&	min, char nuc, uint8_t k){
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
}



//simplest (fastest ?) hash function
uint64_t xorshift64(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}


//A sequence and its reverse complement are represented by the same sequence
string getRepresent (const string& str){
	return (min(str,reversecomplement(str)));
}


string getRepresent2(const string& s){
	for (int i = 0; i < (int)s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]);
		if (s[i] < c) {
			return s;
		} else if (s[i] > c) {
			return reversecomplement(s);
		}
	}
	return s;
}


//return a sketch containing all kmers
vector<minimizer> allHash(uint8_t k,const string& seq){
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


unordered_set <minimizer> allKmerSet(uint8_t k,const string& seq){
	unordered_set<minimizer> sketch;
	for(uint i(0);i+k<=seq.size();++i){
		sketch.insert(seq2int(seq.substr(i,k)));
	}
	return sketch;
}


minimizer seq2int(const string& seq){
	string str(getRepresent(seq));
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(str[i]);
	}
	return res;
}


//compute the sketch of SEQ with H minimizers of size K, the read is separated in PART parts with H/PART minimizers each
vector<minimizer> minHashpart(uint32_t H, uint8_t k,const string& seq, uint8_t part){
	vector<minimizer> result;
	uint size(seq.size()/part);
	for(uint i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}


void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	//	hash<uint32_t> hash;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	//	hashValue=hash(kmer);
	hashValue=xorshift64(kmer);
	for(uint j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}
	for(uint i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		kmer=(min(kmerS,kmerRC));
		hashValue=xorshift64(kmer);
		//		hashValue=hash(kmer);
		for(uint j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift64(hashValue);
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}

//compute H minimizers of size k from seq
vector<minimizer> minHash(uint32_t H, uint8_t k, const string& seq){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	//	hash<uint32_t> hash;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	//hashValue=hash(kmer);
	//here I use a bad hash function for the first hash computation, this COULD lead to bas results,we could put a state of the art hash function as Murmurhash3
	hashValue=xorshift64(kmer);
	for(uint j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}
	for(uint i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		kmer=(min(kmerS,kmerRC));
		hashValue=xorshift64(kmer);
		//hashValue=hash(kmer);
		for(uint j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift64(hashValue);
		}
	}
	return sketchs;
}


void minHash3(uint32_t H, uint8_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;;
	//~ hash<uint32_t> hash;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	hashValue=xorshift64(kmer);
	//~ hashValue=hash(kmer);
	for(uint j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(uint i(1);i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
			hashValue=xorshift64(kmer);
			//~ hashValue=hash(kmer);
			for(uint j(0);j<H;++j){
				if(hashValue<sketch[j]){
					sketch[j]=hashValue;
					sketchs[j]=kmer;
				}
				hashValue=xorshift64(hashValue);
			}
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}


//Compute a sketch of SEQ, with H minimizer of size K distibuted on PART parts, but only kmers in the set FILTER are considered
vector<minimizer> minHashpart2(uint32_t H, uint8_t k, const string& seq, uint8_t part, const unordered_set<minimizer>& filter){
	vector<minimizer> result;
	uint size(seq.size()/part);
	for(uint i(0);i<part;++i){
		//~ minHash3(H/part,k,seq.substr(i*size,size+16),result,filter);
		minHash3(H/part,k,seq.substr(i*size,size),result,filter);
	}
	return result;
}


minimizer seq2intStranded(const string& seq){
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(seq[i]);
	}
	return res;
}


unordered_multimap<string,string> allKmerMapStranded(size_t k,const string& seq, char nuc){
	unordered_multimap<string,string> sketch;
	for(size_t i(0); i+k<=seq.size(); ++i){
		string kmer(seq.substr(i,k));
		sketch.insert({kmer.substr(0,nuc),kmer.substr(nuc)});
	}
	return sketch;
}


bool equalStr(const string& seq1, const string& seq2){
	uint size(min(seq1.size(),seq2.size()));
	return (seq1.substr(0,size))==seq2.substr(0,size);
}


//rewrite with ternary operator
bool isCorrect(const string& seq,const string& ref){
	for(uint i(0); i<seq.size(); ++i){
		if(seq[i]!=ref[i]){
			if(seq[i+1]==ref[i]){
				return equalStr(seq.substr(i+2),ref.substr(i+1));
			}
			if(seq[i]==ref[i+1]){
				return equalStr(seq.substr(i+1),ref.substr(i+2));
			}
			return (seq.substr(i+1)==ref.substr(i+1));
		}
	}
	return true;
}


double percentStrandedErrors(uint8_t k, const string& seq, const unordered_multimap<string, string>& genomicKmers, char nuc){
	double inter(0);
	string kmer;
	kmer.reserve(k);
	uint i(0);
	for(; i+k<=seq.size(); ++i){
		kmer=seq.substr(i,k);
		if(kmer.size()!=k){
			cout<<"wtf"<<endl;
		}
		auto range(genomicKmers.equal_range(kmer.substr(0,nuc)));
		for (auto it(range.first); it!=range.second; it++){
			if(isCorrect(kmer.substr(nuc),it->second)){
				inter++;
				break;
			}else{
			}
		}
	}
	return double(100*inter/(seq.size()-k+1));;
}


uint32_t sketchHammingComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2){
	uint32_t res(0);
	for(uint i(0); i<sketch1.size(); ++i){
		if(sketch1[i]==sketch2[i]){++res;}
	}
	return res;
}


uint32_t sketchComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2){
	uint32_t res(0);
	unordered_set<minimizer> minimizerSet;
	for(uint i(0); i<sketch1.size(); ++i){
		minimizerSet.insert(sketch1[i]);
	}
	for(uint i(0); i<sketch2.size(); ++i){
		if(minimizerSet.count(sketch2[i])!=0){++res;}
	}
	return res;
}
