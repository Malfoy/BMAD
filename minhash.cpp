#include "minhash.h"

using namespace std;


char nuc2int(char c){
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


string reversecomplement (const string &s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc += revcomp(s[i]);
	}
	return rc;
}


uint64_t xorshift64(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}

string getRepresent (const string& str){
	string rc(reversecomplement(str));
	if(rc<str){
		return rc;
	}else{
		return str;
	}
}


string getRepresent2(const string &s){
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


vector<minimizer> allHash(size_t k,const string& seq){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	size_t i(0);
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


unordered_set <minimizer> allKmerSet(size_t k,const string& seq){
	unordered_set<minimizer> sketch;
	for(size_t i(0);i+k<=seq.size();++i){
		sketch.insert(seq2int(seq.substr(i,k)));
	}
	return sketch;
}

minimizer seq2int(const string& seq){
	string str(getRepresent(seq));
	//	cout<<"lol"<<endl;
	//	cin.get();
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(str[i]);
	}
	return res;
}

vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
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
	for(size_t j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}
	for(size_t i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		kmer=(min(kmerS,kmerRC));
		hashValue=xorshift64(kmer);
		//		hashValue=hash(kmer);
		for(size_t j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift64(hashValue);
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}

void updateMinimizer(minimizer&	min, char nuc,size_t k){
	minimizer offset(1<<(2*k));
	min<<=2;
	min+=nuc2int(nuc);
	min%=offset;
}


void updateMinimizerEnd(minimizer&	min, char nuc,size_t k){
	min>>=2;
	min+=(nuc2int(nuc)<<(2*k-2));
}


void updateMinimizerRC(minimizer&	min, char nuc,size_t k){
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
}


void minHash3(size_t H, size_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;;
	//~ hash<uint32_t> hash;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	hashValue=xorshift64(kmer);
	//~ hashValue=hash(kmer);
	for(size_t j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(size_t i(1);i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
			hashValue=xorshift64(kmer);
			//~ hashValue=hash(kmer);
			for(size_t j(0);j<H;++j){
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

vector<minimizer> minHashpart2(size_t H, size_t k,const string& seq, size_t part, const unordered_set<minimizer>& filter){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
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
