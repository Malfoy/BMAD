#include "minhash.h"
#include "xor.h"
#include "utils.h"


using namespace std;


void updateMinimizer(minimizer&	min, char nuc, uint64_t k){
	minimizer offset(1);
	offset<<=(2*k);
	min<<=2;
	min+=nuc2int(nuc);
	min%=offset;
}


void updateMinimizer32(uint64_t &	min, char nuc, uint64_t k){
	minimizer offset(1);
	offset<<=(2*k);
	min<<=2;
	min+=nuc2int(nuc);
	min%=offset;
}


void updateMinimizerEnd(minimizer&	min, char nuc, uint64_t k){
	min>>=2;
	min+=(nuc2int(nuc)<<(2*k-2));
}


void updateMinimizerRC(minimizer&	min, char nuc, uint64_t k){
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
}


void updateMinimizerRC32(uint64_t &	min, char nuc, uint64_t k){
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
}


//return a sketch containing all kmers
vector<minimizer> allKmer(uint64_t k,const string& seq){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint64_t i(0);
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


//return a sketch containing all genomic kmers
vector<minimizer> allGenomicKmers(uint64_t k,const string& seq,const unordered_set <minimizer>& set){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint64_t i(0);
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
unordered_set<minimizer> allKmerSet(uint64_t k,const string& seq){
	unordered_set<minimizer> set;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	uint64_t i(0);
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
vector<minimizer> minHashPart(uint64_t H, uint64_t k,const string& seq, uint64_t part){
	vector<minimizer> result;
	uint64_t size(seq.size()/part);
	for(uint64_t i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}


void minHash2(uint64_t H, uint64_t k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	//	hashValue=hash(kmer);
	hashValue=hash64(kmer);
	for(uint64_t j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}
	for(uint64_t i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		// minimizer kmerS(seq2intStranded(seq.substr(i,k)));
		// minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(i,k))));
		kmer=(min(kmerS,kmerRC));
		hashValue=hash64(kmer);
		//		hashValue=hash(kmer);
		for(uint64_t j(0); j<H; ++j){
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
vector<minimizer> minHash(uint64_t H, uint64_t k, const string& seq){
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

	for(uint64_t j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}
	for(uint64_t i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		// minimizer kmerS(seq2intStranded(seq.substr(i,k)));
		// minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(i,k))));
		kmer=(min(kmerS,kmerRC));
		// int2seq(kmer, 16);
		// cin.get();
		hashValue=hash64(kmer);
		//hashValue=hash(kmer);
		for(uint64_t j(0); j<H; ++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=hash64(hashValue);
		}
	}
	return sketchs;
}


void minHash3(uint64_t H, uint64_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	hashValue=hash64(kmer);
	for(uint64_t j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=hash64(hashValue);
	}

	for(uint64_t i(1);i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
			hashValue=hash64(kmer);
			for(uint64_t j(0);j<H;++j){
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
vector<minimizer> minHashGenomic(uint64_t H, uint64_t k, const string& seq, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	uint64_t i(1);
	while(filter.unordered_set::count(kmer)==0 and i+k<seq.size()){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);
		++i;
	}
	hashValue=hash64(kmer);
	for(uint64_t j(0);j<H;++j){
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
			for(uint64_t j(0);j<H;++j){
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


bool isCorrect(minimizer seq, minimizer ref, uint64_t n){
	// return true;
	for(uint64_t i(1); i<n; ++i){
		unsigned char s(seq>>(2*(n-i))),r(ref>>(2*(n-i)));
		minimizer offset2 (1<<(2*(n-i)));
		seq%=(offset2);
		ref%=(offset2);
		if((s)!=r){
			minimizer offset (1<<(2*(n-i-1)));
			if((seq>>(2*(n-i-1)))==r){
				return (seq%(offset)==(ref>>2));
			}
			if(s==ref>>(2*(n-i-1))){
				return (ref%(offset)==(seq>>2));
			}
			return seq==ref;
		}
	}
	return true;
}

//TODO optimize with array
unordered_multimap<minimizer,minimizer> allKmerMap(const uint64_t k,const string& seq, const  uint64_t nuc){
	unordered_multimap<minimizer,minimizer> map;
	minimizer kmer (seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC (rc(kmer,k));
	minimizer seed(getBegin(kmer, k-nuc));
	minimizer body(getEnd(kmer, k-nuc));
	minimizer seedRC(getBegin(kmerRC, k-nuc));
	minimizer bodyRC(getEnd(kmerRC, k-nuc));
	// int2seq(kmer, k);
	// int2seq(seed, nuc);
	// int2seq(body, k-nuc);
	// int2seq(kmerRC, k);
	// int2seq(seedRC, nuc);
	// int2seq(bodyRC, k-nuc);
	// cin.get();
	for(uint64_t i(0); ; ++i){
		map.insert({seed,body});
		map.insert({seedRC,bodyRC});
		if(i+k<seq.size()){
			updateMinimizer(kmer,seq[i+k],k);
			updateMinimizerRC(kmerRC,seq[i+k],k);
			seed=(getBegin(kmer, k-nuc));
			body=(getEnd(kmer, k-nuc));
			seedRC=(getBegin(kmerRC, k-nuc));
			bodyRC=(getEnd(kmerRC, k-nuc));
			// updateMinimizer(seed,seq[i+nuc],nuc);
			// updateMinimizer(body,seq[i+k],k-nuc);
			// updateMinimizerRC(seedRC,seq[i+nuc],nuc);
			// updateMinimizerRC(bodyRC,seq[i+k],k-nuc);
			// int2seq(kmer, k);
			// int2seq(seed, nuc);
			// int2seq(body, k-nuc);
			// int2seq(kmerRC, k);
			// int2seq(seedRC, nuc);
			// int2seq(bodyRC, k-nuc);
			// cin.get();
		}else{
			return map;
		}
	}
	return map;
}


bool compareMinimizer (minimizer i,minimizer j) { return (hash64(i)<hash64(j)); }



//return a sketch containing all quasi genomic kmers
vector<minimizer> allQuasiGenomicKmers(uint64_t k,const string& seq,const unordered_multimap<minimizer,minimizer>& map,uint64_t nuc){
	// cout<<"begin"<<endl;
	vector<minimizer> sketch;
	vector<minimizer> tmp;
	minimizer seed(seq2intStranded(seq.substr(0,nuc)));
	minimizer body(seq2intStranded(seq.substr(nuc,k-nuc)));
	// int2seq(seed,nuc);
	// int2seq(body,k-nuc);
	// int2seq(seq2intStranded(seq.substr(0,k)),k);
	uint64_t i(0);
	do{
		// cout<<"loopbeg"<<endl;
		tmp={};
		auto range = map.equal_range(seed);
    	for (auto it = range.first; it != range.second; ++it){
			if(isCorrect(body,it->second,k-nuc)){
				tmp.push_back(getRepresent(cat(seed,it->second,k-nuc),k));
				
				// cout<<"accepted"<<endl;
				// int2seq(cat(seed,it->second,nuc), k);
				// int2seq(cat(seed,body,nuc), k);
				// int2seq(body,k-nuc);
				// int2seq(it->second,k-nuc);
			}
		}
		if(!tmp.empty()){
			sort(tmp.begin(),tmp.end());
			sketch.push_back(tmp[0]);
		}
		// cout<<"update"<<endl;
		if(i+k<seq.size()){
			// int2seq(seq2intStranded(seq.substr(i+1,k)),k);
			updateMinimizer(seed,seq[i+nuc],nuc);
			updateMinimizer(body,seq[i+k],k-nuc);
			// int2seq(seed,nuc);
			// int2seq(body,k-nuc);
			// int2seq(cat(seed,body,k-nuc), k);
			// cout<<"next"<<endl;
			// cin.get();
		}else{
			// cout<<"end"<<endl;
			// cout<<sketch.size();
			// cin.get();
			return sketch;
		}
		++i;
	}while(true);
	return sketch;
}


vector<minimizer> smallestGenomicKmers(uint64_t H,uint64_t k, const string& seq,unordered_set <minimizer> set){
	vector<minimizer> all(allGenomicKmers(k,seq,set));
	sort(all.begin(),all.end(),compareMinimizer);
	vector<minimizer> res(H);
	for(uint64_t i(0);i<H;++i){
		res[i]=all[i];
	}
	return res;
}


unordered_map<minimizer,vector<readNumber>> indexReadSet(const string& readFile, const uint64_t k, const uint64_t seedSize,unordered_multimap<minimizer,minimizer> map){
	ifstream readS(readFile,ios::in);
	vector<string> reads(getReads(readS, 1000));
	unordered_map<minimizer,vector<readNumber>> min2reads;
	vector<minimizer> sketch;
	readNumber rn(0);
	for(uint64_t i(0);i<reads.size();++i,++rn){
		sketch=allQuasiGenomicKmers(k,reads[i],map,seedSize);
		for(uint64_t ii(0);ii<sketch.size();++ii){
			(min2reads[sketch[ii]]).push_back(rn);
		}
	}
	return min2reads;
}


unordered_map<minimizer,uint64_t > kmerCounting(const string& readFile, const uint64_t k){
	unordered_map<minimizer,uint64_t > count;
	string seq;
	ifstream readS(readFile);
	vector<string> reads;
	while(!readS.eof()){
		reads=(getReads(readS, 1000));
		// cout<<reads.size()<<endl;
		for(uint64_t ii(0);ii<reads.size();++ii){
			seq=reads[ii];
			uint64_t i(0);
			minimizer kmerS(seq2intStranded((seq.substr(0,k))));
			minimizer kmerRC(rc(kmerS,k));
			minimizer kmer(min(kmerRC,kmerS));
			bool end(false);
			do{
				count[kmer]=min(count[kmer]+1,(uint64_t)200);
				if(i+k<seq.size()){
					updateMinimizer(kmerS, seq[i+k], k);
					updateMinimizerRC(kmerRC, seq[i+k], k);
					kmer=min(kmerRC,kmerS);
					++i;
				}else{
					end=true;
				}
			}while(!end);
		}
	}
	return count;
}


unordered_set<minimizer> getSolidSet(unordered_map<minimizer,uint64_t >& count, uint64_t T){
	unordered_set<minimizer> solid;
	for ( auto it = count.begin(); it != count.end(); ++it ){
		// cout<<it->second<<endl;
		if(it->second>=T){
			// cout<<"lol"<<endl;
			solid.insert(it->first);
		}
	}
	return solid;
}


unordered_multimap<minimizer,minimizer> getSolidMap(unordered_map<minimizer,uint64_t >& count, uint64_t T, uint64_t k, uint64_t nuc){
	unordered_multimap<minimizer,minimizer> map;
	for ( auto it = count.begin(); it != count.end(); ++it ){
		if(it->second>=T){
			minimizer kmer(it->first);
			minimizer seed(getBegin(kmer, k-nuc));
			minimizer body(getEnd(kmer, k-nuc));
			minimizer kmerRC (rc(kmer,k));
			minimizer seedRC(getBegin(kmerRC, k-nuc));
			minimizer bodyRC(getEnd(kmerRC, k-nuc));
			map.insert({seed,body});
			map.insert({seedRC,bodyRC});
		}
	}
	// count={};
	return map;
}


void bench(){
	uint64_t H(1000);
	uint64_t nuc(10);
	uint64_t k(11),part(2);
	ifstream readFile1("reads_virus_10k_0001.ref");
	ifstream readFile2("reads.fa");
	vector<string> V1,Vref;
	Vref=getReads(readFile1, 1);
	string ref(Vref[0]);
	cout<<"ref"<<ref.size()<<endl;
	V1=getReads(readFile2,100);
	for(uint64_t j(0);j<15;j+=2){
	uint64_t i(0),res1(0),res2(0),res3(0),res4(0),res5(0),res6(0),res7(0),res8(0),res9(0),res10(0);
		k=11+j;
		cout<<"H:"<<H<<endl;
		cout<<"k:"<<k<<endl;
		for(;i+1<V1.size();i+=1){
			string seq1(V1[i]),seq2(V1[i+1]);
			// string seq1(V1[i]),seq2(randomSeq(seq1.size()));
			// cout<<seq1.size()<<" "<<seq2.size()<<endl;
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
		<<"minhash vanilla                  : "<<res1/i<<endl
		<<"minhash with minimizer reparition: "<<res3/i<<endl
		<<"minhash with only genomic kmers  : "<<res5/i<<endl
		<<"genomic kmers                    : "<<res7/i<<endl
		<<"quasi genomic kmer               : "<<res9/i<<endl<<endl;
	}
}
