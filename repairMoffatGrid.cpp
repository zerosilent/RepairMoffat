#include <iostream>
#include <vector>
#include <unordered_map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <math.h>
#include "Pair.cpp"
#include "PairFrequency.cpp"
#include "Record.cpp"
#include "PriorityQueue.cpp"
#include "intAlphabet.cpp"
#include "wmPairs.cpp"

using namespace std;
using namespace sdsl;

class Repair{
	public:
		vector<Record> *secuenceSymbols;
		unordered_map<Pair,PairFrequency,PairHasher> *hashTable;
		unordered_map<int,Pair> *dictionary;
		unordered_map<int,int> *symbolOcurrencyHash;
		unordered_map<int,int> *firstSymbolOcurrencyHash;
		wt_gmr<> sizeTextNonTerminalSymbols;
		wt_gmr<> perm;//Permutation for large alphabet
		
		wm_int<> wmatrix;//Grid to store elements
		int_vector<> indexPositionsNonTerminals;//Added by jconcha: keeps index of A[i] of points in wmatrix
		int_vector<> indexPositionsColumn;////Added by jconcha: keeps index of columns of points in wmatrix
		
		sd_vector<> sdb_wmatrix;//Sparse bitvector for restore original coordinates on wmatrix
		int biggerTerminalSymbolNumber;
		int moreFrequentPair;
		int nonTerminalSymbolNumber;
		PriorityQueue *pQueue;

		Repair(){
			secuenceSymbols = new vector<Record>();
			hashTable = new unordered_map<Pair,PairFrequency,PairHasher>;
			hashTable->max_load_factor(0.7);
			dictionary = new unordered_map<int,Pair>;
			symbolOcurrencyHash = new unordered_map<int,int>;//Hash to store symbol and last occurency
			firstSymbolOcurrencyHash = new unordered_map<int,int>;//Hash to store symbol and first occurency
		}

		~Repair(){
			delete secuenceSymbols;
			delete hashTable;
			delete pQueue;
		}
		
		void deleteMemory(){
			delete secuenceSymbols;
			delete hashTable;
			delete pQueue;
		}

		void fillPriorityQueue(){
			unordered_map<Pair,PairFrequency,PairHasher>::iterator got = hashTable->begin();
			int n = secuenceSymbols->size();
			pQueue = new PriorityQueue(n);

			pQueue->freqMoreFrequentPair = moreFrequentPair;
			while(got!= hashTable->end()){
				PairFrequency pf = got->second;
				Pair p = got->first;

				PairWithFrequencies *pwf = new PairWithFrequencies(p.current_symbol_Num,p.next_symbol_Num,pf.first_pair_Pos,pf.frequency);

				if(pf.frequency < ceil(sqrt(n)) - 1){
					//cout << "pwf:"<<pwf->frequency<<endl;
					(*(pQueue->lessFrequentPairs))[pf.frequency].push_back(*pwf);
				}
				else{
					//cout << "pwf:"<<pwf->frequency<<endl;
					pQueue->moreFrequentPairs->push_back(*pwf);
				}
				got++;
			}
		}


		void fillSecuenceSymbolsAndHashTable(vector<int> &symbols){
			int k = 0;
			for(unsigned int i = 0; i < symbols.size() - 1; i++){
				Pair p;

				p.current_symbol_Num = symbols[i];
				p.next_symbol_Num = symbols[i+1];

				Record r;
				r.symbolNumber = symbols[i];

				unordered_map<Pair,PairFrequency,PairHasher>::const_iterator got = hashTable->find(p);
				unordered_map<int,int>::const_iterator gotSymbol = symbolOcurrencyHash->find(symbols[i]);

				if (got == hashTable->end()){
					PairFrequency pf;
					pf.first_pair_Pos = i;
					pf.last_pair_Pos = i;
					pf.frequency = 1;

					if(moreFrequentPair < pf.frequency)
						moreFrequentPair = pf.frequency;

					(*hashTable)[p] = pf;


					r.prevNeighborPair = -1;
					r.nextNeighborPair = -1;
					r.prevSymbolPos = -1;
					r.nextSymbolPos = -1;
					r.prevUsablePos = -1;
					r.nextUsablePos = -1;



					if(gotSymbol != symbolOcurrencyHash->end()){
						r.prevSymbolPos = gotSymbol->second;
						(*secuenceSymbols)[gotSymbol->second].nextSymbolPos = i;
					}
					else{
						(*firstSymbolOcurrencyHash)[symbols[i]] = i;
					}
					(*symbolOcurrencyHash)[symbols[i]] = i;
				}
				else{
					PairFrequency pf = got->second;
					pf.frequency = pf.frequency + 1;
					if(moreFrequentPair < pf.frequency)
						moreFrequentPair = pf.frequency;

					//Update previous Symbol pair
					(*secuenceSymbols)[pf.last_pair_Pos].nextNeighborPair = i;
					r.prevNeighborPair = pf.last_pair_Pos;
					r.nextNeighborPair = -1;
					pf.last_pair_Pos = i;
					
					(*hashTable)[p] = pf;
					//cout << "p:" << p.current_symbol_Num << "," << p.next_symbol_Num << endl;

					//Update previous symbol
					if(gotSymbol != symbolOcurrencyHash->end()){
						r.prevSymbolPos = gotSymbol->second;
						r.nextSymbolPos = -1;
						(*secuenceSymbols)[gotSymbol->second].nextSymbolPos = i;
					}
					(*symbolOcurrencyHash)[symbols[i]] = i;
				}
				r.prevUsablePos = k - 1;
				r.nextUsablePos = k + 1;
				k++;
				secuenceSymbols->push_back(r);
			}
			//Last symbol
			Record r;
			int i = symbols.size() - 1;
			r.symbolNumber = symbols[i];
			r.prevNeighborPair = -1;
			r.nextNeighborPair = -1;
			r.prevSymbolPos = -1;
			r.nextSymbolPos = -1;
			unordered_map<int,int>::const_iterator gotSymbol = symbolOcurrencyHash->find(symbols[i]);
			//Update previous symbol
			if(gotSymbol != symbolOcurrencyHash->end()){
				r.prevSymbolPos = gotSymbol->second;
				r.nextSymbolPos = -1;
				(*secuenceSymbols)[gotSymbol->second].nextSymbolPos = i;
			}
			(*symbolOcurrencyHash)[symbols[i]] = i;
			r.prevUsablePos = k - 1;
			r.nextUsablePos = k + 1;
			k++;
			secuenceSymbols->push_back(r);
			delete(symbolOcurrencyHash);
		}

		void printHashTable(){
			for(unordered_map<Pair,PairFrequency,PairHasher>::iterator it = hashTable->begin(); it != hashTable->end(); ++it) {
			  cout << "H:" << it->first.current_symbol_Num <<","<< it->first.next_symbol_Num << " f:"<< it->second.frequency << endl;
			}
		}

		void printDictionary(){
			cout<<"dictionary:"<<endl;
			for(unordered_map<int,Pair>::iterator it = dictionary->begin(); it != dictionary->end(); ++it) {
			  cout << "n:"<< it->first << " pair:"<< it->second.current_symbol_Num <<","<< it->second.next_symbol_Num << endl;
			}
		}
		
		void printDictionary_perm(){
			cout<<"dictionary:"<<endl;
			
			int k = 0;
			int element = 0;
			//int element = biggerTerminalSymbolNumber + 1;
			while(k + 1 < perm.size()){
				//cout << "n:" << wmatrix.pointWeights[element] << " ";
				cout << "n:" << element + biggerTerminalSymbolNumber + 1 << ",";
				//k+1 before k, because perm stores C[i]B[i] instead of B[i]C[i]
				cout <<"pair:" << perm[k+1] << "," << perm[k] << endl;
				k+=2;
				element++;
			}
		}

		void printAllSecuenceSymbols(){
			for(unsigned int i = 0; i < secuenceSymbols->size();i++){
				int symbol = (*secuenceSymbols)[i].symbolNumber;
				cout << symbol << " ";
			}
			cout << endl;
			cout <<"--Positions--"<<endl;
			for(unsigned int i = 0; i < secuenceSymbols->size();i++){
				int symbol = (*secuenceSymbols)[i].nextUsablePos;
				cout << symbol << " ";
			}
			cout << endl;
		}
		/*
		 * Repair pairs incrementally taking most frequent pair, but we need to rename it to "order it" by alphabetically order
		 * This renaming algorithm takes O(n), using two hashes to rename dictionary and pairs non terminal symbols
		 * */
		void renameNonTerminalsPairsAndDictionary(WmPairs *pairs){
			bit_vector bv_present(nonTerminalSymbolNumber);
			unordered_map<int,int> hashReplacement;
			//unordered_map<int,int> reversedHashReplacement;
			
			for(int i = 0; i < pairs->size();i++){
				bv_present[pairs->v_pair[i].column] = 1;
			}
			//cout << "bitvector present:"<<bv_present<<endl;
			
			//cout << "ones positions:" << endl;
			
			size_t ones = bit_vector::rank_1_type(&bv_present)(bv_present.size());
			
			//TODO:Special care must be taken to not replace with a number lower than nonTerminals.
			bit_vector::select_1_type sb_sel(&bv_present);
			/*
			for(int i = 1; i <= ones;i++){
				cout << sb_sel(i) << endl;
			}
			*/
			//cout << endl;
			
			//cout << "renaming:" << endl;
			
			int i = 1;
			
			int j = 0;
			
			int current;
			int previous = pairs->v_pair[j].column;
			int currentReplacement = sb_sel(i);
			do{
				current = pairs->v_pair[j].column;
				//Increments i to make selection faster :-)
				if(current != previous){
					i++;
					if(i <= ones){
						currentReplacement = sb_sel(i);//sel in this case is log log u, but is executed only r times (r log log u)
					}
				}
				hashReplacement[pairs->v_pair[j].column] = currentReplacement;
				//cout << pairs->v_pair[j].column <<","<<currentReplacement<< endl;
				previous = current;
				j++;
			}while(j < pairs->size());
			
			//Replacements in pairs:
			for(int i = 0; i < pairs->size();i++){
				/*
				pairs->v_pair[i].row = 
				bv_present[pairs->v_pair[i].column] = 1;
				*/
				unordered_map<int,int>::iterator gotRow = hashReplacement.find(pairs->v_pair[i].row);
				unordered_map<int,int>::iterator gotColumn = hashReplacement.find(pairs->v_pair[i].column);
				//unordered_map<int,int>::iterator gotPw = hashReplacement.find(wmatrix.pointWeights[i]);
				
				if(gotRow != hashReplacement.end()){
					pairs->v_pair[i].row = gotRow->second;
				}
				if(gotColumn != hashReplacement.end()){
					pairs->v_pair[i].column = gotColumn->second;
				}
				/*if(gotPw != hashReplacement.end()){
					wmatrix.pointWeights[i] = gotPw->second;
				}*/
			}
			
			//TODO: Renaming also dictionary
			//Represent permutation using Golynski and Rao permutations for large alphabets
			convertRulesDictToIntVector(hashReplacement);
			//delete dictionary;
			//pairs->print();
		}
		
		void fillWMatrixGrid(){
			WmPairs *pairs = new WmPairs();
						
			int k = 0;
			for(int i = biggerTerminalSymbolNumber + 1; i < nonTerminalSymbolNumber;i++){
				pairs->insert((*dictionary)[i].next_symbol_Num,(*dictionary)[i].current_symbol_Num,i);
			}			
			cout << "---ORIGINAL---" << endl;
			pairs->print();
			
			//wmatrix.fillPointWeights(pairs);
			pairs->fillindexPositions();
			
			
			//sort_columnOnlyByNumber
			
			pairs->sort_columnOnlyByNumber(dictionary,biggerTerminalSymbolNumber);
			cout << "---SORTED BY COLUMN---" << endl;
			pairs->print();	
			
			
			pairs->sort(dictionary,biggerTerminalSymbolNumber);
			//renameNonTerminalsPairsAndDictionary(pairs);
			//cout << "---"<<endl;
			cout << "---SORTED---" << endl;
			pairs->print();	
			
			convertDictionaryToPerm();
			/*
			renameNonTerminalsPairsAndDictionary(pairs);
			cout << "---RENAMED---" << endl;
			pairs->print();
			cout << "lala" << endl;
			*/
			
			int numberPoints = pairs->v_pair.size();
			int numRows = pairs->getNumberRows();

			//cout << "pairs size:"<<pairs->size()<<endl;
			bit_vector bv_wmatrix = bit_vector(2*numRows+numberPoints);//Doesn't matter the extra space for the moment, it will be resized.
			
			pairs->fillBitvectorB(bv_wmatrix);
			//cout << "bv_wmatrix:"<<bv_wmatrix<<endl;
			//pairs->printPosZeroesBitVectorB(bv_wmatrix);
			sdb_wmatrix = bv_wmatrix; 

			pairs->renumberPairs();
			cout << "---RENUMBERED PAIRS---" << endl;
			pairs->print();
			
			cout << "sdb_wmatrix:" << sdb_wmatrix << endl;
			
			int minColumn = INT_MAX;
			
			for(int i = 0; i < pairs->v_pair.size(); i++){
				if(minColumn > pairs->v_pair[i].column){
					minColumn = pairs->v_pair[i].column;
				}
				//cout <<"i:"<<i<< " pair:" << WmPairs::getOldCoordinate(i+1,v[i],bv_wmatrix).column << "," << WmPairs::getOldCoordinate(i+1,v[i],bv_wmatrix).row << endl;
				//cout << v[i] << " ";
			}		
			
			indexPositionsNonTerminals = pairs->indexPositionsNonTerminals;
			indexPositionsColumn = pairs->indexPositionsColumn;
			
			cout << "indexPositionsColumn:" << indexPositionsColumn << endl;
			cout << "indexPositionsNonTerm:" << indexPositionsNonTerminals << endl;
			
			for(int i = 0; i < pairs->v_pair.size(); i++){
				//cout << pairs->v_pair[i].row << ","<< pairs->v_pair[i].column << endl;
				//TODO: Add to getOldCoordinateSdb int_vector with C[i]B[i] values, instead of accesing column, position of column should be accessed 
				
				//cout << wmatrix.pointWeights[pairs->v_pair[i].column - 1] << ", col:" << pairs->v_pair[i].row << endl;
				//cout << "A[i]:" << wmatrix.indexPositionsNonTerminals[i] + biggerTerminalSymbolNumber + 1;
				
				
				//WmPair p = pairs->getOldCoordinateSdb(wmatrix.pointWeights[pairs->v_pair[i].column], pairs->v_pair[i].row,sdb_wmatrix);
				
				//cout << "," << pairs->indexPositionsColumn[i] << endl;
				WmPair p = WmPairs::getOldCoordinateSdb(pairs->indexPositionsColumn[i] + 1 , pairs->v_pair[i].row,sdb_wmatrix);
				cout << "o: "<< indexPositionsNonTerminals[i] + biggerTerminalSymbolNumber + 1 <<" pos:"<< p.row<<","<<p.column<<endl;
				
			}
			cout << endl;
			
			
			
			//cout << "pointweights:"<< wmatrix.pointWeights << endl;
			
			//int_vector<> v(pairs->v_pair.size());
		
			/*cout << "order nonTerm: "<<endl;
			for(int i = 0; i < pairs->v_pair.size(); i++){
				cout << pairs->v_pair[i].nonTerminalSymbol << " ";
			}
			cout << endl;
		*/
			//cout << "matrix asdf:"<<endl;
			/*for(int i = 0; i < pairs->v_pair.size(); i++){
				v[i] = pairs->v_pair[i].row;
				//v[i] = pairs->indexPositionsColumn[i];
				//cout <<"i:"<<i<< " pair:" << WmPairs::getOldCoordinate(i+1,v[i],bv_wmatrix).column << "," << WmPairs::getOldCoordinate(i+1,v[i],bv_wmatrix).row << endl;
				//cout << v[i] << " ";
			}*/
			delete(pairs);//Also with delete call destructor
			/*
			cout << "lala" << endl;
			cout << v[1] << endl;
			cout << WmPairs::getOldCoordinate(1,v[1],bv_wmatrix).row << "," << WmPairs::getOldCoordinate(1,v[1],bv_wmatrix).column << endl;
			//cout << endl;
			cout << "endllala" << endl;
			*/
			//wm_int<> wmatrix;
			construct_im(wmatrix, indexPositionsNonTerminals);
			cout << "wmatrix:" << wmatrix << endl;
			/*
			 * Careful, range_search_2d returns column,row instead of row,column
			 * */
			 
			WmPair p = WmPairs::getNewNumerationSdb(105,112,sdb_wmatrix);
			cout << "lalalalalala:" << p.row << endl;
			/*
			
			WmPair p = WmPairs::getOldCoordinateSdb(pairs->indexPositionsColumn[i] + 1 , pairs->v_pair[i].row,sdb_wmatrix);
			pair<long unsigned int, vector<pair<long unsigned int, long unsigned int> > > resultQuery = wmatrix.range_search_2d(110,115,110,115,true);
			
			
			
			cout << resultQuery.first << endl;
			cout << "minCol:"<<minColumn << endl;
			for(int i = 0; i < resultQuery.second.size();i++){
				//Mincolumn should be added in the rename	
				cout << resultQuery.second[i].second << "," <<minColumn + resultQuery.second[i].first<<  endl;
				WmPair p = pairs->getOldCoordinateSdb(minColumn + resultQuery.second[i].first - 1, resultQuery.second[i].second,sdb_wmatrix);
				cout << "o:"<< p.column<<","<<p.row<<endl;
			}
			*/

			//Free unused structures
			//v.~int_vector();//Destructor
			
		}
		
		void convertDictionaryToPerm(){
			int size = dictionary->size() * 2;
			
			int_vector<> v(size);
			int_vector<> len_dict_nonTermSym(dictionary->size());
			
			int k = 0;
			for(int i = biggerTerminalSymbolNumber + 1; i < nonTerminalSymbolNumber;i++){
				int firstNum = (*dictionary)[i].current_symbol_Num;
				int secondNum = (*dictionary)[i].next_symbol_Num;
				
				v[k] = firstNum;
				k++;
				v[k] = secondNum;
				k++;
				
				int sizeLeft = firstNum - (biggerTerminalSymbolNumber + 1);
				int sizeRight = secondNum - (biggerTerminalSymbolNumber + 1);
				
				int sizeTextNonTerminalSymbol; 
				len_dict_nonTermSym[i - (biggerTerminalSymbolNumber + 1)] = ((sizeLeft >= 0)? len_dict_nonTermSym[sizeLeft] : 1) + ((sizeRight >= 0)? len_dict_nonTermSym[sizeRight] : 1);
				
			}
			
			util::bit_compress(v);
			util::bit_compress(len_dict_nonTermSym);
			construct_im(sizeTextNonTerminalSymbols,len_dict_nonTermSym,0);
			construct_im(perm, v, 0);
		}
		
		void convertRulesDictToIntVector(unordered_map<int,int> &hashReplacement){
			int size = dictionary->size() * 2;
			
			int_vector<> v(size);
			int_vector<> len_dict_nonTermSym(dictionary->size());
			
			int k = 0;
			for(int i = biggerTerminalSymbolNumber + 1; i < nonTerminalSymbolNumber;i++){
				int firstNum = (*dictionary)[i].current_symbol_Num;
				int secondNum = (*dictionary)[i].next_symbol_Num;
				
				unordered_map<int,int>::iterator gotFirst = hashReplacement.find(firstNum);
				unordered_map<int,int>::iterator gotSecond = hashReplacement.find(secondNum);
				if(gotFirst != hashReplacement.end()){
					firstNum = gotFirst->second;
				}
				if(gotSecond != hashReplacement.end()){
					secondNum = gotSecond->second;
				}
				
				v[k] = firstNum;
				k++;
				v[k] = secondNum;
				k++;
				
				int sizeLeft = firstNum - (biggerTerminalSymbolNumber + 1);
				int sizeRight = secondNum - (biggerTerminalSymbolNumber + 1);
				
				int sizeTextNonTerminalSymbol; 
				len_dict_nonTermSym[i - (biggerTerminalSymbolNumber + 1)] = ((sizeLeft >= 0)? len_dict_nonTermSym[sizeLeft] : 1) + ((sizeRight >= 0)? len_dict_nonTermSym[sizeRight] : 1);
				//sizeTextNonTerminalSymbols->push_back(sizeTextNonTerminalSymbol);
			}
			
			util::bit_compress(v);
			util::bit_compress(len_dict_nonTermSym);
			construct_im(sizeTextNonTerminalSymbols,len_dict_nonTermSym,0);
			construct_im(perm, v, 0);
		}
		
		void convertRulesDictToIntVector(){
			int size = dictionary->size() * 2;
			
			int_vector<> v(size);
			int_vector<> len_dict_nonTermSym(dictionary->size());
			
			int k = 0;
			for(int i = biggerTerminalSymbolNumber + 1; i < nonTerminalSymbolNumber;i++){
				int firstNum = (*dictionary)[i].current_symbol_Num;
				int secondNum = (*dictionary)[i].next_symbol_Num;
				v[k] = firstNum;
				k++;
				v[k] = secondNum;
				k++;
				
				int sizeLeft = firstNum - (biggerTerminalSymbolNumber + 1);
				int sizeRight = secondNum - (biggerTerminalSymbolNumber + 1);
				
				
				
				int sizeTextNonTerminalSymbol; 
				len_dict_nonTermSym[i - (biggerTerminalSymbolNumber + 1)] = ((sizeLeft >= 0)? len_dict_nonTermSym[sizeLeft] : 1) + ((sizeRight >= 0)? len_dict_nonTermSym[sizeRight] : 1);
				//sizeTextNonTerminalSymbols->push_back(sizeTextNonTerminalSymbol);
			}
			
			util::bit_compress(v);
			util::bit_compress(len_dict_nonTermSym);
			construct_im(sizeTextNonTerminalSymbols,len_dict_nonTermSym,0);
			construct_im(perm, v, 0);	
		}
		
		void printSecuenceSymbols(){
			unsigned int currentPos = 0;

			while(currentPos < secuenceSymbols->size()){
				cout << (*secuenceSymbols)[currentPos].symbolNumber << " ";
				currentPos = (*secuenceSymbols)[currentPos].nextUsablePos;
			}
			cout << endl;
		}
		
		void printwmatrix(){
			cout <<"wmatrix"<<endl;
			for(int i = 0; i < wmatrix.size();i++){
				cout << wmatrix[i] << " ";
			}
			cout << endl;
		}

		void printPositions(int &a, int &b){
			Pair p;
			p.current_symbol_Num = a;
			p.next_symbol_Num = b;
			unordered_map<Pair,PairFrequency,PairHasher>::const_iterator got = hashTable->find(p);
			if (got != hashTable->end()){
				PairFrequency pf = got->second;
				cout << "Freq:" << pf.frequency << endl;

				int currentPos = pf.first_pair_Pos;
				for(int i = 0; i < pf.frequency; i++){
					cout << currentPos << endl;
					currentPos = (*secuenceSymbols)[currentPos].nextNeighborPair;
				}
			}
		}

		void printSymbolPositions(int &a){
			unordered_map<int,int>::const_iterator got = firstSymbolOcurrencyHash->find(a);
			cout << "Symbol pos:" << a << endl;
			if(got != firstSymbolOcurrencyHash->end()){
				int currentPos = got->second;
				while(currentPos != -1){
					cout << currentPos << endl;
					currentPos = (*secuenceSymbols)[currentPos].nextSymbolPos;
				}
			}

		}

		void printPriorityQueue(){
			cout << "less freqPairs"<<endl;
			for(unsigned int i = 0; i < pQueue->lessFrequentPairs->size(); i++){
				list<PairWithFrequencies>::iterator it = (*(pQueue->lessFrequentPairs))[i].begin();
				while(it != (*(pQueue->lessFrequentPairs))[i].end()){
					cout << "i:" << i << " pair: "<< (*it).current_symbol_Num <<","<< (*it).next_symbol_Num << " pos:"<< (*it).first_pair_Pos << " freq:" << (*it).frequency << endl;
					it++;
				}
			}
			cout << "more freqPairs" << endl;
			list<PairWithFrequencies>::iterator it = (*(pQueue->moreFrequentPairs)).begin();
			while(it != (*(pQueue->moreFrequentPairs)).end()){
				cout << "+ pair: "<< (*it).current_symbol_Num <<","<< (*it).next_symbol_Num << " pos:"<< (*it).first_pair_Pos << " freq:" << (*it).frequency << endl;
				it++;
			}

		}

		void fillAllFromTextFile(string &fileName){
			intAlphabet *ia = new intAlphabet();

			ia->computeIntSymbolsFromFile2(fileName);
			nonTerminalSymbolNumber = ia->numSymbol;
			//cout << ia->numSymbol << endl;
			fillSecuenceSymbolsAndHashTable(*(ia->intSymbols));
			//cout << "lala"<<endl;
			//ia->printHash();
			biggerTerminalSymbolNumber = ia->maxTerminalSymbol;
			delete ia;

			fillPriorityQueue();
			//printPriorityQueue();
		}


		void fillAllStr(string &strSymbols){
			intAlphabet *ia = new intAlphabet();
			ia->computeIntSymbols2(strSymbols);
			nonTerminalSymbolNumber = ia->numSymbol;
			//cout << ia->numSymbol << endl;
			fillSecuenceSymbolsAndHashTable(*(ia->intSymbols));
			//cout << "lala"<<endl;
			//ia->printHash();
			biggerTerminalSymbolNumber = ia->maxTerminalSymbol;
			delete ia;

			fillPriorityQueue();
			printPriorityQueue();
		}

		void writeCRFiles(string &strFile){
			//Writes file C (includes alphabet and positions in C)


			//Compressed String, alphabet and positions, C file
			ofstream cFile(strFile+".c");

			//Writes positions of c in bitvector.
			bit_vector b_positions_c_in_t;
			b_positions_c_in_t = bit_vector(secuenceSymbols->size());
			
			//Writes bigger non terminal and terminal symbol
			cFile << (nonTerminalSymbolNumber - 1) << endl;
			cFile << biggerTerminalSymbolNumber << endl;
			/*
			//Symbols
			unsigned int currentPos = 0;
			while(currentPos < secuenceSymbols->size()){
				cFile << (*secuenceSymbols)[currentPos].symbolNumber;
				b_positions_c_in_t[currentPos] = 1;
				currentPos = (*secuenceSymbols)[currentPos].nextUsablePos;
				if(currentPos < secuenceSymbols->size()){
					cFile << ",";
				}
			}
			cFile << endl;
			*/
			//Serialize bitvector pos C
			sd_vector<> sdb(b_positions_c_in_t);
			cFile << size_in_bytes(sdb) << endl;
			sdb.serialize(cFile);
			cFile << endl;
			//Serialize golynski perms
			cFile << size_in_bytes(perm) << endl;
			perm.serialize(cFile);
			
			cout << "perm: ";
			for(int i = 0; i < perm.size(); i++){
				cout << perm[i] << " ";
			}
			cout << endl;
			
			cFile << endl;
			//Serialize w_matrix grid
			cFile << size_in_bytes(wmatrix) << endl;
			wmatrix.serialize(cFile);
			cFile << endl;
			
			//Serialize both indexes arrays in wmatrix (todo: serialize inside)
			cFile << size_in_bytes(indexPositionsNonTerminals) << endl;
			indexPositionsNonTerminals.serialize(cFile);
			cFile << endl;
			
			
			cFile << size_in_bytes(indexPositionsColumn) << endl;
			indexPositionsColumn.serialize(cFile);
			cFile << endl;
			
			
			//Serialize golynski sizeTextNonTerminalSymbols
			cFile << size_in_bytes(sizeTextNonTerminalSymbols) << endl;
			sizeTextNonTerminalSymbols.serialize(cFile);
			cFile << endl;
			//Serialize w_matrix bitvector
			//sd_vector sdb_wmatrix(bv_wmatrix);
			cFile << size_in_bytes(sdb_wmatrix) << endl;
			sdb_wmatrix.serialize(cFile);
			cFile << endl;
			
			cFile.close();
		}

		
		void replaceMostFrequentPair(){
			PairWithFrequencies pfq = pQueue->getMoreFrequentPair(secuenceSymbols->size());

			Pair p;
			p.current_symbol_Num = pfq.current_symbol_Num;
			p.next_symbol_Num = pfq.next_symbol_Num;
			
			/*if(pfq.frequency == 1)
				std::cout << "FinalPair:" << p.current_symbol_Num << "," << p.next_symbol_Num << " f:"<< pfq.frequency << std::endl;
			*/
			if(pfq.current_symbol_Num == -1 && pfq.next_symbol_Num == -1){
				return;
			}

			(*dictionary)[nonTerminalSymbolNumber] = p;

			//std::cout << "PairDelete before:" << pfq.current_symbol_Num << "," << pfq.next_symbol_Num << " f:"<< pfq.frequency << std::endl;

			int currentPos = (*hashTable)[p].first_pair_Pos;

			unordered_map<Pair,PairFrequency,PairHasher> *newHash = new unordered_map<Pair,PairFrequency,PairHasher>();
			newHash->max_load_factor(0.75);
			

			while(currentPos != -1){
				//Edit neighbors
				//Neighbors
				//cout << "currentPos:"<<currentPos<<" pair:"<<p.current_symbol_Num<<","<<p.next_symbol_Num<<endl;
				int leftNeighborSymbolPos = -1;
				int rightNeighborSymbolPos = -1;

				leftNeighborSymbolPos = (*secuenceSymbols)[currentPos].prevUsablePos;
				rightNeighborSymbolPos = (*secuenceSymbols)[currentPos].nextUsablePos;
				int rightNeighbotSymbolPosPlusOne = -1;
				
				if(rightNeighborSymbolPos != -1){
					rightNeighbotSymbolPosPlusOne = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
				}
				
				if(leftNeighborSymbolPos < 0){
					leftNeighborSymbolPos = -1;
				}
				if(rightNeighborSymbolPos >= secuenceSymbols->size()){
					rightNeighborSymbolPos = -1;
				}

				if((*secuenceSymbols)[currentPos].symbolNumber == p.current_symbol_Num && rightNeighborSymbolPos != -1 && (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber == p.next_symbol_Num){
					/*
					 *  if ab is the pair to replace, let x be the left character of a and let y the right character of b
					 *  then xaby is the context, where pos of a = currentPos, b = currentPos + 1,x = currentPos -1, y = currentPos + 2 (rightNeighborSymbolPos)
					 *
					 *  So:
					 *    pair xa starts at currentPos -1 (leftNeighborSymbolPos)
					 *    pair ab starts at currentPos
					 *    pair by starts at currentPos + 1 (rightNeighborSymbolPos)
					 *
					 *    pair xA starts at currentPos - 1 (leftNeighborSymbolPos)
					 *    pair bY starts at currentPos + 1 (rightNeighborSymbolPos)
					 * */

					//xa (if x is the pos of a non terminal symbol number, means that it was counted in xA previous step, where previous y is current a)
					if(leftNeighborSymbolPos >= 0 && (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber != nonTerminalSymbolNumber){

						Pair xa;
						xa.current_symbol_Num = (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber;
						xa.next_symbol_Num = (*secuenceSymbols)[currentPos].symbolNumber;
						//cout << "xa:"<< xa.current_symbol_Num<<","<<xa.next_symbol_Num<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(xa);
						PairFrequency xaPf;

						if(got != hashTable->end()){
							xaPf.frequency = (*newHash)[xa].frequency - 1;
						}
						else{
							xaPf = (*hashTable)[xa];
							xaPf.frequency = xaPf.frequency - 1;
						}

						if(xaPf.first_pair_Pos == currentPos){
							xaPf.first_pair_Pos = (*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair;
						}
						if(xaPf.last_pair_Pos == currentPos){
							xaPf.last_pair_Pos = (*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair;
						}

						(*newHash)[xa] = xaPf;
						
						if((*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair].nextNeighborPair = (*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair;
						}
						if((*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair].prevNeighborPair = (*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair;
						}
					}

					//by
					if(rightNeighbotSymbolPosPlusOne >= 0 && rightNeighbotSymbolPosPlusOne < secuenceSymbols->size() && (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber != -1){

						Pair by;
						by.current_symbol_Num = (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber;
						by.next_symbol_Num = (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber;
						//cout << "by:"<< by.current_symbol_Num<<","<<by.next_symbol_Num<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(by);
						PairFrequency byPf;

						if(got != hashTable->end()){
							byPf.frequency = (*newHash)[by].frequency - 1;
						}
						else{
							byPf = (*hashTable)[by];
							byPf.frequency = byPf.frequency - 1;
						}

						if(byPf.first_pair_Pos == currentPos){
							byPf.first_pair_Pos = (*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair;
						}
						if(byPf.last_pair_Pos == currentPos){
							byPf.last_pair_Pos = (*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair;
						}

						(*newHash)[by] = byPf;

						if((*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair].nextNeighborPair = (*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair;
						}
						if((*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair].prevNeighborPair = (*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair;
						}
					}
					
					//xA
					if(leftNeighborSymbolPos >= 0 && rightNeighborSymbolPos >= 0 && rightNeighborSymbolPos < secuenceSymbols->size() && rightNeighbotSymbolPosPlusOne >= 0 && rightNeighbotSymbolPosPlusOne < secuenceSymbols->size() && (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber != -1){
						Pair xA;
						xA.current_symbol_Num = (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber;
						xA.next_symbol_Num = nonTerminalSymbolNumber;
						//cout << "xA:"<< xA.current_symbol_Num<<","<<xA.next_symbol_Num<< " pos:"<< leftNeighborSymbolPos << " rsp:"<<rightNeighbotSymbolPosPlusOne<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(xA);
						PairFrequency xAPf;

						if(got != newHash->end()){
							xAPf = (*newHash)[xA];
							
							(*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair = xAPf.last_pair_Pos;
							if(xAPf.last_pair_Pos != -1){
								(*secuenceSymbols)[xAPf.last_pair_Pos].nextNeighborPair = leftNeighborSymbolPos;
							}
							
							xAPf.last_pair_Pos = leftNeighborSymbolPos;
							xAPf.frequency = xAPf.frequency + 1;
						}
						else{
							xAPf.first_pair_Pos = leftNeighborSymbolPos;
							xAPf.last_pair_Pos = leftNeighborSymbolPos;
							xAPf.frequency = 1;
						}
						(*newHash)[xA] = xAPf;
					}
					//Ay
					if(rightNeighbotSymbolPosPlusOne >= 0 && rightNeighbotSymbolPosPlusOne < secuenceSymbols->size() && (*secuenceSymbols)[currentPos].nextNeighborPair != rightNeighbotSymbolPosPlusOne){
						Pair Ay;
						Ay.current_symbol_Num = nonTerminalSymbolNumber;
						Ay.next_symbol_Num = (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber;
						//cout << "Ay:"<< Ay.current_symbol_Num<<","<<Ay.next_symbol_Num<<endl;
						//if(rightNeighborSymbolPos != -1 && pfq.next_symbol_Num == (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber){
							unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(Ay);
							PairFrequency AyPf;

							if(got != newHash->end()){
								AyPf = (*newHash)[Ay];
								
								(*secuenceSymbols)[currentPos].prevNeighborPair = AyPf.last_pair_Pos;
								if(AyPf.last_pair_Pos != -1){
									(*secuenceSymbols)[AyPf.last_pair_Pos].nextNeighborPair = currentPos;
								}
								
								AyPf.last_pair_Pos = currentPos;
								AyPf.frequency = AyPf.frequency + 1;
							}
							else{
								AyPf.first_pair_Pos = currentPos;
								AyPf.last_pair_Pos = currentPos;
								AyPf.frequency = 1;
							}
							(*newHash)[Ay] = AyPf;
						//}
					}
					
					if(rightNeighborSymbolPos >= 0 && rightNeighborSymbolPos < secuenceSymbols->size()){
						(*secuenceSymbols)[currentPos].symbolNumber = nonTerminalSymbolNumber;

						//Update positions of "emulated" linked list
						int prevPosUse = (*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos;
						int nextPosUse = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
						if(prevPosUse > -1 && prevPosUse < secuenceSymbols->size()){
							(*secuenceSymbols)[prevPosUse].nextUsablePos = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
						}
						if(nextPosUse > -1 && nextPosUse < secuenceSymbols->size()){
							(*secuenceSymbols)[nextPosUse].prevUsablePos = (*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos;
						}

						(*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair = -1;
					}
				}
				else if((*secuenceSymbols)[currentPos].symbolNumber != -1 && rightNeighborSymbolPos >= 0 && rightNeighborSymbolPos < secuenceSymbols->size() && (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber != -1){
					//If it's different than -1,-1 and if it's not on new hash table, add it to newhash
					Pair np;
					np.current_symbol_Num = (*secuenceSymbols)[currentPos].symbolNumber;
					np.next_symbol_Num = (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber;
					unordered_map<Pair,PairFrequency,PairHasher>::iterator it2 = newHash->find(np);
					
					if(it2 == newHash->end()){
						PairFrequency npf;
						npf.frequency = 1;
						npf.first_pair_Pos = currentPos;
						npf.last_pair_Pos = currentPos;
						(*newHash)[np] = npf;
					}
					else{
						PairFrequency npf = it2->second;
						npf.frequency = npf.frequency + 1;
						npf.last_pair_Pos = currentPos;
						(*newHash)[np] = npf;
					}
				}
				
				
				int nextPos = (*secuenceSymbols)[currentPos].nextNeighborPair;
				
				currentPos = nextPos;
			}
			(*newHash)[p] = (*hashTable)[p];
			(*newHash)[p].frequency = -1;

			//Adds modifications into priority queue
			pQueue->addChangesToPriorityQueue(newHash, hashTable, secuenceSymbols->size());

			//Adds all frequency modifications from hashTableSymbolsToModify to hashTable
			unordered_map<Pair,PairFrequency,PairHasher>::iterator it = newHash->begin();
			while(it != newHash->end()){
				if((it->second).frequency > 0){
					(*hashTable)[it->first] = it->second;
				}
				else{
					unordered_map<Pair,PairFrequency,PairHasher>::iterator it2 = hashTable->find(it->first);
					if(it2 != hashTable->end())
						hashTable->erase(it2);
				}
				//cout << "mod htable:"<< (it->first).current_symbol_Num << "," << (it->first).next_symbol_Num <<" f:" << (it->second).frequency << endl;
				it++;
			}
			if(newHash->size() > 0)
				nonTerminalSymbolNumber++;
			delete newHash;			
		}
		
		void replaceMostFrequentPair2(){
			PairWithFrequencies pfq = pQueue->getMoreFrequentPair(secuenceSymbols->size());

			Pair p;
			p.current_symbol_Num = pfq.current_symbol_Num;
			p.next_symbol_Num = pfq.next_symbol_Num;
			
			/*if(pfq.frequency == 1)
				std::cout << "FinalPair:" << p.current_symbol_Num << "," << p.next_symbol_Num << " f:"<< pfq.frequency << std::endl;
			*/
			if(pfq.current_symbol_Num == -1 && pfq.next_symbol_Num == -1){
				return;
			}

			(*dictionary)[nonTerminalSymbolNumber] = p;

			//std::cout << "PairDelete before:" << pfq.current_symbol_Num << "," << pfq.next_symbol_Num << " f:"<< pfq.frequency << std::endl;

			int currentPos = (*hashTable)[p].first_pair_Pos;

			unordered_map<Pair,PairFrequency,PairHasher> *newHash = new unordered_map<Pair,PairFrequency,PairHasher>();
			newHash->max_load_factor(0.75);
			

			while(currentPos != -1){
				//Edit neighbors
				//Neighbors
				//cout << "currentPos:"<<currentPos<<" pair:"<<p.current_symbol_Num<<","<<p.next_symbol_Num<<endl;
				int leftNeighborSymbolPos = -1;
				int rightNeighborSymbolPos = -1;
				int rightNeighbotSymbolPosPlusOne = -1;

				leftNeighborSymbolPos = (*secuenceSymbols)[currentPos].prevUsablePos;
				rightNeighborSymbolPos = (*secuenceSymbols)[currentPos].nextUsablePos;				
				
				if(rightNeighborSymbolPos != -1){
					rightNeighbotSymbolPosPlusOne = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
				}
				
				if(leftNeighborSymbolPos < 0){
					leftNeighborSymbolPos = -1;
				}
				if(rightNeighborSymbolPos >= secuenceSymbols->size()){
					rightNeighborSymbolPos = -1;
				}
				if(rightNeighbotSymbolPosPlusOne >= secuenceSymbols->size()){
					rightNeighbotSymbolPosPlusOne = -1;
				}
				
				int nextPos;

				if((*secuenceSymbols)[currentPos].symbolNumber == p.current_symbol_Num && rightNeighborSymbolPos != -1 && (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber == p.next_symbol_Num){
					/*
					 *  if ab is the pair to replace, let x be the left character of a and let y the right character of b
					 *  then xaby is the context, where pos of a = currentPos, b = currentPos + 1,x = currentPos -1, y = currentPos + 2 (rightNeighborSymbolPos)
					 *
					 *  So:
					 *    pair xa starts at currentPos -1 (leftNeighborSymbolPos)
					 *    pair ab starts at currentPos
					 *    pair by starts at currentPos + 1 (rightNeighborSymbolPos)
					 *
					 *    pair xA starts at currentPos - 1 (leftNeighborSymbolPos)
					 *    pair bY starts at currentPos + 1 (rightNeighborSymbolPos)
					 * */

					//xa (if x is the pos of a non terminal symbol number, means that it was counted in Ay previous step, where previous y is current a)
					if(leftNeighborSymbolPos >= 0 && (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber != -1){//&& (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber != nonTerminalSymbolNumber

						Pair xa;
						xa.current_symbol_Num = (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber;
						xa.next_symbol_Num = (*secuenceSymbols)[currentPos].symbolNumber;
						//cout << "xa:"<< xa.current_symbol_Num<<","<<xa.next_symbol_Num<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(xa);
						PairFrequency xaPf;

						if(got != hashTable->end()){
							xaPf.frequency = (*newHash)[xa].frequency - 1;
						}
						else{
							xaPf = (*hashTable)[xa];
							xaPf.frequency = xaPf.frequency - 1;
						}

						if(xaPf.first_pair_Pos == currentPos){
							xaPf.first_pair_Pos = (*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair;
						}
						if(xaPf.last_pair_Pos == currentPos){
							xaPf.last_pair_Pos = (*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair;
						}

						(*newHash)[xa] = xaPf;
						
						if((*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair].nextNeighborPair = (*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair;
						}
						if((*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[leftNeighborSymbolPos].nextNeighborPair].prevNeighborPair = (*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair;
						}
					}

					//by
					if(rightNeighbotSymbolPosPlusOne >= 0 &&  (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber != -1){

						Pair by;
						by.current_symbol_Num = (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber;
						by.next_symbol_Num = (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber;
						//cout << "by:"<< by.current_symbol_Num<<","<<by.next_symbol_Num<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(by);
						PairFrequency byPf;

						if(got != hashTable->end()){
							byPf.frequency = (*newHash)[by].frequency - 1;
						}
						else{
							byPf = (*hashTable)[by];
							byPf.frequency = byPf.frequency - 1;
						}

						if(byPf.first_pair_Pos == currentPos){
							byPf.first_pair_Pos = (*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair;
						}
						if(byPf.last_pair_Pos == currentPos){
							byPf.last_pair_Pos = (*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair;
						}

						(*newHash)[by] = byPf;

						if((*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair].nextNeighborPair = (*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair;
						}
						if((*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair != -1){
							(*secuenceSymbols)[(*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair].prevNeighborPair = (*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair;
						}
					}
					
					//xA
					if(leftNeighborSymbolPos >= 0 && rightNeighborSymbolPos >= 0 && rightNeighbotSymbolPosPlusOne >= 0 && (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber != -1 && (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber != -1){
						Pair xA;
						xA.current_symbol_Num = (*secuenceSymbols)[leftNeighborSymbolPos].symbolNumber;
						xA.next_symbol_Num = nonTerminalSymbolNumber;
						//cout << "xA:"<< xA.current_symbol_Num<<","<<xA.next_symbol_Num<< " pos:"<< leftNeighborSymbolPos << " rsp:"<<rightNeighbotSymbolPosPlusOne<<endl;
						unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(xA);
						PairFrequency xAPf;

						if(got != newHash->end()){
							xAPf = (*newHash)[xA];
							
							(*secuenceSymbols)[leftNeighborSymbolPos].prevNeighborPair = xAPf.last_pair_Pos;
							if(xAPf.last_pair_Pos != -1){
								(*secuenceSymbols)[xAPf.last_pair_Pos].nextNeighborPair = leftNeighborSymbolPos;
							}
							
							xAPf.last_pair_Pos = leftNeighborSymbolPos;
							xAPf.frequency = xAPf.frequency + 1;
						}
						else{
							xAPf.first_pair_Pos = leftNeighborSymbolPos;
							xAPf.last_pair_Pos = leftNeighborSymbolPos;
							xAPf.frequency = 1;
						}
						(*newHash)[xA] = xAPf;
					}
					//Ay
					if(rightNeighbotSymbolPosPlusOne >= 0 &&  (*secuenceSymbols)[currentPos].nextNeighborPair != rightNeighbotSymbolPosPlusOne){
						Pair Ay;
						Ay.current_symbol_Num = nonTerminalSymbolNumber;
						Ay.next_symbol_Num = (*secuenceSymbols)[rightNeighbotSymbolPosPlusOne].symbolNumber;
						//cout << "Ay:"<< Ay.current_symbol_Num<<","<<Ay.next_symbol_Num<<endl;
						//if(rightNeighborSymbolPos != -1 && pfq.next_symbol_Num == (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber){
							unordered_map<Pair,PairFrequency,PairHasher>::iterator got = newHash->find(Ay);
							PairFrequency AyPf;

							if(got != newHash->end()){
								AyPf = (*newHash)[Ay];
								
								(*secuenceSymbols)[currentPos].prevNeighborPair = AyPf.last_pair_Pos;
								if(AyPf.last_pair_Pos != -1){
									(*secuenceSymbols)[AyPf.last_pair_Pos].nextNeighborPair = currentPos;
								}
								
								AyPf.last_pair_Pos = currentPos;
								AyPf.frequency = AyPf.frequency + 1;
							}
							else{
								AyPf.first_pair_Pos = currentPos;
								AyPf.last_pair_Pos = currentPos;
								AyPf.frequency = 1;
							}
							(*newHash)[Ay] = AyPf;
						//}
					}
					
					if(rightNeighborSymbolPos >= 0){
						(*secuenceSymbols)[currentPos].symbolNumber = nonTerminalSymbolNumber;

						//Update positions of "emulated" linked list
						int prevPosUse = (*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos;
						int nextPosUse = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
						if(prevPosUse > -1 && prevPosUse < secuenceSymbols->size()){
							(*secuenceSymbols)[prevPosUse].nextUsablePos = (*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos;
						}
						if(nextPosUse > -1 && nextPosUse < secuenceSymbols->size()){
							(*secuenceSymbols)[nextPosUse].prevUsablePos = (*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos;
						}

						(*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].prevUsablePos = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].nextUsablePos = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].prevNeighborPair = -1;
						(*secuenceSymbols)[rightNeighborSymbolPos].nextNeighborPair = -1;
					}
					nextPos = (*secuenceSymbols)[currentPos].nextNeighborPair;
				}
				else if((*secuenceSymbols)[currentPos].symbolNumber != -1 && rightNeighborSymbolPos >= 0 && (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber != -1){
					//If it's different than -1,-1 and if it's not on new hash table, add it to newhash
					nextPos = (*secuenceSymbols)[currentPos].nextNeighborPair;
					Pair np;
					np.current_symbol_Num = (*secuenceSymbols)[currentPos].symbolNumber;
					np.next_symbol_Num = (*secuenceSymbols)[rightNeighborSymbolPos].symbolNumber;
					unordered_map<Pair,PairFrequency,PairHasher>::iterator it2 = newHash->find(np);
					
					if(it2 == newHash->end()){
						PairFrequency npf;
						npf.frequency = 1;
						npf.first_pair_Pos = currentPos;
						npf.last_pair_Pos = currentPos;
						(*secuenceSymbols)[currentPos].prevNeighborPair = -1;
						(*secuenceSymbols)[currentPos].nextNeighborPair = -1;
						(*newHash)[np] = npf;
					}
					else{
						PairFrequency npf = it2->second;
						npf.frequency = npf.frequency + 1;
						int prevNeighborPair = (*secuenceSymbols)[currentPos].prevNeighborPair;
						(*secuenceSymbols)[currentPos].prevNeighborPair = npf.last_pair_Pos;
						npf.last_pair_Pos = currentPos;
						if(prevNeighborPair >=0){
							(*secuenceSymbols)[prevNeighborPair].nextNeighborPair = currentPos;
						}
						(*newHash)[np] = npf;
					}
				}
				else{
					//Just remove pointers of all, not a valid pair (replaced in previous step, happens when rightNeigborSymbolPos = nextNeighborPair)
					nextPos = (*secuenceSymbols)[currentPos].nextNeighborPair;
					int prevNeighborPair = (*secuenceSymbols)[currentPos].prevNeighborPair;
					int nextNeighborPair = (*secuenceSymbols)[currentPos].nextNeighborPair;
					if(prevNeighborPair >= 0 && nextNeighborPair >= 0){
						(*secuenceSymbols)[prevNeighborPair].nextNeighborPair = nextNeighborPair;
						(*secuenceSymbols)[nextNeighborPair].prevNeighborPair = prevNeighborPair;
					}
					(*secuenceSymbols)[currentPos].nextNeighborPair = -1;
					(*secuenceSymbols)[currentPos].prevNeighborPair = -1;
				}
				
				
				
				
				currentPos = nextPos;
			}
			(*newHash)[p] = (*hashTable)[p];
			(*newHash)[p].frequency = -1;

			//Adds modifications into priority queue
			pQueue->addChangesToPriorityQueue(newHash, hashTable, secuenceSymbols->size());

			//Adds all frequency modifications from hashTableSymbolsToModify to hashTable
			unordered_map<Pair,PairFrequency,PairHasher>::iterator it = newHash->begin();
			while(it != newHash->end()){
				if((it->second).frequency > 0){
					(*hashTable)[it->first] = it->second;
				}
				else{
					unordered_map<Pair,PairFrequency,PairHasher>::iterator it2 = hashTable->find(it->first);
					if(it2 != hashTable->end())
						hashTable->erase(it2);
				}
				//cout << "mod htable:"<< (it->first).current_symbol_Num << "," << (it->first).next_symbol_Num <<" f:" << (it->second).frequency << endl;
				it++;
			}
			if(newHash->size() > 0)
				nonTerminalSymbolNumber++;
			delete newHash;			
		}
		
		void printGmr(wt_gmr<> &perm){
			for(int i = 0; i < perm.size(); i++){
				cout << perm[i]<<" ";
			}
			cout << endl;
		}
		
		void printPositionsBitVector(bit_vector * b_positions_c_in_t){
			for(int i = 0; i < b_positions_c_in_t->size();i++){
				if((*b_positions_c_in_t)[i] == 1){
					cout << i << " ";
				}
			}
			cout << endl;
		}
		
		bool replaceAllPairsFreqOne(){
			//Hard way, step by step replace each pair of frequency one
			Pair p;
			int currentPos = 0;
			int nextPos = (*secuenceSymbols)[currentPos].nextUsablePos;
			
			bool replacedPair = false;
			
			while(currentPos >= 0 &&  currentPos < secuenceSymbols->size() && nextPos >= 0 && nextPos < secuenceSymbols->size()){
				
				p.current_symbol_Num = (*secuenceSymbols)[currentPos].symbolNumber;
				p.next_symbol_Num = -1;
				
				if(nextPos >= 0 && nextPos < secuenceSymbols->size())
					p.next_symbol_Num = (*secuenceSymbols)[nextPos].symbolNumber;
				
				
				int nextUsablePos = (*secuenceSymbols)[nextPos].nextUsablePos;
				
				if(p.current_symbol_Num != -1 && p.next_symbol_Num != -1){
					replacedPair = true;
					(*dictionary)[nonTerminalSymbolNumber] = p;
					(*secuenceSymbols)[currentPos].symbolNumber = nonTerminalSymbolNumber;
					(*secuenceSymbols)[nextPos].symbolNumber = -1;
					
					//TODO: Change positions of pointers prevNeighborPair and prevUsablePos
					int leftPos = (*secuenceSymbols)[currentPos].prevUsablePos; 
					
					if(nextUsablePos != -1 && nextUsablePos < secuenceSymbols->size()){
						(*secuenceSymbols)[nextUsablePos].prevUsablePos = currentPos;
					}
					
					if(nextPos != -1 && nextPos < secuenceSymbols->size()){
						(*secuenceSymbols)[currentPos].nextUsablePos = (*secuenceSymbols)[nextPos].nextUsablePos;						
						(*secuenceSymbols)[nextPos].prevNeighborPair = -1;
						(*secuenceSymbols)[nextPos].nextNeighborPair = -1;
						(*secuenceSymbols)[nextPos].nextUsablePos = -1;
						(*secuenceSymbols)[nextPos].prevUsablePos = -1;
					}
					
					nonTerminalSymbolNumber++;
				}
					
				currentPos = nextUsablePos;
				
				
				nextPos = -1;
				if(nextUsablePos >= 0 &&  nextUsablePos < secuenceSymbols->size())
					nextPos = (*secuenceSymbols)[nextUsablePos].nextUsablePos;
				
				
			}
			return replacedPair;
		}
		
		void replaceMostFrequentPairs(string &textFile){
			//fillAllStr(textFile);
			fillAllFromTextFile(textFile);
			//printSecuenceSymbols();
			//printHashTable();
			//printPositions(3,1);
			//printSymbolPositions(1);
			//printPriorityQueue();
			/*
			while(pQueue->freqMoreFrequentPair > 0){
				PairWithFrequencies pfq = pQueue->getMoreFrequentPair(secuenceSymbols->size());
				cout << pfq.current_symbol_Num << "," << pfq.next_symbol_Num << " f:"<< pfq.frequency << endl;
			}
			*/

			//printSecuenceSymbols();
			//printAllSecuenceSymbols();
			
			
			//cout << "----------START-----------" << endl;
			while(pQueue->freqMoreFrequentPair > 1){
				//cout <<"-------------Iteration: Start----------------"<<endl;
				//printSecuenceSymbols();
				//cout<<"--Queue before---"<<endl;
				//printPriorityQueue();
				//printAllSecuenceSymbols();
				replaceMostFrequentPair2();
				//cout<<"--Queue after---"<<endl;
				//printPriorityQueue();
				//cout << "max:"<< pQueue->freqMoreFrequentPair<<endl;
				//printHashTable();				
				//cout <<"-------------Iteration: End----------------"<<endl;
			}
			
			//A loop that continues until no pairs of freq = 1 remains
			do{
			}
			while(replaceAllPairsFreqOne());	
			
			fillWMatrixGrid();
			
			//cout << "wt_gmr_size:" << perm.size() << endl;
			//printGmr(perm);
			
			//printPriorityQueue();
			cout <<"writing cr symbols"<<endl;
			
			//printSecuenceSymbols();
			
			//printAllSecuenceSymbols();
			//printDictionary();
			writeCRFiles(textFile);
			//printPositionsBitVector();
		}
		
		

};

int main(int argc, char **argv){
	if(argc > 1){
		Repair *rp = new Repair();
		string textFile = argv[1];
		rp->replaceMostFrequentPairs(textFile);
		/*
		for(int i = 0; i < rp->sizeTextNonTerminalSymbols.size();i++){
			cout << (rp->sizeTextNonTerminalSymbols)[i] << " ";
		}
		cout << endl;
		*/
		//rp->printAllSecuenceSymbols();
		rp->printDictionary();
		rp->printDictionary_perm();
		//rp->printwmatrix();
		rp->printGmr(rp->perm);
		delete rp;
	}
	else{
		cout << "Use: repairMoffat word" << endl;
	}
}
