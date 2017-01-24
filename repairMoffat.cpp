#include <iostream>
#include <vector>
#include <unordered_map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <math.h>
#include "Pair.cpp"
#include "PairFrequency.cpp"
#include "Record.cpp"
#include "PriorityQueue.cpp"
#include "intAlphabet.cpp"

using namespace std;
using namespace sdsl;

class Repair{
	public:
		vector<Record> *secuenceSymbols;
		unordered_map<Pair,PairFrequency,PairHasher> *hashTable;
		unordered_map<int,Pair> *dictionary;
		unordered_map<int,int> *symbolOcurrencyHash;
		unordered_map<int,int> *firstSymbolOcurrencyHash;
		int biggerTerminalSymbolNumber;
		int moreFrequentPair;
		int nonTerminalSymbolNumber;
		PriorityQueue *pQueue;

		Repair(){
			secuenceSymbols = new vector<Record>();
			hashTable = new unordered_map<Pair,PairFrequency,PairHasher>;
			hashTable->max_load_factor(0.75);
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


		void fillSecuenceSymbolsAndHashTable(vector<int> symbols){
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
			cout<<"dictionary"<<endl;
			for(unordered_map<int,Pair>::iterator it = dictionary->begin(); it != dictionary->end(); ++it) {
			  cout << "n:"<< it->first << " pair:"<< it->second.current_symbol_Num <<","<< it->second.next_symbol_Num << endl;
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

		void printSecuenceSymbols(){
			unsigned int currentPos = 0;

			while(currentPos < secuenceSymbols->size()){
				cout << (*secuenceSymbols)[currentPos].symbolNumber << " ";
				currentPos = (*secuenceSymbols)[currentPos].nextUsablePos;
			}
			cout << endl;
		}

		void printPositions(int a, int b){
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

		void printSymbolPositions(int a){
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

		void fillAllFromTextFile(string fileName){
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


		void fillAllStr(string strSymbols){
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

		void writeCRFiles(string strFile){
			//Writes file C and R (includes alphabet and positions in C)

			//Dictionary, R file
			ofstream rFile(strFile+".r");
			if(rFile.is_open()){
				//Alphabet :-)
				rFile << (nonTerminalSymbolNumber - 1) << endl;
				rFile <<biggerTerminalSymbolNumber<<endl;
				
				for(unordered_map<int,Pair>::iterator it = dictionary->begin(); it != dictionary->end(); ++it) {
				  rFile << it->first << ","<< it->second.current_symbol_Num <<","<< it->second.next_symbol_Num << endl;
				}
			}
			
			rFile.close();

			//Compressed String, alphabet and positions, C file
			ofstream cFile(strFile+".c");

			//Writes positions of c in bitvector.
			bit_vector b_positions_c_in_t;
			b_positions_c_in_t = bit_vector(secuenceSymbols->size());

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
			sd_vector<> sdb(b_positions_c_in_t);
			sdb.serialize(cFile);
			cFile << b_positions_c_in_t.serialize(cFile) << endl;

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
		
		void printPositionsBitVector(bit_vector * b_positions_c_in_t){
			for(int i = 0; i < b_positions_c_in_t->size();i++){
				if((*b_positions_c_in_t)[i] == 1){
					cout << i << " ";
				}
			}
			cout << endl;
		}
		
		void replaceMostFrequentPairs(string textFile){
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
				replaceMostFrequentPair();
				//cout<<"--Queue after---"<<endl;
				//printPriorityQueue();
				//cout << "max:"<< pQueue->freqMoreFrequentPair<<endl;
				//printHashTable();				
				//cout <<"-------------Iteration: End----------------"<<endl;
			}
			printPriorityQueue();
			cout <<"writing cr symbols"<<endl;
			
			printSecuenceSymbols();
			
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
		delete rp;
	}
	else{
		cout << "Use: repairMoffat word" << endl;
	}
}
