#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include "PairWithFrequencies.cpp"



using namespace std;

class PriorityQueue{
	public:
		std::vector< std::list<PairWithFrequencies> > *lessFrequentPairs;
		std::list<PairWithFrequencies> *moreFrequentPairs;
		int freqMoreFrequentPair;
		
		PriorityQueue(){
			lessFrequentPairs = new std::vector< std::list<PairWithFrequencies> >();
			moreFrequentPairs = new std::list<PairWithFrequencies>();
			freqMoreFrequentPair = 0;
		}
		
		PriorityQueue(int n){
			lessFrequentPairs = new std::vector< std::list<PairWithFrequencies> >();
			for(int i = 0; i < n; i++){
				std::list<PairWithFrequencies> *p = new std::list<PairWithFrequencies>();
				lessFrequentPairs->push_back(*p);
			}
			
			moreFrequentPairs = new std::list<PairWithFrequencies>();
		}
		
		~PriorityQueue(){
			delete(moreFrequentPairs);
			delete(lessFrequentPairs);
		}
		
		void deleteFromPriorityQueue(PairWithFrequencies &pairToDelete, int n){
			std::list<PairWithFrequencies>::iterator it;
			
			std::cout << "PairDelete:" << pairToDelete.current_symbol_Num << "," << pairToDelete.next_symbol_Num << std::endl;
			if(pairToDelete.frequency < ceil(sqrt(n)) - 1){
				it = (*lessFrequentPairs)[pairToDelete.frequency].begin();
				bool found = false;
				while(found == false && it != (*lessFrequentPairs)[pairToDelete.frequency].end()){
					if(it->current_symbol_Num == pairToDelete.current_symbol_Num && it->next_symbol_Num == pairToDelete.next_symbol_Num){
						found = true;
						(*lessFrequentPairs)[pairToDelete.frequency].erase(it);
					}
					else{
						it++;
					}
					
					//std::cout << "it pfq:" << it->current_symbol_Num << "," << it->next_symbol_Num << std::endl;
				}	
			}
			else{
				it = moreFrequentPairs->begin();
				bool found = false;
				while(found == false && it != moreFrequentPairs->end()){
					if(it->current_symbol_Num == pairToDelete.current_symbol_Num && it->next_symbol_Num == pairToDelete.next_symbol_Num){
						found = true;
						(*lessFrequentPairs)[pairToDelete.frequency].erase(it);
					}
					else{
						it++;
					}
					
					//std::cout << "it pfq:" << it->current_symbol_Num << "," << it->next_symbol_Num << std::endl;
				}
			}			
			updateFreqMoreFreqPair(n);
		}
		
		void addChangesToPriorityQueue(std::unordered_map<Pair,PairFrequency,PairHasher> *hashTableSymbolsToModify, std::unordered_map<Pair,PairFrequency,PairHasher> *hashTable, int n){
			updateFreqMoreFreqPair(n);
			std::unordered_map<Pair,PairFrequency,PairHasher>::iterator it = hashTableSymbolsToModify->begin();
			while(it != hashTableSymbolsToModify->end()){
				int freq = (it->second).frequency;
				int originalFreq = -1;
				
				originalFreq = (*hashTable)[it->first].frequency;
				
				//cout << "lala:"<<(it->first).current_symbol_Num <<","<<(it->first).next_symbol_Num << " f:"<<originalFreq << " newf:" << freq << " limit:"<<(ceil(sqrt(n)) - 1)<<endl;
				
				PairWithFrequencies pfq;
				pfq.current_symbol_Num = (it->first).current_symbol_Num;
				pfq.next_symbol_Num = (it->first).next_symbol_Num;
				pfq.first_pair_Pos = (it->second).first_pair_Pos;
				//pfq.last_pair_Pos = (it->second).last_pair_Pos;
				pfq.frequency = (it->second).frequency;
				//std::cout << "hashMod:"<<pfq.current_symbol_Num<<","<<pfq.next_symbol_Num<<" f:"<<pfq.frequency << std::endl;
				
				if(originalFreq < ceil(sqrt(n)) - 1){					
					if(freq > 0){
						//Insert it or modify it
						std::list<PairWithFrequencies>::iterator it2 = (*lessFrequentPairs)[originalFreq].begin();
						bool find = false;
						while(find == false && it2 != (*lessFrequentPairs)[originalFreq].end()){
							if(it2->current_symbol_Num == pfq.current_symbol_Num && it2->next_symbol_Num == pfq.next_symbol_Num){
								find = true;
							}
							else{
								it2++;
							}
						}
						if(find == true){
							(*lessFrequentPairs)[originalFreq].erase(it2);
						}
						//Check if new frequency is greater than ceil(sqrt(n)) - 1
						if(freq < ceil(sqrt(n)) - 1){
							(*lessFrequentPairs)[freq].push_back(pfq);
						}
						else{
							moreFrequentPairs->push_back(pfq);
						}
						if(freq > freqMoreFrequentPair){
							freqMoreFrequentPair = freq;
						}
					}
					else{
						//Delete it
						std::list<PairWithFrequencies>::iterator it2 = (*lessFrequentPairs)[originalFreq].begin();
						bool find = false;
						while(find == false && it2 != (*lessFrequentPairs)[originalFreq].end()){
							if(it2->current_symbol_Num == pfq.current_symbol_Num && it2->next_symbol_Num == pfq.next_symbol_Num){
								(*lessFrequentPairs)[originalFreq].erase(it2);
								find = true;
							}
							else{
								it2++;
							}
						}
					}			
				}
				else{
					if(freq > 0){
						//Insert it or modify it
						std::list<PairWithFrequencies>::iterator it2 = moreFrequentPairs->begin();
						bool find = false;
						while(find == false && it2 != moreFrequentPairs->end()){
							if(it2->current_symbol_Num == pfq.current_symbol_Num && it2->next_symbol_Num == pfq.next_symbol_Num){
								find = true;
							}
							else{
								it2++;
							}
						}
						if(find == true){
							moreFrequentPairs->erase(it2);							
						}
						//Move it if neccesary to lessFreqPairs where correspond
						if(freq < ceil(sqrt(n)) - 1){
							(*lessFrequentPairs)[freq].push_back(pfq);
						}
						else{
							moreFrequentPairs->push_back(pfq);
						}
						
						if(freq > freqMoreFrequentPair){
							freqMoreFrequentPair = freq;
						}
					}
					else{
						//Delete it
						std::list<PairWithFrequencies>::iterator it2 = moreFrequentPairs->begin();
						bool find = false;
						while(find == false && it2 != moreFrequentPairs->end()){
							if(it2->current_symbol_Num == pfq.current_symbol_Num && it2->next_symbol_Num == pfq.next_symbol_Num){
								moreFrequentPairs->erase(it2);
								find = true;
							}
							else{
								it2++;
							}
						}
					}
				}
				
				it++;
			}
			updateFreqMoreFreqPair(n);
		}
		
		void updateFreqMoreFreqPair(int n){
			std::list<PairWithFrequencies>::iterator it;
			
			if(freqMoreFrequentPair < ceil(sqrt(n)) - 1){
				//Checks current frequencies
				int maxFreq = 0;
				if(freqMoreFrequentPair >= 0){
					it = (*lessFrequentPairs)[freqMoreFrequentPair].begin();
					while(it != (*lessFrequentPairs)[freqMoreFrequentPair].end()){
						if(it->frequency > maxFreq){
							maxFreq = it->frequency;
						}
						it++;
					}
					
					if((*lessFrequentPairs)[freqMoreFrequentPair].size() == 0 || maxFreq == 1){
						freqMoreFrequentPair--;
					}
				}
			}
			else{
				int maxFreq = 0;
				if(freqMoreFrequentPair >= 0){
					it = moreFrequentPairs->begin();
					while(it != moreFrequentPairs->end()){
						if(it->frequency > maxFreq){
							maxFreq = it->frequency;
						}
						it++;
					}
					
					if(moreFrequentPairs->size() == 0 || maxFreq == 1){
						freqMoreFrequentPair = ceil(sqrt(n)) - 2;
					}
					else{
						freqMoreFrequentPair = maxFreq;
					}
				}				
			}
		}
		
		PairWithFrequencies getMoreFrequentPair(int n){
			PairWithFrequencies pf;
			int maxFreq = 0;
			//Looks in lessFrequentPairs
			
			if(freqMoreFrequentPair < ceil(sqrt(n)) - 1){

				if(freqMoreFrequentPair > 0 && (*lessFrequentPairs)[freqMoreFrequentPair].size() > 0){
					std::list<PairWithFrequencies>::iterator itDelete;
				
					std::list<PairWithFrequencies>::iterator it = (*lessFrequentPairs)[freqMoreFrequentPair].begin();
					
					while(it != (*lessFrequentPairs)[freqMoreFrequentPair].end()){
						if((*it).frequency > maxFreq){
							maxFreq = (*it).frequency;
							itDelete = it;
						}
						it++;
					}					
					
					if(maxFreq > 0){
						pf = (*itDelete);
						moreFrequentPairs->erase(itDelete);
					}
				}
				
			}
			//Looks in moreFrequentPairs
			else{		
				std::list<PairWithFrequencies>::iterator itDelete;
				
				std::list<PairWithFrequencies>::iterator it = moreFrequentPairs->begin();
				
				while(it != moreFrequentPairs->end()){
					if((*it).frequency > maxFreq){
						maxFreq = (*it).frequency;
						itDelete = it;
					}
					it++;
				}
				
				if(maxFreq > 0){
					pf = (*itDelete);
					moreFrequentPairs->erase(itDelete);
				}
			}
			updateFreqMoreFreqPair(n);
			return pf;
		}
};
