#include <iostream>
#include <vector>
#include <time.h>  
#include <limits.h> 
#include "wmPair.cpp"

using namespace std;
using namespace sdsl;

class WmPairs{
	public:
	
		std::vector<WmPair> v_pair;
		int_vector<> indexPositionsColumn;//Indexes of original positions of pairs before being sorted (this vector is the one to index in wmatrix)
		int_vector<> indexPositionsNonTerminals;//Indexes of original positions of nonterminals in pairs before being sorted

		~WmPairs(){
			v_pair.clear();//Deletes elements
			vector<WmPair>().swap(v_pair);//This free memory
		}
		
		void fillindexPositions(){
			indexPositionsColumn.resize(size());
			for(int i = 0; i < size();i++){
				indexPositionsColumn[i] = i;
			}
			indexPositionsNonTerminals.resize(size());
			for(int i = 0; i < size();i++){
				indexPositionsNonTerminals[i] = i;
			}
		}
		
		int size(){
			return v_pair.size();
		}

		void print(){
			int numRows = v_pair.size();
			for(int i = 0; i < numRows; i++){
				cout << v_pair[i].row << ","<<v_pair[i].column << endl;
			}
			cout << endl;
		}
		/*
		 * a < b -> 0
		 * a == b -> 1
		 * a > b -> 2
		 * 
		 * 
		 * */
		bool isTerminal(int &number, int &biggerTerminalSymbolNumber){
			return(number <= biggerTerminalSymbolNumber);
		}
		int compare(int &a, int &b, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			/*if(a == b){
				return 1;
			}*/
			int a_1 = a;
			int a_2 = a;
			int b_1 = b;	
			int b_2 = b;
			
			if(!isTerminal(a,biggerTerminalSymbolNumber)){
				a_1 = (*dictionary)[a].current_symbol_Num;
				a_2 = (*dictionary)[a].next_symbol_Num;
			}
			if(!isTerminal(b,biggerTerminalSymbolNumber)){
				b_1 = (*dictionary)[b].current_symbol_Num;
				b_2 = (*dictionary)[b].next_symbol_Num;
			}
			
			if(isTerminal(a,biggerTerminalSymbolNumber) && isTerminal(b,biggerTerminalSymbolNumber)){
				if(a == b)
					return 1;
				if(a < b)
					return 0;
				if(a > b)
					return 2;
			}
			else if(!(isTerminal(a,biggerTerminalSymbolNumber)) && isTerminal(b,biggerTerminalSymbolNumber)){
				int value = compare(a_1,b,dictionary,biggerTerminalSymbolNumber);
				if(a_1 != b)
					return value;
				else
					return 2;
			}
			else if(isTerminal(a,biggerTerminalSymbolNumber) && !(isTerminal(b,biggerTerminalSymbolNumber))){
				int value = compare(a,b_1,dictionary,biggerTerminalSymbolNumber);
				if(a != b_1)
					return value;
				else
					return 0;
			}
			else{
				int value = compare(a_1,b_1,dictionary,biggerTerminalSymbolNumber);
				if(value == 1){
					return compare(a_2,b_2,dictionary,biggerTerminalSymbolNumber);
				}
				else{
					return value;
				}
			}
		}
		/*
		int compare2(int a, int b, unordered_map<int,Pair> *dictionary, int biggerTerminalSymbolNumber){
						
			if(a == b){
				return 1;
			}
			
			if(a <= biggerTerminalSymbolNumber && b <= biggerTerminalSymbolNumber){
				if(a < b){
					return 0;
				}
				else if(a > b){
					return 2;
				}
			}
			
			int a_1 = a;
			int a_2 = a;
			int b_1 = b;	
			int b_2 = b;
			
			if(a > biggerTerminalSymbolNumber){
				a_1 = (*dictionary)[a].current_symbol_Num;
				a_2 = (*dictionary)[a].next_symbol_Num;
			}
			
			if(b > biggerTerminalSymbolNumber){
				b_1 = (*dictionary)[b].current_symbol_Num;	
				b_2 = (*dictionary)[b].next_symbol_Num;
			}
			
			if(a_1 == b_1){
				if(a_2 == b_2)
					return 1;
				return compare2(a_2,b_2,dictionary,biggerTerminalSymbolNumber);
			}
			//a_1 != b_1
			else{
				if(a_1 < b_1){
					return 0;
				}
				else{
					return 2;
				}
			}		
		}
		*/
		
		//Compare from right to left
		int compareReverse(int &a, int &b, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			/*if(a == b){
				return 1;
			}*/
			int a_1 = a;
			int a_2 = a;
			int b_1 = b;	
			int b_2 = b;
			
			if(!isTerminal(a,biggerTerminalSymbolNumber)){
				a_1 = (*dictionary)[a].current_symbol_Num;
				a_2 = (*dictionary)[a].next_symbol_Num;
			}
			if(!isTerminal(b,biggerTerminalSymbolNumber)){
				b_1 = (*dictionary)[b].current_symbol_Num;
				b_2 = (*dictionary)[b].next_symbol_Num;
			}
			
			if(isTerminal(a,biggerTerminalSymbolNumber) && isTerminal(b,biggerTerminalSymbolNumber)){
				if(a == b)
					return 1;
				if(a < b)
					return 0;
				if(a > b)
					return 2;
			}
			else if(!(isTerminal(a,biggerTerminalSymbolNumber)) && isTerminal(b,biggerTerminalSymbolNumber)){
				int value = compare(a_2,b,dictionary,biggerTerminalSymbolNumber);
				if(a_2 != b)
					return value;
				else
					return 2;
			}
			else if(isTerminal(a,biggerTerminalSymbolNumber) && !(isTerminal(b,biggerTerminalSymbolNumber))){
				int value = compare(a,b_2,dictionary,biggerTerminalSymbolNumber);
				if(a != b_2)
					return value;
				else
					return 0;
			}
			else{
				int value = compare(a_2,b_2,dictionary,biggerTerminalSymbolNumber);
				if(value == 1){
					return compare(a_1,b_1,dictionary,biggerTerminalSymbolNumber);
				}
				else{
					return value;
				}
			}
		}
		
		void sort_columnOnlyByNumber(unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			/* initialize random seed: */
			srand (time(NULL));
			quicksort(0,v_pair.size()-1,2,dictionary,biggerTerminalSymbolNumber);
		}
		
		void sort(unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			/* initialize random seed: */
			srand (time(NULL));
			quicksort(0,v_pair.size()-1,0,dictionary,biggerTerminalSymbolNumber);
			quicksort(0,v_pair.size()-1,1,dictionary,biggerTerminalSymbolNumber);
		}
		
		void insert(int row, int col){
			WmPair p;
			p.row = row;
			p.column = col;
			v_pair.push_back(p);
		}
		
		void insert(int row, int column, int nonTerminalSymbol){
			WmPair p;
			p.row = row;
			p.column = column;
			p.nonTerminalSymbol = nonTerminalSymbol;
			v_pair.push_back(p);
		}

		void swap_NoIdx(int i, int j){
			int r = v_pair[i].row;
			int c = v_pair[i].column;
			int nts = v_pair[i].nonTerminalSymbol;
			//int idx_Pos_i = indexPositionsColumn[i];
			int idx_NonTerm_Pos_i = indexPositionsNonTerminals[i];
			/*
			int pwValue = pointWeights[i];
			
			pointWeights[i] = pointWeights[j];
			pointWeights[j] = pwValue;
			*/
			//indexPositionsColumn[i] = indexPositionsColumn[j];
			//indexPositionsColumn[j] = idx_Pos_i;
			
			indexPositionsNonTerminals[i] = indexPositionsNonTerminals[j];
			indexPositionsNonTerminals[j] = idx_NonTerm_Pos_i;
			
			v_pair[i].row = v_pair[j].row;
			v_pair[i].column = v_pair[j].column;
			v_pair[i].nonTerminalSymbol = v_pair[j].nonTerminalSymbol;
			
			v_pair[j].row = r;
			v_pair[j].column = c;
			v_pair[j].nonTerminalSymbol = nts;
		}
		
		void swap(int i,int j){
			int r = v_pair[i].row;
			int c = v_pair[i].column;
			int nts = v_pair[i].nonTerminalSymbol;
			int idx_Pos_i = indexPositionsColumn[i];
			int idx_NonTerm_Pos_i = indexPositionsNonTerminals[i];
			/*
			int pwValue = pointWeights[i];
			
			pointWeights[i] = pointWeights[j];
			pointWeights[j] = pwValue;
			*/
			indexPositionsColumn[i] = indexPositionsColumn[j];
			indexPositionsColumn[j] = idx_Pos_i;
			
			indexPositionsNonTerminals[i] = indexPositionsNonTerminals[j];
			indexPositionsNonTerminals[j] = idx_NonTerm_Pos_i;
			
			v_pair[i].row = v_pair[j].row;
			v_pair[i].column = v_pair[j].column;
			v_pair[i].nonTerminalSymbol = v_pair[j].nonTerminalSymbol;
			
			v_pair[j].row = r;
			v_pair[j].column = c;
			v_pair[j].nonTerminalSymbol = nts;
		}
		
		int getNumberColumns(){
			int numCol = 0;
			for(int i = 0; i < v_pair.size();i++){
				if(numCol < v_pair[i].column)
					numCol = v_pair[i].column;
			}
			return numCol;
		}
		
		int getNumberRows(){
			int numRows = 0;
			for(int i = 0; i < v_pair.size();i++){
				if(numRows < v_pair[i].row)
					numRows = v_pair[i].row;
			}
			return numRows;
		}
		
		int partition_inplace_column_byNumber(int &pivot, int &begin, int &end, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){		
			pivot = v_pair[end].column;
			int row_pivot = v_pair[end].row;
			
			int i = begin;
			for(int j = begin; j < end; j++){
				//v_pair[j].column < pivot
				if(v_pair[j].column < pivot){
					swap_NoIdx(i,j);
					i++;
				}
				//A column with same value as pivot, but with smaller row
				//v_pair[j].column == pivot && v_pair[j].row < row_pivot
				else if(v_pair[j].column == pivot && v_pair[j].row < row_pivot){
					swap_NoIdx(i,j);
					i++;
				}
			}
			swap_NoIdx(i,end);	
			
			return i;		
		}
		
		int partition_inplace_column(int &pivot, int &begin, int &end, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){		
			pivot = v_pair[end].column;
			int row_pivot = v_pair[end].row;
			
			int i = begin;
			for(int j = begin; j < end; j++){
				//v_pair[j].column < pivot
				if(compare(v_pair[j].column,pivot,dictionary,biggerTerminalSymbolNumber) == 0){
					swap(i,j);
					i++;
				}
				//A column with same value as pivot, but with smaller row
				//v_pair[j].column == pivot && v_pair[j].row < row_pivot
				else if(compare(v_pair[j].column,pivot,dictionary,biggerTerminalSymbolNumber) == 1 && compare(v_pair[j].row,row_pivot,dictionary,biggerTerminalSymbolNumber) == 0){
					swap(i,j);
					i++;
				}
			}
			swap(i,end);	
			
			return i;		
		}

		//TODO: This should be ordered with compare reverse
		int partition_inplace_row(int &pivot, int &begin, int &end, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){		
			pivot = v_pair[end].row;
			int col_pivot = v_pair[end].column;
			
			int i = begin;
			for(int j = begin; j < end; j++){
				//v_pair[j].row < pivot
				if(compare(v_pair[j].row,pivot,dictionary,biggerTerminalSymbolNumber) == 0){
					swap(i,j);
					i++;
				}
				//A row with same value as pivot, but with smaller column
				//v_pair[j].row == pivot && v_pair[j].column < col_pivot
				else if(compare(v_pair[j].row,pivot,dictionary,biggerTerminalSymbolNumber) == 1 && compare(v_pair[j].column,col_pivot,dictionary,biggerTerminalSymbolNumber) == 0){
					swap(i,j);
					i++;
				}
			}
			swap(i,end);	
			
			return i;		
		}
		
		void quicksort(int begin, int end, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			quicksort(begin,end,1,dictionary,biggerTerminalSymbolNumber);
		}
		
		void quicksort(int begin, int end, int type, unordered_map<int,Pair> *dictionary, int &biggerTerminalSymbolNumber){
			if(begin < end){				
				int pivot;
				
				switch(type){
					case 0:
						pivot = partition_inplace_row(pivot,begin,end,dictionary,biggerTerminalSymbolNumber);
						break;
					case 1:
						pivot = partition_inplace_column(pivot,begin,end,dictionary,biggerTerminalSymbolNumber);
						break;
					case 2:
						pivot = partition_inplace_column_byNumber(pivot,begin,end,dictionary,biggerTerminalSymbolNumber);
						break;
				}
				 
				quicksort(begin,pivot-1,type,dictionary,biggerTerminalSymbolNumber);
				quicksort(pivot+1,end,type,dictionary,biggerTerminalSymbolNumber);
				//print();				
			}			
		}
		
		void renumberPairs2(int biggerTerminalSymbolNumber){
			int minRow = INT_MAX;
			
			for(int i = 0; i < v_pair.size();i++){
				if(minRow > v_pair[i].row)
					minRow = v_pair[i].row;
			}
			
			for(int i = 0; i < v_pair.size();i++){
				v_pair[i].row = i+minRow;
			}
		}
		
		void renumberPairs(){
			for(int i = 0; i < v_pair.size();i++){
				v_pair[i].column = i+1;
			}
			/*int minCol = INT_MAX;
			
			for(int i = 0; i < v_pair.size();i++){
				if(minCol > v_pair[i].column)
					minCol = v_pair[i].column;
			}
			
			for(int i = 0; i < v_pair.size();i++){
				v_pair[i].column = i+minCol;
			}*/
		}
		
		static WmPair getNewNumeration(int x1, int x2, bit_vector &b){
			WmPair p;
			//rank_support_v<1> b_rank(&b);
			bit_vector::select_1_type b_sel(&b);
			//Another one should be added, because sdsl begins from 0 and not from 1
			p.row = b_sel(x1) - x1 + 1 + 1;
			p.column = b_sel(x2 + 1) - (x2 + 1) + 1;
			return p;
		} 
		//-- DEPRECATED, use the one with sdb instead (very sparse bitvector)
		static WmPair getOldCoordinate(int row, int column, bit_vector &b){
			WmPair p;
			//rank_support_v<1> b_rank(&b);
			bit_vector::select_0_type b_sel(&b);
			//Another one should be added, because sdsl begins from 0 and not from 1 
			p.row = b_sel(row) - row + 1;
			p.column = column;
			return p;
		} 
        
        static WmPair getNewNumerationSdb(int x1, int x2, sd_vector<> &b){
			WmPair p;
			//rank_support_v<1> b_rank(&b);
			sd_vector<>::select_1_type b_sel(&b);
			//Another one should be added, because sdsl begins from 0 and not from 1
			p.row = b_sel(x1) - x1 + 1 + 1;
			p.column = b_sel(x2 + 1) - (x2 + 1) + 1;
			return p;
		} 
		
		static WmPair getOldCoordinateSdb(int row, int column, sd_vector<> &b){
			WmPair p;
			//rank_support_v<1> b_rank(&b);
			sd_vector<>::select_0_type b_sel(&b);
			//Another one should be added, because sdsl begins from 0 and not from 1
			p.row = b_sel(row) - row + 1;
			p.column = column;
			return p;
		} 
		
		void printPosZeroesBitVectorB(bit_vector &b){
			for(int i = 0; i < b.size();i++){
				if(b[i] == 0)
					cout << i << " ";
			}
			cout << endl;
		}
		
		void fillBitvectorB(bit_vector &b){
			int bitVectorPointer = 0;

			unordered_map<int,int> hash;
			//Fill hashmap value with num ocurrences
			for(int i = 0; i < v_pair.size();i++){

				int key = v_pair[i].column;
				unordered_map<int,int>::const_iterator got = hash.find(key);
				if (got == hash.end()){
					hash[key] = 1;
				}
				else{
					hash[key] = got->second + 1;
				}
			}

			//Fill bitvector with ocurrences
			int numColumns = v_pair.size();
			
			for(int i = 0; i < v_pair.size();i++){
				if(numColumns < v_pair[i].column)
					numColumns = v_pair[i].column;
			}
			
			int k = 0;

			for(int i = 1; i <= numColumns; i++){
				unordered_map<int,int>::const_iterator got = hash.find(i);
				b[k] = 1;
				k++;
				if (got != hash.end()){
					int zeroes = got->second;//Total number of zeroes is total number of points
					for(int j = 0; j < zeroes; j++){
						b[k] = 0;
						k++;
					}
					if(zeroes == 0){
						b[k] = 1;
						k++;
					}
				}
			}


			b.resize(k);
		}
		
		//Deprecated
		void fillBitvectorBOld(bit_vector &b){
			int bitVectorPointer = 0;

			unordered_map<int,int> hash;
			//Fill hashmap value with num ocurrences
			int numColumns = 0;
			
			for(int i = 0; i < v_pair.size();i++){

				int key = v_pair[i].column;
				unordered_map<int,int>::const_iterator got = hash.find(key);
				if (got == hash.end()){
					hash[key] = 1;
				}
				else{
					hash[key] = got->second + 1;
				}
				if(numColumns < key)
					numColumns = key;
			}
			
			//Print hashmap
			/*cout << "hash:" << endl;
			
			for(unordered_map<int,int>::const_iterator it = hash.begin(); it != hash.end();it++){
				cout << it->first << ","<< it->second << endl;
			}*/
			

			//Fill bitvector with ocurrences
			//int numColumns = v_pair[v_pair.size()-1].column;
			
			//cout << numColumns << endl;

			int k = 0;

			for(int i = 1; i <= numColumns; i++){
				unordered_map<int,int>::const_iterator got = hash.find(i);
				b[k] = 1;
				k++;
				if (got != hash.end()){
					int zeroes = got->second;//Total number of zeroes is total number of points
					for(int j = 0; j < zeroes; j++){
						b[k] = 0;
						k++;
					}
				}
			}


			b.resize(k);
		}
};

