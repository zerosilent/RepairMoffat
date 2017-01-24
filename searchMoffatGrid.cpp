#include <iostream>
#include <vector>
#include <unordered_map>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sstream>
#include <math.h>
#include "Pair.cpp"
#include "wmPairs.cpp"
#include "PairFrequency.cpp"
#include "Record.cpp"
#include "PriorityQueue.cpp"
#include "intAlphabet.cpp"

using namespace std;
using namespace sdsl;

class Depair{
	public:
		vector<int> *secuenceSymbols;
		vector<char> *text;//text
		vector<int> *carray;//Array C
		wt_gmr<> perm;//Permutation for large alphabet
		
		wm_int<> wmatrix;//Grid to store elements
		int_vector<> indexPositionsNonTerminals;//Added by jconcha: keeps index of A[i] of points in wmatrix
		int_vector<> indexPositionsColumn;////Added by jconcha: keeps index of columns of points in wmatrix
		
		wt_gmr<> sizeTextNonTerminalSymbols;
		sd_vector<> sdb_positions_c_in_t;
		sd_vector<> sdb_wmatrix;//Sparse bitvector for restore original coordinates on wmatrix
		int biggerTerminalSymbolNumber;
		int nonTerminalSymbolNumber;
		int currentCharPos;
		ofstream outputFile;

		Depair(){
			secuenceSymbols = new vector<int>();
			text = new vector<char>();
			carray = new vector<int>();
			currentCharPos = 0;
		}

		~Depair(){
			delete secuenceSymbols;
			delete text;
			delete carray;
		}
		
		Pair accessDictionary(int positionDictionary){
			Pair p;
			int element = biggerTerminalSymbolNumber + 1;
			int position = positionDictionary - element;
			if(2*position + 1 < perm.size()){
				p.current_symbol_Num = perm[2*position];
				p.next_symbol_Num = perm[2*position + 1];
			}
			return p;
		}

		void printDictionary(){
			cout<<"dictionary"<<endl;
			
			int k = 0;
			int element = biggerTerminalSymbolNumber + 1;
			while(k + 1 < perm.size()){
				cout << "n:" << element << " ";
				cout <<"pair:" << perm[k] << "," << perm[k+1] << endl;
				k+=2;
				element++;
			}
		}
		
		void printDictionaryChars(){
			cout<<"dictionary"<<endl;
			
			int k = 0;
			int element = biggerTerminalSymbolNumber + 1;
			while(k + 1 < perm.size()){
				cout << "n:" << char(element) << " ";
				cout <<"pair:" << char(perm[k]) << "," << char(perm[k+1]) << endl;
				k+=2;
				element++;
			}
		}
		
		void printwmatrix(){
			cout <<"wmatrix"<<endl;
			for(int i = 0; i < wmatrix.size();i++){
				cout << wmatrix[i] << " ";
			}
			cout << endl;
		}

		void printSecuenceSymbols(){
			int currentPos = 0;

			while(currentPos < secuenceSymbols->size()){
				cout << (*secuenceSymbols)[currentPos] << " ";
				currentPos++;
			}
			cout << endl;
		}
		
		void printCArray(){
			int currentPos = 0;

			while(currentPos < carray->size()){
				cout << (*carray)[currentPos] << " ";
				currentPos++;
			}
			cout << endl;
		} 
		
		void printGmr(wt_gmr<> &perm){
			for(int i = 0; i < perm.size(); i++){
				cout << perm[i]<<" ";
			}
			cout << endl;
		}
		
		void printText(){
			for(int i=0; i < text->size();i++){
				cout << (*text)[i];
			}
		}
		
		void writeText(string textFile){
			ofstream rFile(textFile+".out");
			if(rFile.is_open()){
				//Alphabet :-)
				for(int i=0; i < text->size();i++){
					rFile << (*text)[i];
				}
			}
			
			rFile.close();
		}
		
		void fillSymbolsText(int num_c){
			
			if(num_c <= biggerTerminalSymbolNumber){
				//secuenceSymbols->push_back(num_c);
				char c= num_c;
				outputFile << c;
				//text->push_back(c);
			}
			else{
				Pair p = accessDictionary(num_c);
				fillSymbolsText(p.current_symbol_Num);
				fillSymbolsText(p.next_symbol_Num);
			}
		}
		
		void fillSymbolsText(int num_c, int begin, int end){
			
			if(num_c <= biggerTerminalSymbolNumber){
				//secuenceSymbols->push_back(num_c);
				if(currentCharPos >= begin && currentCharPos <= end){
					char c= num_c;
					outputFile << c;
				}
				currentCharPos++;
				//text->push_back(c);
			}
			else{
				//currentCharPos++;
				Pair p = accessDictionary(num_c);
				fillSymbolsText(p.current_symbol_Num,begin,end);
				fillSymbolsText(p.next_symbol_Num,begin,end);
			}
		}
		
		void fillSymbols(int num_c){
			
			if(num_c <= biggerTerminalSymbolNumber){
				secuenceSymbols->push_back(num_c);
				
				char c= num_c;
				text->push_back(c);
			}
			else{
				Pair p = accessDictionary(num_c);
				fillSymbols(p.current_symbol_Num);
				fillSymbols(p.next_symbol_Num);
			}
		}
		
		void reconstructSecuenceSymbols(){
			for(int i = 0; i < carray->size();i++){
				fillSymbols((*carray)[i]);
			}
		}
		
		void reconstructSecuenceSymbolsText(string textFile){
			outputFile.open(textFile+".out", std::ofstream::out);
			for(int i = 0; i < carray->size();i++){
				fillSymbolsText((*carray)[i]);
			}
			outputFile.close();
		}
		
		void reconstructSecuenceSymbolsTextRange(string textFile, int begin, int end){
			outputFile.open(textFile+".out", std::ofstream::out);
			currentCharPos = 0;
			//rank_support_v<> rb(&sdb_positions_c_in_t);
			sd_vector<>::rank_1_type sdb_rank(&sdb_positions_c_in_t);
			int c_begin;
			int c_end;
			
			if(end > sdb_positions_c_in_t.size()){
				end = sdb_positions_c_in_t.size();
			}
			
			c_begin = sdb_rank(begin);
			c_end = sdb_rank(end);
			cout << "b:" << c_begin << endl;
			cout << "e:" << c_end << endl;
			for(int i = c_begin; i <= c_end;i++){
				fillSymbolsText((*carray)[i],begin,end);
			}
			outputFile.close();
		}
		
		// You could also take an existing vector as a parameter.
		vector<string> split(string str, char delimiter) {
		  vector<string> internal;
		  stringstream ss(str); // Turn the string into a stream.
		  string tok;
		  
		  while(getline(ss, tok, delimiter)) {
			internal.push_back(tok);
		  }
		  
		  return internal;
		}
		
		Pair getElementDictionary(int &key, int &biggerTerminalSymbolNumber){
			Pair p;
			p.current_symbol_Num = perm[2*(key - biggerTerminalSymbolNumber + 1)];
			p.next_symbol_Num = perm[2*(key - biggerTerminalSymbolNumber + 1) + 1];
			
			return p;
		}
		/*
		void readDictionary(int &biggerTerminalSymbolNumber, int &nonTerminalSymbolNumber){
			cout << perm.size() << endl;
			for(int i = 0; i < perm.size() - 1;i++){
				int b = perm[i];
				int c = perm[i+1];
				(*dictionary)[biggerTerminalSymbolNumber + 1 + 2*i].current_symbol_Num = b;
				(*dictionary)[biggerTerminalSymbolNumber + 1 + 2*i + 1].current_symbol_Num = c;
			}
			
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
		}
		*/
		void readPermutationAndCArray(string crFileStr){
			//
			string line;
			ifstream cFile (crFileStr+".c");
			
			if (cFile.is_open()){
				getline(cFile,line);
				//cout << line << endl;
				nonTerminalSymbolNumber = atoi(line.c_str());
				cout << "nonTerminalSymbolNumber:" << nonTerminalSymbolNumber << endl;
				
				getline(cFile,line);
				//cout << line << endl;
				biggerTerminalSymbolNumber = atoi(line.c_str());
				cout << "biggerTerminalSymbolNumber:" << biggerTerminalSymbolNumber << endl;
				
				/*
				getline(cFile,line);
				
				stringstream ss(line); // Turn the string into a stream.
				string tok;
				//cout "line:"<<line << endl;
				while(getline(ss, tok, ',')) {
					int num = atoi(tok.c_str());
					carray->push_back(num);
					//cout << num << endl;
				}*/
				carray->push_back(nonTerminalSymbolNumber);
				
				getline(cFile,line);
				int bytesToRead = atoi(line.c_str());
				
				cout << "bytesToRead:" << bytesToRead << endl;
				
				line = "";
				char c;
				int i = 0;
				//Reads one more to read \n
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				//cout << "line:" << line << endl;
				
				istringstream is(line);
				//cout << "a:"<< line << endl;
				sdb_positions_c_in_t.load(is);
				cout << "sdb_positions_c_in_t:"<<sdb_positions_c_in_t<<endl;
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:" << bytesToRead << endl;
				//cout << "B:" <<bytesToRead << "line:" << line << endl;
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
								
				istringstream is_perm(line);
				
				
				perm.load(is_perm);
				//perm.print();
				//cout << "perm:"<<perm<<endl;
				for(int i = 0; i < perm.size(); i++){
					cout << perm[i] << " ";
				}
				cout << endl;
				
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:" << bytesToRead << endl;
				line = "";
				i = 0;
				
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_wmatrix(line);
				
				wmatrix.load(is_wmatrix);
				cout << "wmatrix:" << wmatrix << endl;
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:" << bytesToRead << endl;
				line = "";
				i = 0;
				
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_indexPosNonTerm(line);
				indexPositionsNonTerminals.load(is_indexPosNonTerm);
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:" << bytesToRead << endl;
				line = "";
				i = 0;
				
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_indexPosCol(line);
				indexPositionsColumn.load(is_indexPosCol);
				
				cout << "wmatrix indexPositionsNonTerminals:"<< indexPositionsNonTerminals << endl;
				
				cout << "wmatrix indexPositionsColumn:"<< indexPositionsColumn << endl;
				
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:"<<bytesToRead << endl;
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_sizeNonTerm(line);
				sizeTextNonTerminalSymbols.load(is_sizeNonTerm);
				//cout << "sizeTextNonTerminalSymbols:" << sizeTextNonTerminalSymbols<<endl;
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				cout << "bytesToRead:" << bytesToRead << endl;
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_sdbwmatrix(line);
				sdb_wmatrix.load(is_sdbwmatrix);
				
				cout << sdb_wmatrix << endl;
				
				cFile.close();
				//cout << "d:" << endl;
			}
		}
		
		void printPositionsBitVector(){
			for(int i = 0; i < sdb_positions_c_in_t.size();i++){
				if(sdb_positions_c_in_t[i] == 1){
					cout << i << " ";
				}
			}
			cout << endl;
		}
		
		void printCoordinatesWmatrix(){
			cout <<"coord wmatrix:"<<endl;
			cout << "size b:"<<sdb_wmatrix.size() << endl;
			cout << sdb_wmatrix << endl;
			cout << "wmatrix size:"<<wmatrix.size() << endl;
			for(int i = 0; i < wmatrix.size();i++){
				//WmPairs.
				WmPair p;
				
				p = WmPairs::getOldCoordinateSdb(i,wmatrix[i],sdb_wmatrix);
				cout << p.row << ","<< p.column << endl;
			}
			cout << endl;
		}
		
		int_vector<> secondaries(int nonTerminal, int offset){
		
		}
		
		Pair getRangeBinarySearchColumn(int begin, int end, string &patternFile){
			Pair p;
			
			return p;
		}
		
		Pair getRangeBinarySearchRow(int begin, int end, string &patternFile){
			Pair p;
			
			return p;
		}
		
		int_vector<> search(string patternFile){
			//Read the pattern from file into a int vector
			vector<int> *pattern = new vector<int>();
			string line;
			char c;
			ifstream pFile (patternFile);
			while(pFile.get(c)){
				int numSymbolChar = c;
				pattern->push_back(numSymbolChar);
			}
			
			//Split pattern in all type of partitions
			int m = pattern->size();
			
			if(m == 1){
				return secondaries(patternFile);
			}
			else if(m > 1){
				for(int t = 0; t < m - 1; t++){
					//TODO: Find range of rows of p<
				}			
			}
			else{
				//Secondary ocurrence
			}
			//Free memory
			delete pattern;
		}
};

int main(int argc, char **argv){
	Depair *dp = new Depair();
	//intAlphabet *ia = new intAlphabet();
	if(argc > 1){
		string textFile = argv[1];//Prefix name with .c and .r files
		dp->readPermutationAndCArray(textFile);
		//dp->printDictionary();
		//dp->printCArray();
		dp->reconstructSecuenceSymbols();
		//dp->printSecuenceSymbols();
		//dp->printText();
		
		
		if(argc == 2){
			dp->reconstructSecuenceSymbolsText(textFile);
		}
		else{
			int begin = atoi(argv[2]);
			int end = atoi(argv[3]);
			cout << begin << endl;
			cout << end << endl;
			dp->reconstructSecuenceSymbolsTextRange(textFile,begin,end);
		}
		//dp->printPositionsBitVector();
		
		//dp->printGmr(dp->perm);
		dp->printDictionary();
		//dp->printwmatrix();
		//dp->printCoordinatesWmatrix();
		//dp->printDictionaryChars();
	}
	else{
		cout << "Use: depairMoffat word" << endl;
	}
	delete dp;
}
