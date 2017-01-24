#include <iostream>
#include <vector>
#include <unordered_map>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sstream>
#include <math.h>
#include "Pair.cpp"
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
			cout<<"dictionary:"<<endl;
			
			int k = 0;
			int element = biggerTerminalSymbolNumber + 1;
			while(k + 1 < perm.size()){
				cout << "n:" << element << " ";
				cout <<"pair:" << perm[k] << "," << perm[k+1] << endl;
				k+=2;
				element++;
			}
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
		
		
		
		void readPermutationAndCArray(string crFileStr){
			//
			string line;
			ifstream cFile (crFileStr+".c");
			
			if (cFile.is_open()){
				getline(cFile,line);
				//cout << line << endl;
				nonTerminalSymbolNumber = atoi(line.c_str());
				
				getline(cFile,line);
				//cout << line << endl;
				biggerTerminalSymbolNumber = atoi(line.c_str());
				
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
				
				line = "";
				char c;
				int i = 0;
				//Reads one more to read \n
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is(line);
				//cout << "a:"<< line << endl;
				sdb_positions_c_in_t.load(is);
				
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
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
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_wmatrix(line);
				
				wmatrix.load(is_wmatrix);
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_sizeNonTerm(line);
				sizeTextNonTerminalSymbols.load(is_sizeNonTerm);
				
				getline(cFile,line);
				bytesToRead = atoi(line.c_str());
				line = "";
				i = 0;
				while(i <= bytesToRead){
					cFile.get(c);
					line = line + c;
					i++;
				}
				
				istringstream is_sdbwmatrix(line);
				sdb_wmatrix.load(is_sdbwmatrix);
				
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
		/*dp->printPositionsBitVector();
		
		dp->printGmr(dp->perm);
		dp->printDictionary();*/
	}
	else{
		cout << "Use: depairMoffat word" << endl;
	}
	delete dp;
}
