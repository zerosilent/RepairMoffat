#include <iostream>
#include <vector>
#include <unordered_map>
#include <sdsl/sd_vector.hpp>
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
		unordered_map<int,Pair> *dictionary;
		sd_vector<> sdb_positions_c_in_t;
		int biggerTerminalSymbolNumber;
		int nonTerminalSymbolNumber;
		int currentCharPos;
		ofstream outputFile;

		Depair(){
			secuenceSymbols = new vector<int>();
			text = new vector<char>();
			carray = new vector<int>();
			dictionary = new unordered_map<int,Pair>;
			currentCharPos = 0;
		}

		~Depair(){
			delete secuenceSymbols;
			delete text;
			delete carray;
			delete dictionary;
		}

		void printDictionary(){
			cout<<"dictionary"<<endl;
			for(unordered_map<int,Pair>::iterator it = dictionary->begin(); it != dictionary->end(); ++it) {
			  cout << "n:"<< it->first << " pair:"<< it->second.current_symbol_Num <<","<< it->second.next_symbol_Num << endl;
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
				fillSymbolsText((*dictionary)[num_c].current_symbol_Num);
				fillSymbolsText((*dictionary)[num_c].next_symbol_Num);
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
				fillSymbolsText((*dictionary)[num_c].current_symbol_Num,begin,end);
				fillSymbolsText((*dictionary)[num_c].next_symbol_Num,begin,end);
			}
		}
		
		void fillSymbols(int num_c){
			
			if(num_c <= biggerTerminalSymbolNumber){
				secuenceSymbols->push_back(num_c);
				
				char c= num_c;
				text->push_back(c);
			}
			else{
				fillSymbols((*dictionary)[num_c].current_symbol_Num);
				fillSymbols((*dictionary)[num_c].next_symbol_Num);
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
		
		void readDictionaryAndCArray(string crFileStr){
			//
			string line;
			ifstream rFile (crFileStr+".r");
			if (rFile.is_open()){
				getline(rFile,line);
				//cout << line << endl;
				nonTerminalSymbolNumber = atoi(line.c_str());
				
				getline(rFile,line);
				//cout << line << endl;
				biggerTerminalSymbolNumber = atoi(line.c_str());
				
				
				while (getline (rFile,line)){
					Pair p;
					int key;
					vector<string> strSplit = split(line,',');
					if(strSplit.size() > 2){
						key = atoi(strSplit[0].c_str());
						p.current_symbol_Num = atoi(strSplit[1].c_str());
						p.next_symbol_Num = atoi(strSplit[2].c_str());
					}
					(*dictionary)[key]=p;
					//cout << line << endl;
				}
				rFile.close();
			}
			
			ifstream cFile (crFileStr+".c");
			
			if (cFile.is_open()){
				getline(cFile,line);
				
				stringstream ss(line); // Turn the string into a stream.
				string tok;
				//cout "line:"<<line << endl;
				while(getline(ss, tok, ',')) {
					int num = atoi(tok.c_str());
					carray->push_back(num);
					//cout << num << endl;
				}
				
				char c;
				line = "";
				while(cFile.get(c)){
					line = line + c;
				}
				//getline(cFile,line);
				istringstream is(line);
				sdb_positions_c_in_t.load(is);
				
				cFile.close();
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
		dp->readDictionaryAndCArray(textFile);
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
		dp->printPositionsBitVector();
	}
	else{
		cout << "Use: depairMoffat word" << endl;
	}
	delete dp;
}
