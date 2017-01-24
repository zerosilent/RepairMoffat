#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <stdlib.h>     /* atoi */

using namespace std;

class intAlphabet{
	public:
		vector<int> *intSymbols;
		int numSymbol;
		int maxTerminalSymbol;
		
		intAlphabet(){
			intSymbols = new vector<int>;
		}
		
		vector<char> getChars(){
			vector<char> chars;
			for(unsigned int i = 1; i <= intSymbols->size();i++){
				char c = (*intSymbols)[i];
				chars.push_back(c);
			}
			return chars;
		}
		
		void deleteMemory(){
			intSymbols->clear();
			delete intSymbols;
		}
		
		~intAlphabet(){
			deleteMemory();
		}
		
		void computeIntSymbols2(string &strSymbols){
			numSymbol = 0;
			for(unsigned int i=0; i < strSymbols.size();i++){
				//Fill hashunordered_map value with num ocurrences
				char key = strSymbols[i];
				int numSymbolChar = key;//Ascii number, it will consume a little more space, but will work for all symbols
				intSymbols->push_back(numSymbolChar);
				if(numSymbol < numSymbolChar){
					numSymbol = numSymbolChar;
					maxTerminalSymbol = numSymbolChar;
				}
			}
			numSymbol++;
		}
		
		void computeIntSymbolsFromFile2(string &fileName){
			numSymbol = 1;		
			ifstream myfile(fileName);
			if(myfile.is_open()){
				//Fill hashunordered_map value with num ocurrences
				char key;
				while (myfile.get(key)){          // loop getting single characters
					int numSymbolChar = key;//Ascii number, it will consume a little more space, but will work for all symbols
					intSymbols->push_back(numSymbolChar);
					if(numSymbol < numSymbolChar){
						numSymbol = numSymbolChar;
						maxTerminalSymbol = numSymbolChar;
					}
				}
				myfile.close();
				numSymbol++;
			}
		}
		
		string getStrSymbols(){
			string strSymbols = "";
			for(unsigned int i = 0; i < intSymbols->size(); i++){
				char c = (*intSymbols)[i];
				strSymbols = strSymbols + c;
			}
			return strSymbols;
		}

		void print(){
			for(unsigned int i = 0; i < intSymbols->size(); i++){
				cout<< (*intSymbols)[i] << " ";
			}
			cout << endl;
		}
};
/*
int main(){
	intAlphabet ia;
	ia.computeIntSymbols("abracadabra");
	ia.print();
	cout << "sizeAlphabet:"<<ia.numSymbol<<endl;
	cout << ia.getStrSymbols() << endl;
}
*/
