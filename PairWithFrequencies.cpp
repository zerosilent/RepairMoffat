class PairWithFrequencies{
	public:
		int current_symbol_Num;
		int next_symbol_Num;
		int first_pair_Pos;
		int frequency;
		
		PairWithFrequencies(){
			current_symbol_Num = -1;
			next_symbol_Num = -1;
			first_pair_Pos = -1;
			frequency = 0;
		}

		PairWithFrequencies(int &current_symbol_Num, int &next_symbol_Num){
			this->current_symbol_Num = current_symbol_Num;
			this->next_symbol_Num = next_symbol_Num;
		}

		PairWithFrequencies(int &current_symbol_Num, int &next_symbol_Num, int &first_pair_Pos, int &frequency){
			this->current_symbol_Num = current_symbol_Num;
			this->next_symbol_Num = next_symbol_Num;
			this->first_pair_Pos = first_pair_Pos;
			this->frequency = frequency;
		}

		static int computeBitShifting(int a, int b){
			return (a << 16) | b;
		}
		
		bool operator<(const Pair &p)  const {
			int a = computeBitShifting(current_symbol_Num,next_symbol_Num);
			int b = computeBitShifting(p.current_symbol_Num,p.next_symbol_Num);
			return a < b;
		}

};
