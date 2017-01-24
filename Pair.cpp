

class Pair{
	public:
		int current_symbol_Num;
		int next_symbol_Num;
		
		Pair(){
			current_symbol_Num = -1;
			next_symbol_Num = -1;
		}
		Pair(int a, int b){
			current_symbol_Num = a;
			next_symbol_Num = b;
		}

		static int computeBitShifting(int a, int b){
			//(8*(sizeof(int))/2) = 16 bits, because sizeof(int) = 4 bytes * 8 = 32 bits (16 bits is the half)
			return (a << (8*(sizeof(int))/2)) | b;
		}
		
		bool operator<(const Pair &p)  const {
			int a = computeBitShifting(current_symbol_Num,next_symbol_Num);
			int b = computeBitShifting(p.current_symbol_Num,p.next_symbol_Num);
			return a < b;
		}
		
		bool operator==(const Pair &p) const {
			return ((p.current_symbol_Num == current_symbol_Num) && (p.next_symbol_Num == next_symbol_Num));
		}
};

struct PairHasher
{
  int operator()(const Pair& k) const
  {
    return (std::hash<int>()(k.computeBitShifting(k.current_symbol_Num,k.next_symbol_Num)));
  }
};
