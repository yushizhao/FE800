# include <fstream>
# include <iostream>
# include <stdint.h>
# include <typeinfo> 

using namespace std;

int main(){
  	ifstream bi;
	bi.open("ANU");
	
	uint32_t a;
	const int size = sizeof(a);

	int counter = 0;

	while(bi.read((char*)&a, size)) {
//		bi.read(reinterpret_cast<char *>(&a), sizeof(a));
//		cout << a << "\n";
//		cout << typeid(a).name() << "\n";
		counter++;
	}

	cout << counter << "\n";
	return 0;
}