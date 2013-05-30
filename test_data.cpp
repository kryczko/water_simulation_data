#include <iostream>
#include <fstream>


using namespace std;

int main()
{

ofstream output;

output.open("test.dat");
double x = 1;
for (int i = 0; i < 192*2; i ++)
{
	output << "  " << i*x << "  " << i*x << "  " << i*x << endl;
}

return 0;
}
