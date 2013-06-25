#include <iostream>
#include <vector>

using namespace std;

int main()
{

vector <int> array;
int j = 7;
array[0].push_back(j);
array[1].push_back(j);

cout << array[0][0] << array[1][0] << endl;
return 0;
}
