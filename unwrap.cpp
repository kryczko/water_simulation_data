#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;


int main()
{

ifstream inputfile;
ofstream outputfile;

string infile;
int timesteps, number_of_atoms;
double lattice;


cout << "Please enter your file: ";
cin >> infile;
cout << "Please enter the number of atoms: ";
cin >> number_of_atoms;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
cout << "Program running...please wait a moment.\n\n";

inputfile.open(infile.c_str());

double x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!inputfile.eof())
{
        inputfile >> x[counter] >> y[counter] >> z[counter];
        counter ++;

}
double dx[number_of_atoms*timesteps], dy[number_of_atoms*timesteps], dz[number_of_atoms*timesteps]; 



for (int i = 1; i < timesteps; i ++)
{
	for (int j = 0; j < number_of_atoms; j ++)
	{
		dx[j] = 0;
		dy[j] = 0;
		dz[j] = 0;

		double dRx = x[j + i*number_of_atoms] - x[j + (i - 1)*number_of_atoms];
		double dRy = y[j + i*number_of_atoms] - y[j + (i - 1)*number_of_atoms];
		double dRz = z[j + i*number_of_atoms] - z[j + (i - 1)*number_of_atoms];

		dx[j + i*number_of_atoms] = dx[j + (i - 1)*number_of_atoms] + (dRx - round(dRx));
		dy[j + i*number_of_atoms] = dy[j + (i - 1)*number_of_atoms] + (dRy - round(dRy));
		dz[j + i*number_of_atoms] = dz[j + (i - 1)*number_of_atoms] + (dRz - round(dRz));

		x[j + i*number_of_atoms] = x[j] + dx[j + i*number_of_atoms];
		y[j + i*number_of_atoms] = y[j] + dy[j + i*number_of_atoms];
		z[j + i*number_of_atoms] = z[j] + dz[j + i*number_of_atoms];

	}
}

outputfile.open("unwrapped.dat");

for (int i = 0; i < number_of_atoms*timesteps; i ++)
{
	outputfile << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
}


inputfile.close();
outputfile.close();

return 0;
}
