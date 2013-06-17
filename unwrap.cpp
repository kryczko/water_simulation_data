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
double xyzdistance[number_of_atoms*(timesteps - 1)][3];

for (int i = 1; i < timesteps; i ++)
{
	for (int j = 0; j < number_of_atoms; j ++)
	{
		double dx = abs(x[j + i*number_of_atoms] - x[j + (i - 1)*number_of_atoms]);
		double dy = abs(y[j + i*number_of_atoms] - y[j + (i - 1)*number_of_atoms]);
		double dz = abs(z[j + i*number_of_atoms] - z[j + (i - 1)*number_of_atoms]);

		if (dx > 0.5)
		{
			x[j + i*number_of_atoms] = 1.0 - x[j + i*number_of_atoms];
		}

		if (dy > 0.5)
                {
                        y[j + i*number_of_atoms] = 1.0 - y[j + i*number_of_atoms];
                }

		if (dz > 0.5)
                {
                        z[j + i*number_of_atoms] = 1.0 - z[j + i*number_of_atoms];
                }

	}
}

outputfile.open("newcoords.dat");

for (int i = 0; i < number_of_atoms*timesteps; i ++)
{
	outputfile << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
}


inputfile.close();
outputfile.close();

return 0;
}
