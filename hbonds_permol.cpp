#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

int pbc_round(double input)
{
	int i  = input;

	if (abs(input - i) >= 0.5)
	{
		if (input > 0) {i += 1;}
		if (input < 0) {i -= 1;}
	}
return i;
}

double min_distance(double array[], int length)
{
double min = array[0];
for (int i = 1; i < length; i ++)
{
	if (array[i] < min && array[i] != 0)
	{
		min = array[i];
	}
}
return min;
}


int main()
{
// file streams
ifstream inputfile;
ofstream hbonds_outputfile;

string infile;
int timesteps, nooa, noha;
double lattice;

cout << "Please enter the filename of your file: ";
cin >> infile;
cout << "Please enter the number of oxygen atoms: ";
cin >> nooa;
cout << "Please enter the number of hydrogen atoms: ";
cin >> noha;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
cout << "Program running...please wait a moment.\n\n";
// end of main menu*/

//read the inputfile
inputfile.open(infile.c_str());

int noa = nooa + noha;

double x[noa*timesteps], y[noa*timesteps], z[noa*timesteps];
int counter = 0;

while (!inputfile.eof())
{
	inputfile >> x[counter] >> y[counter] >> z[counter];	
	counter ++;

}


double oxyz[nooa][3], hxyz[noha][3], ohindices[nooa][4];

hbonds_outputfile.open("hbonds_histogram.dat");

for (int i = 0; i < timesteps; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		oxyz[j][0] = lattice*x[j + i*noa];
		oxyz[j][1] = lattice*y[j + i*noa];
		oxyz[j][2] = lattice*z[j + i*noa];
	}
	for (int j = 0; j < noha; j ++)
	{
		hxyz[j][0] = lattice*x[nooa + j + i*noa];
		hxyz[j][1] = lattice*y[nooa + j + i*noa];
		hxyz[j][2] = lattice*z[nooa + j + i*noa];
	}
	for (int j = 0; j < nooa; j++)
	{
		int count = 0;

		for (int k = 0; k < noha; k ++)
		{
			double dx = oxyz[j][0] - hxyz[k][0];
			double dy = oxyz[j][1] - hxyz[k][1];
			double dz = oxyz[j][2] - hxyz[k][2];
	
			dx -= lattice*pbc_round(dx/lattice);
			dy -= lattice*pbc_round(dy/lattice);
			dz -= lattice*pbc_round(dz/lattice);			
	
			double ohdist = sqrt ( dx*dx + dy*dy + dz*dz );
	
			if (ohdist < 1.1)
			{
				ohindices[j][count] = k;
				count ++;
			}
		}
	}

	for (int i = 0; i < nooa; i ++)
	{
		hbonds_outputfile << ohindices[i][0] << "\t" << ohindices[i][1] << "\t" << ohindices[i][2] << "\t" << ohindices[i][3] << endl;
	}	 
	
/*	for (int j = 0; j < nooa; j ++)
	{
		vector <int> ooindex, ohindex;

		for (int k = 0; k < nooa; k ++)
		{
			double odx = oxyz[k][0] - oxyz[j][0];
			double ody = oxyz[k][1] - oxyz[j][1];
			double odz = oxyz[k][2] - oxyz[j][2];

			double oodist = sqrt (odx*odx + ody*ody + odz*odz );

			if (oodist > 0.0 && oodist < 3.6)
			{
				ooindex.push_back(k);
			}
		}

		

		for (int k = 0; k < ooindex.size(); k ++)
		{
			double h1dx = 
	}*/
}

inputfile.close();
hbonds_outputfile.close();

return 0;
} 
		
		






















