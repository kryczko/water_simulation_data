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
ofstream oh_outputfile;

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

inputfile.open(infile.c_str());

double x[(nooa + noha)*timesteps], y[(nooa + noha)*timesteps], z[(nooa + noha)*timesteps];
int counter = 0;

while (!inputfile.eof())
{
        inputfile >> x[counter] >> y[counter] >> z[counter];
        counter ++;

}

oh_outputfile.open("oh_histogram.dat");

double oxyz[nooa][3], hxyz[noha][3], bin[200] = {};
int ohindices[nooa][4]; 
vector <double> ohdistance;

for (int i = 0; i < nooa; i ++)
{
	for (int y = 0; y < 4; y ++)
	{
		ohindices[i][y] = -1;
	}
}


for (int i = 0; i < timesteps; i ++)
{
	for (int j = 0; j < nooa; j ++)
        {
	        oxyz[j][0] = lattice*x[j + i*(nooa + noha)];
                oxyz[j][1] = lattice*y[j + i*(nooa + noha)];
                oxyz[j][2] = lattice*z[j + i*(nooa + noha)];
       	}

	for (int j = 0; j < noha; j ++)
	{
		hxyz[j][0] = lattice*x[nooa + j + i*(nooa + noha)];
		hxyz[j][1] = lattice*y[nooa + j + i*(nooa + noha)];
		hxyz[j][2] = lattice*z[nooa + j + i*(nooa + noha)];
	}

	for (int j = 0; j < nooa; j ++)
	{
		for (int k = 0; k < noha; k ++)
		{
			double dx = oxyz[j][0] - hxyz[k][0];
			double dy = oxyz[j][1] - hxyz[k][1];
			double dz = oxyz[j][2] - hxyz[k][2];
		
			dx -= lattice*pbc_round(dx/lattice);
			dy -= lattice*pbc_round(dy/lattice);
			dz -= lattice*pbc_round(dz/lattice);

			double distance = sqrt(dx*dx + dy*dy + dz*dz); 
			
			if (distance <= 1.1)
			{	
				ohdistance.push_back (distance);
			}
		}
	}
}

for (int i = 0; i < ohdistance.size(); i ++)
{
	int bin_num = ohdistance[i]*100;
	bin[bin_num] ++;
}

for (int i = 0; i < 200; i ++)
{
	oh_outputfile << i/100. << "\t" << bin[i]/(timesteps*noha) << endl;
	oh_outputfile << (i+1)/100. << "\t" << bin[i]/(timesteps*noha) << endl;
}

cout << "\n\nYour O-H distance histogram data has been placed in \"oh_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
oh_outputfile.close();
return 0;
}


















