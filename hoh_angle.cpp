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
ofstream angle_outputfile;

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

int number_of_atoms = nooa + noha;

double x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!inputfile.eof())
{
        inputfile >> x[counter] >> y[counter] >> z[counter];
        counter ++;

}


angle_outputfile.open("angle_histogram.dat");

double oxyz[nooa][3], hxyz[noha][3];

vector <double> angles;

for (int i = 0; i < timesteps; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		oxyz[j][0] = lattice*x[j + i*number_of_atoms];
		oxyz[j][1] = lattice*y[j + i*number_of_atoms];
		oxyz[j][2] = lattice*z[j + i*number_of_atoms];
	}
	
	for (int j = 0; j < noha; j ++)
	{
		hxyz[j][0] = lattice*x[nooa + j + i*number_of_atoms];
		hxyz[j][1] = lattice*y[nooa + j + i*number_of_atoms];
		hxyz[j][2] = lattice*z[nooa + j + i*number_of_atoms];
	}
		
	for (int j = 0; j < nooa; j ++)
	{
		int count = 0;
		vector <int> ohindex;
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
				ohindex.push_back(k);
				count ++;
			}
		}
			if (count == 2)
			{
				double dx1 = oxyz[j][0] - hxyz[ohindex[0]][0];
				double dx2 = oxyz[j][0] - hxyz[ohindex[1]][0];
				double dy1 = oxyz[j][1] - hxyz[ohindex[0]][1];
                                double dy2 = oxyz[j][1] - hxyz[ohindex[1]][1];
				double dz1 = oxyz[j][2] - hxyz[ohindex[0]][2];
                                double dz2 = oxyz[j][2] - hxyz[ohindex[1]][2];

				dx1 -= lattice*pbc_round(dx1/lattice);
                        	dy1 -= lattice*pbc_round(dy1/lattice);
                        	dz1 -= lattice*pbc_round(dz1/lattice);
				dx2 -= lattice*pbc_round(dx2/lattice);
                        	dy2 -= lattice*pbc_round(dy2/lattice);
                        	dz2 -= lattice*pbc_round(dz2/lattice);
	
				double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
				double dist1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
				double dist2 = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2 ); 
				
				double angle = acos ( dot / (dist1*dist2) ) * 57.2957795;
				angles.push_back(angle);
			}
		
	}		
}
double bin[2000] = {};
double sum = 0;
for (int i = 0; i < angles.size(); i ++)
{
	sum += angles[i];
	int bin_num = angles[i]*10.;
	bin[bin_num] ++;
}
cout << "The average is " << sum/(nooa*timesteps) << " degrees." << endl;
for (int i = 0; i < 2000; i ++)
{
	angle_outputfile << i/10. << "\t" << bin[i]/(nooa*timesteps) << endl;
	angle_outputfile << (i + 1)/10. << "\t" << bin[i]/(nooa*timesteps) << endl;
}


       cout << "\n\nYour H-O-H angle histogram data has been placed in \"angle_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
angle_outputfile.close();

return 0;
}






























