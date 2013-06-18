#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

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

	double ox[nooa], oy[nooa], oz[nooa], hx[noha], hy[noha], hz[noha];
	double dx1, dx2, dy1, dy2, dz1, dz2, dot[nooa], distance[noha], angle[nooa*timesteps], bin[20000] = {};
	int bin_number;

	for (int i = 0; i < timesteps; i ++)
	{
		 for (int j = 0; j < nooa; j ++)
                {
                        ox[j] = lattice*x[j + i*number_of_atoms];
                        oy[j] = lattice*y[j + i*number_of_atoms];
                        oz[j] = lattice*z[j + i*number_of_atoms];
                }
                for (int k = 0; k < noha; k ++)
                {
                        hx[k] = lattice*x[nooa + k + i*number_of_atoms];
                        hy[k] = lattice*y[nooa + k + i*number_of_atoms];
                        hz[k] = lattice*z[nooa + k + i*number_of_atoms];
                }
                for (int n = 0; n < nooa; n ++)
                {
                        dx1 = ox[n] - hx[2*n];
                        dx2 = ox[n] - hx[2*n + 1];
                        dy1 = oy[n] - hy[2*n];
                        dy2 = oy[n] - hy[2*n + 1];
                        dz1 = oz[n] - hz[2*n];
                        dz2 = oz[n] - hz[2*n + 1];

                        dx1 -= lattice*pbc_round(dx1/lattice);
                        dy1 -= lattice*pbc_round(dy1/lattice);
                        dz1 -= lattice*pbc_round(dz1/lattice);

                        dx2 -= lattice*pbc_round(dx2/lattice);
                        dy2 -= lattice*pbc_round(dy2/lattice);
                        dz2 -= lattice*pbc_round(dz2/lattice);

			dot[n] = dx1*dx2 + dy1*dy2 + dz1*dz2;

                        distance[2*n] = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
                        distance[2*n + 1] = sqrt ( dx2*dx2 + dy2*dy2 + dz2*dz2 );
                }
                for (int g = 0; g < nooa; g ++)
                {
                        angle[g + i*nooa] = acos( dot[g]/(distance[2*g]*distance[2*g + 1]) ) * 57.2957795;
                }

	}
	double sum(0);
	int n(0);
	for (int i = 0; i < nooa*timesteps; i ++)
	{
		bin_number = angle[i]*100;
		bin[bin_number] += 1;
		sum += angle[i];
		n ++;
	}
	angle_outputfile << "# The average is " << sum/n << "\n\n";
	for (int j = 0; j < 20000; j ++)
	{
		angle_outputfile << j/100. << "\t" << bin[j]/128000. << endl;
		angle_outputfile << (j + 1)/100. << "\t" << bin[j]/128000. << endl;
	}

        cout << "\n\nYour H-O-H angle histogram data has been placed in \"angle_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
angle_outputfile.close();

return 0;
}
