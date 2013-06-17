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
ifstream inputfile;
ofstream rmsd_outputfile;

string infile;
int timesteps, number_of_atoms;
double lattice;

cout << "Please enter the filename of your file: ";
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

 rmsd_outputfile.open("msd.dat");

        int nooa = number_of_atoms/3;
        double oxyz[nooa*timesteps][3], distance[nooa*timesteps], rmsdistance[timesteps], sum[timesteps], COM[timesteps][3];

        for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < nooa; j ++)
                {
                        oxyz[j + i*nooa][0] = lattice*x[j + i*number_of_atoms];
                        oxyz[j + i*nooa][1] = lattice*y[j + i*number_of_atoms];
                        oxyz[j + i*nooa][2] = lattice*z[j + i*number_of_atoms];
                }
        }

        for(int i = 0; i < timesteps; i ++)
        {
                double dx(0), dy(0), dz(0);
                for (int j = 0; j < nooa; j ++)
                {
                        dx += oxyz[j + i*nooa][0];
                        dy += oxyz[j + i*nooa][1];
                        dz += oxyz[j + i*nooa][2];
                }

                dx /= nooa;
                dy /= nooa;
                dz /= nooa;

                COM[i][0] = dx;
                COM[i][1] = dy;
                COM[i][2] = dz;
        }

        for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < nooa; j ++)
                {
                        oxyz[j + i*nooa][0] -= COM[i][0];
                        oxyz[j + i*nooa][1] -= COM[i][1];
                        oxyz[j + i*nooa][2] -= COM[i][2];

                        oxyz[j + i*nooa][0] -= lattice*pbc_round(oxyz[j + i*nooa][0]/lattice);
                        oxyz[j + i*nooa][1] -= lattice*pbc_round(oxyz[j + i*nooa][1]/lattice);
                        oxyz[j + i*nooa][2] -= lattice*pbc_round(oxyz[j + i*nooa][2]/lattice);
                }
        }

        for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < nooa; j ++)
                {
                        double dxa = oxyz[j + i*nooa][0] - oxyz[0][0];
                        double dya = oxyz[j + i*nooa][1] - oxyz[0][1];
                        double dza = oxyz[j + i*nooa][2] - oxyz[0][2];

                        distance[j + i*nooa] = dxa*dxa + dya*dya + dza*dza ;
                }
        }
	 for (int i = 0; i < timesteps; i ++)
        {
                double add(0);
                for (int j = 0; j < nooa; j ++)
                {
                        add += distance[j + i*nooa];

                }
                sum[i] = add/nooa;
        }
        for (int i = 0; i < timesteps; i ++)
        {
                rmsd_outputfile <<  sum[i] << endl;
        }

cout << "\n\nYour rms displacement data has been placed in \"msd.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
rmsd_outputfile.close();
return 0;
}

