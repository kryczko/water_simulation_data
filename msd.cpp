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

double xyz[number_of_atoms*timesteps][3];
int counter = 0;

while (!inputfile.eof())
{
        inputfile >> xyz[counter][0] >> xyz[counter][1] >> xyz[counter][2];
        counter ++;

}

 rmsd_outputfile.open("msd.dat");

        int noa = number_of_atoms;
        double distance[noa*(timesteps - 1)], rmsdistance[timesteps - 1], sum[timesteps - 1], COM[timesteps][3];

        for(int i = 0; i < timesteps; i ++)
        {
                double dx(0), dy(0), dz(0);
                for (int j = 0; j < noa; j ++)
                {
                        dx += xyz[j + i*noa][0];
                        dy += xyz[j + i*noa][1];
                        dz += xyz[j + i*noa][2];
                }

                dx /= noa;
                dy /= noa;
                dz /= noa;

                COM[i][0] = dx;
                COM[i][1] = dy;
                COM[i][2] = dz;
        }

        for (int i = 1; i < timesteps; i ++)
        {
                for (int j = 0; j < noa; j ++)
                {
                        double dx1 = xyz[j][0] - COM[0][0];
			double dy1 = xyz[j][1] - COM[0][1];
			double dz1 = xyz[j][2] - COM[0][2];
			double dx2 = xyz[j + i*noa][0] - COM[i][0];
			double dy2 = xyz[j + i*noa][1] - COM[i][1];
			double dz2 = xyz[j + i*noa][2] - COM[i][2];			
			
			double dxa = dx2 - dx1;
			double dya = dy2 - dy1;
			double dza = dz2 - dz1;

                        distance[j + (i-1)*noa] = dxa*dxa + dya*dya + dza*dza ;
                }
        }
	 for (int i = 0; i < timesteps - 1; i ++)
        {
                double add(0);
                for (int j = 0; j < noa; j ++)
                {
                        add += distance[j + i*noa];

                }
                sum[i] = add/noa;
        }
        for (int i = 0; i < timesteps - 1; i ++)
        {
                rmsd_outputfile <<  sum[i] << endl;
        }

cout << "\n\nYour rms displacement data has been placed in \"msd.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
rmsd_outputfile.close();
return 0;
}

