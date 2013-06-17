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
ofstream oh_outputfile;

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

oh_outputfile.open("oh_histogram.dat");

        double ox[number_of_atoms/3], oy[number_of_atoms/3], oz[number_of_atoms/3], hx[2*number_of_atoms/3], hy[2*number_of_atoms/3], hz[2*number_of_atoms/3];
        double dx1, dx2, dy1, dy2, dz1, dz2, distance[2*number_of_atoms/3], final_distance[2*number_of_atoms/3*timesteps], bin[200]={}  ;
        int bin_number;

        for (int i = 0; i < timesteps; i ++)
        {
                 for (int j = 0; j < number_of_atoms/3; j ++)
                {
                        ox[j] = lattice*x[j + i*number_of_atoms];
                        oy[j] = lattice*y[j + i*number_of_atoms];
                        oz[j] = lattice*z[j + i*number_of_atoms];
                }
                for (int k = 0; k < 2*number_of_atoms/3; k ++)
                {
                        hx[k] = lattice*x[number_of_atoms/3 + k + i*number_of_atoms];
                        hy[k] = lattice*y[number_of_atoms/3 + k + i*number_of_atoms];
                        hz[k] = lattice*z[number_of_atoms/3 + k + i*number_of_atoms];
                }
                for (int n = 0; n < number_of_atoms/3; n ++)
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

                        distance[2*n] = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
                        distance[2*n + 1] = sqrt ( dx2*dx2 + dy2*dy2 + dz2*dz2 );
                }
                for (int g = 0; g < 2*number_of_atoms/3; g ++)
                {
                        final_distance[g + i*2*number_of_atoms/3] = distance[g];
                }

        }
        double sum(0);
        int n(0);
        for (int i = 0; i < 2*number_of_atoms/3*timesteps; i ++)
        {

                bin_number = final_distance[i]*100;
                bin[bin_number] += 1;
                sum += final_distance[i];
                n ++;
        }
        oh_outputfile << "# The average is " << sum/n << "\n\n";

for (int j = 0; j < 200; j ++)
        {
                oh_outputfile << j/100. << "\t" << bin[j]/256000. << endl;
                oh_outputfile << (j + 1)/100. << "\t" << bin[j]/256000. << endl;
        }
        cout << "\n\nYour O-H distance histogram data has been placed in \"oh_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
oh_outputfile.close();
return 0;
}

