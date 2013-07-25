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
ofstream oo_outputfile;

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

//read the inputfile
inputfile.open(infile.c_str());

double x[(nooa+noha)*timesteps], y[(nooa+noha)*timesteps], z[(nooa+noha)*timesteps];
int counter = 0;

while (!inputfile.eof())
{
        inputfile >> x[counter] >> y[counter] >> z[counter];
        counter ++;

}

// done reading inputfile

 oo_outputfile.open("oo_histogram.dat");
        double ox[nooa], oy[nooa], oz[nooa];
        double dx, dy, dz, lowest;
        int bin_number;
        double last_difference[(nooa-1)*timesteps], bin[100] = {};

        for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < nooa; j ++)
                {
                        ox[j] = lattice*x[j + i*(nooa+noha)];
                        oy[j] = lattice*y[j + i*(nooa+noha)];
                        oz[j] = lattice*z[j + i*(nooa+noha)];
                }

                for (int k = 0; k < nooa - 1; k ++)
                {
                        double difference[nooa - (k+1)];
                        for (int n = k + 1; n < nooa; n ++)
                        {
                                dx = ox[n] - ox[k];
                                dy = oy[n] - oy[k];
                                dz = oz[n] - oz[k];

                                dx -= lattice*pbc_round(dx/lattice);
                                dy -= lattice*pbc_round(dy/lattice);
                                dz -= lattice*pbc_round(dz/lattice);

                                double distance = sqrt( dx*dx + dy*dy + dz*dz );
                                difference[n - (k + 1)] = distance;

                        }

                last_difference[k + i*nooa] = min_distance(difference, nooa - (k + 1));
                }



        }
        int n = 0;
        double sum(0);
        for ( int i = 0; i < (nooa-1)*timesteps; i ++)
        {
                if ((last_difference[i] != 0 || last_difference[i] != 0.0) && (last_difference[i] < 5.0))
                {
                        bin_number = last_difference[i]*10;
                        bin[bin_number] += 1;
                        sum += last_difference[i];
                        n ++;
                }
        }
        oo_outputfile << "# The average is " << sum/n << "\n\n";
        for(int i = 0; i < 100; i ++)
        {
                oo_outputfile << i/10. << "\t" << (bin[i]/((nooa/2)*timesteps))*100 << endl;
                oo_outputfile << (i + 1)/10. << "\t" << (bin[i]/((nooa/2)*timesteps))*100 << endl;
        }
        cout << "\n\nYour O-O histogram data has been placed in \"oo_histogram.dat\" and can now be easily plotted with gnuplot.\n\n";

inputfile.close();
oo_outputfile.close();
return 0;
}


