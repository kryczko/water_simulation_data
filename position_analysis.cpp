#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

int main()
{
// file streams
ifstream inputfile;
ofstream oo_outputfile, oh_outputfile, angle_outputfile;

char oodistance, ohdistance, angle; 
string infile;
int timesteps, number_of_atoms;

// main menu for the program
cout << "\n\n     Welcome to the bulk-water data analysis program!\n\n";
cout << "     If you would like to plot the average O-O distance with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> oodistance;
cout << "     If you would like to plot the average O-H distance with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> ohdistance;
cout << "     If you would like to plot the average H-O-H angle with respect to timesteps enter \"y\", if not enter \"n\": ";
cin >> angle;
cout << "     Please enter the filename of your file: ";
cin >> infile;
cout << "     Please enter the number of atoms: ";
cin >> number_of_atoms;
cout << "     Please enter the number of timesteps: ";
cin >> timesteps;
// end of main menu

//read the inputfile
inputfile.open(infile.c_str());

double x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!inputfile.eof())
{
	inputfile >> x[counter] >> y[counter] >> z[counter];
}
// done reading inputfile


//decision for oodistance
if (oodistance == 'y')
{
	oo_outputfile.open("oo_avg_distance.dat");


	cout << "\n\nYour O-O average distance data with respect to timesteps has been placed in \"oo_avg_distance.dat\" and can now be easily plotted with gnuplot.\n\n"; 
}
else
{
	cout << "\n\nYou either entered the wrong key or you do not want to plot the average O-O distance.\n\n";
}

//decision for ohdistance
if (ohdistance == 'y')
{
	oh_outputfile.open("oh_avg_distance.dat");


        cout << "\n\nYour O-H average distance data with respect to timesteps has been placed in \"oh_avg_distance.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the average O-H distance.\n\n";
}

//decision for angle
if (angle == 'y')
{
        angle_outputfile.open("avg_angle.dat");


        cout << "\n\nYour H-O-H average angle data with respect to timesteps has been placed in \"avg_angle.dat\" and can now be easily plotted with gnuplot.\n\n";
}
else
{
        cout << "\n\nYou either entered the wrong key or you do not want to plot the average H-O-H angle.\n\n";
}









return 0;
}
